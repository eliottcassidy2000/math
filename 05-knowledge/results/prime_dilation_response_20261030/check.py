"""Joint prime-dilation responses: exact lattice, clocks, and hostile probes."""
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
GATES = 0


def need(value, why):
    global GATES
    GATES += 1
    if not value:
        raise RuntimeError(why)


@lru_cache(None)
def factors(n):
    out, p = {}, 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


@lru_cache(None)
def sigma(n):
    return prod((p**(a+1)-1)//(p-1) for p,a in factors(n).items())


def aliquot(n):
    return sigma(n)-n


def vp(n, p):
    if n == 0:
        raise ValueError('finite valuation requested at zero')
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def source(n, primes):
    core, exponents = n, []
    for p in primes:
        a = 0
        while core % p == 0:
            core //= p
            a += 1
        exponents.append(a)
    return core, exponents, sigma(core)


def response(n, primes, t):
    core, exponents, c = source(n, primes)
    return c*prod((p**(a+b+1)-1)//(p-1)
                  for p,a,b in zip(primes,exponents,t))-n*prod(p**b for p,b in zip(primes,t))


def cube(primes):
    return [prod(p for i,p in enumerate(primes) if mask >> i & 1)
            for mask in range(1 << len(primes))]


def hankel(n, primes):
    q = cube(primes)
    return [[aliquot(n*a*b) for b in q] for a in q]


def determinant(matrix):
    a = [row[:] for row in matrix]
    d, sign, previous = len(a), 1, 1
    if d == 0:
        return 1
    for k in range(d-1):
        pivot = next((j for j in range(k,d) if a[j][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[pivot], a[k] = a[k], a[pivot]
            sign *= -1
        for i in range(k+1,d):
            for j in range(k+1,d):
                value = a[i][j]*a[k][k]-a[i][k]*a[k][j]
                if value % previous:
                    raise RuntimeError('nonintegral Bareiss division')
                a[i][j] = value//previous
            a[i][k] = 0
        previous = a[k][k]
    return sign*a[-1][-1]


def direct_smith(matrix):
    """Independent determinantal divisors, exhaustive minors, no formula input."""
    d, previous, out = len(matrix), 1, []
    for k in range(1,d+1):
        divisor = 0
        subsets = list(combinations(range(d),k))
        for rows in subsets:
            for cols in subsets:
                value = abs(determinant([[matrix[i][j] for j in cols] for i in rows]))
                divisor = gcd(divisor,value)
                if divisor == 1:
                    break
            if divisor == 1:
                break
        need(divisor > 0 and divisor % previous == 0, 'determinantal divisibility')
        out.append(divisor//previous)
        previous = divisor
    return out


def predicted_smith(n, primes):
    r, d = len(primes), 1 << len(primes)
    core, exponents, c = source(n, primes)
    x = prod(p**(a+1) for p,a in zip(primes,exponents))
    gamma = c*x-n*prod(p-1 for p in primes)
    need(gamma > 0, 'positive top coefficient')
    if r == 1:
        first = gcd(n,c)
        return [first,c*gamma//first]
    support = set(primes) | set(factors(c)) | set(factors(gamma)) | set(factors(n))
    entries = [1]*d
    for ell in support:
        cv, nv, gv = vp(c,ell), vp(n,ell), vp(gamma,ell)
        if ell in primes:
            valuations = [cv]*(d//2)+[nv]+[cv+nv+1]*(d//2-2)+[cv+1+gv]
        else:
            low = min(cv,nv)
            valuations = [low]+[cv]*(d-2)+[cv+gv-low]
        for i,a in enumerate(sorted(valuations)):
            entries[i] *= ell**a
    need(all(b % a == 0 for a,b in zip(entries,entries[1:])), 'predicted Smith chain')
    return entries


def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def tensor_local(local):
    d = 1 << len(local)
    return [[prod(local[i][(a>>i)&1][(b>>i)&1] for i in range(len(local)))
             for b in range(d)] for a in range(d)]


def check_arrow(n, primes):
    core, exponents, c = source(n,primes)
    r, d = len(primes), 1 << len(primes)
    u, v = [], []
    for p,a in zip(primes,exponents):
        power = p**(a+1)
        rep = (power-1)//(p-1)
        u.append([[-p,1],[p*rep+1,-rep]])
        v.append([[1,-1],[0,1]])
    left, right = tensor_local(u), tensor_local(v)
    actual = matmul(matmul(left,hankel(n,primes)),right)
    diagonal = [c*prod(p**(a+1) for i,(p,a) in enumerate(zip(primes,exponents))
                       if mask>>i&1) for mask in range(d)]
    wanted = [[diagonal[i] if i == j else 0 for j in range(d)] for i in range(d)]
    for mask in range(d):
        wanted[-1][mask] -= n*prod(p-1 for i,p in enumerate(primes) if mask>>i&1)
    need(actual == wanted, 'literal unimodular arrowhead transform')
    need(abs(determinant(left)) == abs(determinant(right)) == 1, 'integer basis unimodular')
    gamma = wanted[-1][-1]
    expected = (-1)**(r*2**(r-1))*c**(d-1)*gamma*prod(
        p**((a+1)*(2**(r-1)-1)) for p,a in zip(primes,exponents))
    need(determinant(hankel(n,primes)) == expected, 'full Hankel determinant')
    need(prod(predicted_smith(n,primes)) == abs(expected), 'Smith determinant product')
    return expected


def matrix_checks():
    counts = Counter()
    for p in (2,3,5,7,11):
        for n in range(1,81):
            primes = (p,)
            check_arrow(n,primes)
            need(direct_smith(hankel(n,primes)) == predicted_smith(n,primes), 'one-port Smith')
            counts['r1'] += 1
    for primes in ((2,3),(3,7),(7,11),(7,13),(11,13),(5,19)):
        for n in range(1,81):
            check_arrow(n,primes)
            need(direct_smith(hankel(n,primes)) == predicted_smith(n,primes), 'two-port Smith')
            counts['r2'] += 1
    examples = []
    for primes in ((2,3,5),(3,7,11),(7,11,13)):
        for n in (1,57,343,578,1001):
            det = check_arrow(n,primes)
            actual = direct_smith(hankel(n,primes))
            need(actual == predicted_smith(n,primes), 'three-port exhaustive-minor Smith')
            examples.append({'n':n,'primes':primes,'smith':actual,'determinant':det})
            counts['r3'] += 1
    for primes in ((2,3,5,7),(7,11,13,19)):
        for n in (1,57,343,1001):
            check_arrow(n,primes)
            counts['r4_arrow_det_only'] += 1
    return {'universes':dict(counts),'three_port_examples':examples,
            'small_example':{'primes':[2,3],'n':1,'H':hankel(1,(2,3)),
                             'smith':predicted_smith(1,(2,3))}}


def shift(state, primes, i, modulus):
    p, bit = primes[i], 1 << i
    return tuple(((p+1)*state[j]-p*state[j^bit] if j&bit else state[j|bit]) % modulus
                 for j in range(len(state)))


def profile(n, primes, t, modulus):
    return tuple(response(n,primes,[b+((mask>>i)&1) for i,b in enumerate(t)]) % modulus
                 for mask in range(1 << len(primes)))


def clock_order(p, ell, k):
    need((p-1) % ell == 0, 'unipotent prime guard')
    if ell != 2 or p % 4 == 1:
        return ell**k
    return 2**max(1,k-vp(p+1,2)+1)


def orbit_states(n, primes, modulus):
    initial = tuple(aliquot(n*d) % modulus for d in cube(primes))
    seen, queue = {initial}, deque([initial])
    while queue:
        state = queue.popleft()
        for i in range(len(primes)):
            following = shift(state,primes,i,modulus)
            if following not in seen:
                seen.add(following)
                queue.append(following)
    return seen


def clock_checks():
    rows = []
    cases = [(2,(7,11,13),k) for k in range(1,6)]
    cases += [(3,(7,13),k) for k in range(1,4)]
    cases += [(3,(7,13,19),2),(5,(11,31),2)]
    for ell,primes,k in cases:
        modulus = ell**k
        periods = [clock_order(p,ell,k) for p in primes]
        states = orbit_states(1,primes,modulus)
        direct = {profile(1,primes,t,modulus) for t in product(*(range(x) for x in periods))}
        need(states == direct, 'literal arithmetic grid equals closed response orbit')
        need(len(states) == prod(periods), 'minimal deterministic unipotent clock carrier')
        for i,period in enumerate(periods):
            state = tuple(aliquot(d) % modulus for d in cube(primes))
            initial = state
            for t in range(1,period+1):
                state = shift(state,primes,i,modulus)
                need((state == initial) == (t == period), 'exact local clock period')
        rows.append({'ell':ell,'k':k,'primes':primes,'periods':periods,'states':len(states),
                     'linear_module_size':str(modulus**(1 << len(primes)))})
    dark = []
    for k in range(1,6):
        modulus = 2**k
        for n in (343,1001,57,6):
            actual = orbit_states(n,(7,11,13),modulus)
            smith = predicted_smith(n,(7,11,13))
            module_size = prod(modulus//gcd(modulus,a) for a in smith)
            dark.append({'n':n,'k':k,'states':len(actual),'linear_module_size':module_size})
    return {'unit_core_cases':rows,'source_comparison':dark}


def compiler_checks():
    probes = 0
    for n in range(1,301):
        # sigma has an independent divisor-pair control at these sources.
        direct = sum(d for d in range(1,n+1) if n % d == 0)
        need(direct == sigma(n), 'direct divisor sigma')
        for p in (2,3,7):
            x,y = aliquot(n),aliquot(p*n)
            a = vp(y-x,p)
            recovered = ((p**(a+1)-1)//(p-1))*(y-p*x)-x
            need(recovered == n, 'inherited prime-dilation reconstruction')
            for multiplier in (1,2,3,6,10,21,49,1001):
                total = x+recovered
                for q,b in factors(multiplier).items():
                    e = vp(recovered,q)
                    numerator,denominator = q**(e+b+1)-1,q**(e+1)-1
                    need(total*numerator % denominator == 0, 'exact exponent-update division')
                    total = total*numerator//denominator
                need(total-multiplier*n == aliquot(multiplier*n), 'all-prime compiled response')
                probes += 1
    carry_cases = 0
    for primes in ((2,3),(7,11,13),(3,5,7,11)):
        for n in (1,6,25,57,343,578,1001):
            core,exponents,c = source(n,primes)
            observed = sum(prod(-p for i,p in enumerate(primes) if not (mask>>i&1))*aliquot(n*d)
                           for mask,d in enumerate(cube(primes)))
            need(observed == c, 'full mixed-carry extraction')
            for t in product(range(3),repeat=len(primes)):
                need(response(n,primes,t) == aliquot(n*prod(p**b for p,b in zip(primes,t))),
                     'factorized response matches actual arithmetic')
            carry_cases += 1
    return {'compiled_queries':probes,'full_mixed_carries':carry_cases}


def hostile_checks():
    all_primes = (7,11,13,17,19,23)
    families = []
    for r in range(1,7):
        primes = all_primes[:r]
        n, outside = prod(primes), 3
        a = [aliquot(n*d) % 2 for d in cube(primes)]
        b = [aliquot(outside*d) % 2 for d in cube(primes)]
        need(a[:-1] == b[:-1] == [1]*((1<<r)-1), 'every proper-subset response agrees')
        need(a[-1] == 0 and b[-1] == 1, 'full joint response separates actual sources')
        families.append({'primes':primes,'source':n,'outside_source':outside,
                         'agreeing_proper_probes':(1<<r)-1})
    for primes in ((3,7),(7,11,13),(3,5,7,11)):
        for t in product(range(4),repeat=len(primes)):
            need(response(25,primes,t) % 2 == response(49,primes,t) % 2,
                 'full-dilation trace equivalence control')
    need(aliquot(aliquot(25)) % 2 != aliquot(aliquot(49)) % 2,
         'actual aliquot step breaks dilation-only carrier')
    chain = []
    for n in (1001,343,57):
        matrix = hankel(n,(7,11,13))
        smith = predicted_smith(n,(7,11,13))
        rank2 = sum(a % 2 != 0 for a in smith)
        chain.append({'n':n,'aliquot':aliquot(n),'binary_Hankel_rank':rank2,
                      'binary_states':len(orbit_states(n,(7,11,13),2))})
    need([x['binary_Hankel_rank'] for x in chain] == [8,8,1], '1001 chain rank change')
    return {'all_order_joint_hostile_controls':families,'chain':chain,
            'new_operation_hostile':{'sources':[25,49],'aliquot_sources':[6,8],
                                     'next_observation_parities':[0,1]}}


def replay_tournament_control():
    path = ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec = importlib.util.spec_from_file_location('tournament_control',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.hostile_checks()


def main():
    results = {'matrices':matrix_checks(),'clocks':clock_checks(),'compiler':compiler_checks(),
               'hostiles':hostile_checks(),'hostile_81_27':replay_tournament_control()}
    results['local_gates'] = GATES
    results['inherited_tournament_gates'] = 7
    print(json.dumps(results,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
