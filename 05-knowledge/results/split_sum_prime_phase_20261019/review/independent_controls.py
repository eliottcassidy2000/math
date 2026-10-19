#!/usr/bin/env python3
"""Independent finite controls for both prime packets (standard library only).

No draft imports, floating-point characters, or assert statements.  Uses literal
binary addresses, trial division, signed multiplicative orders, and literal
split-factorization rather than the supplied recurrence-based test functions.
All write targets are in this review directory.
"""
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement
from math import isqrt, prod
from pathlib import Path
import hashlib
import json

CHECKS = 0


def need(test, context):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(context)


def factor(n):
    if type(n) is not int or n < 1:
        raise ValueError(n)
    out = []
    d = 2
    while d*d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            out.append([d, e])
        d = 3 if d == 2 else d+2
    if n > 1:
        out.append([n, 1])
    return out


def is_prime(n):
    return n >= 2 and factor(n) == [[n, 1]]


def divisors(ps):
    # Literal subsets, independent of the supplied divisor-doubling helper.
    return [prod(c) for size in range(len(ps)+1) for c in combinations(ps, size)]


def literal_kernel(ts, h):
    c = Counter(sum(ts[i] for i in range(len(ts)) if mask >> i & 1) % h
                for mask in range(1 << len(ts)))
    return tuple(c[x] for x in range(h))


def signed_prime_data(ps, q):
    n = prod(ps)
    ds = divisors(ps)
    counts = Counter(min(d % q, q-d % q) for d in ds)
    counts = tuple(counts[x] for x in range(1, (q+1)//2))
    hits = sorted(d for d in ds if (d+n//d) % q == 0)
    squares = {x*x % q for x in range(1, q)}
    phase = 1 if (-n) % q in squares else -1
    return n, counts, phase, hits


def audit_cyclic():
    cases = 0
    by_kernel = defaultdict(int)
    stats = Counter()
    for h in range(1, 25):
        for k in range(5):
            for ts in combinations_with_replacement(range(h), k):
                c = literal_kernel(ts, h)
                total = sum(ts) % h
                periods = {a for a in range(h) if c[a:]+c[:a] == c}
                # The reflection x -> a-x has array [c[a], c[a-1], ...].
                rev = c[:1] + c[:0:-1]
                centers = {a for a in range(h) if rev[-a:]+rev[:-a] == c}
                need(centers == {(total+a) % h for a in periods}, ('centers',h,ts))
                need((1 << k) % len(periods) == 0, ('orbit divisibility',h,ts))
                need(len(periods) & (len(periods)-1) == 0, ('2-group periods',h,ts))
                uniform = len(set(c)) == 1
                power2 = h & (h-1) == 0
                # Layer membership by gcd, independently of discrete logs/v_2.
                required = {1 << j for j in range(h.bit_length()-1)}
                observed = set()
                for t in ts:
                    if t:
                        g = 1
                        while t % (2*g) == 0:
                            g *= 2
                        observed.add(g)
                need(uniform == (power2 and required <= observed), ('balance',h,ts))
                if h % 2 == 0:
                    sylow = 1
                    while h % (2*sylow) == 0:
                        sylow *= 2
                    need((len({a % 2 for a in centers}) == 1) == (len(periods) < sylow),
                         ('parity criterion',h,ts))
                    by_kernel[(h,k,c)] |= 1 << (total % 2)
                if power2 and h > 1:
                    if not uniform:
                        need({a % 2 for a in centers} == {total % 2}, ('phase recovery',h,ts))
                    elif c[0] == 1:
                        need(total % 2 == 1, ('minimal blocked phase',h,ts))
                cases += 1
                stats['uniform' if uniform else 'nonuniform'] += 1
    ambiguous = []
    for (h,k,c), phases in by_kernel.items():
        if h & (h-1) == 0:
            expected = 3 if len(set(c)) == 1 and c[0] % 2 == 0 else (2 if len(set(c)) == 1 else phases)
            need(phases == expected, ('same-length ambiguity iff',h,k,c,phases))
            if len(set(c)) != 1:
                need(phases in (1,2), ('nonuniform same-kernel phase uniqueness',h,k,c))
        if phases == 3:
            ambiguous.append([h,k,list(c)])
    return {'all_multisets_h_1_to_24_k_0_to_4': cases,
            'counts': dict(stats), 'opposite_phase_kernel_classes':len(ambiguous),
            'sample_mixed_order_ambiguities': [x for x in ambiguous if x[0] & (x[0]-1)][:8]}


def audit_arithmetic():
    q = 65537
    need(is_prime(q), '65537 primality by trial division')
    need(pow(3,32768,q) == q-1, '3 has full 2-power multiplicative order')
    primes = [p for p in range(2,q) if is_prime(p)]
    layers = Counter()
    firsts = {}
    for p in primes:
        # Order of [p] in F_q^*/{+1,-1}, found by squaring instead of logarithms.
        x, exponent = p, 0
        while x not in (1,q-1):
            x = x*x % q
            exponent += 1
        layer = 15-exponent if exponent else 'zero'
        layers[layer] += 1
        firsts.setdefault(layer,p)
    expected = [3260,1631,817,444,203,93,49,26,9,4,2,1,1,2,0]
    need([layers[i] for i in range(15)] == expected, 'independent signed-order layer census')
    roots_minus_one = [x for x in range(1,q) if x*x % q == q-1]
    need(roots_minus_one == [256,65281], 'unique layer-14 full residues')
    repair = 327941
    candidates = sorted({r+j*q for r in roots_minus_one for j in range(repair//q+1)
                         if r+j*q <= repair})
    candidate_factors = [{'n':n,'factorization':factor(n)} for n in candidates]
    need([n for n in candidates if is_prime(n)] == [repair], 'minimum repair prime: every candidate exhausted')
    need(len(primes) == 6542 and len(primes)+1-15 == 6528, 'primorial balance exponent')
    basis = [3,19,17,13,233,839,431,157,6373,10163,257,2,64513,61441,327941]
    need(all(is_prime(p) for p in basis+[5]), 'all source factors prime')
    need(len(set(basis+[5])) == 16, 'squarefree source')
    n0,c0,phase0,hits0 = signed_prime_data(basis,q)
    n1,c1,phase1,hits1 = signed_prime_data(basis+[5],q)
    need(set(c0)=={1} and phase0 == -1 and not hits0, 'literal minimal balanced blocked source')
    need(set(c1)=={2} and phase1 == 1 and len(hits1)==2, 'literal extremal generating source')
    need(n1 == 36054720350329695914767176044436363972270, 'source integer')
    need(hits1 == [135625235,265840795411928288380603919649882], 'two literal witnesses')
    need(hits1[0]*hits1[1] == n1, 'witness complement pairing')
    quotient = sum(hits1)//q
    need(sum(hits1) == q*4056346726458768152045471341, 'exact split quotient')
    return {'method':'Trial division and signed multiplicative orders; no sieve or logarithm table',
            'prime_count_below_65537':len(primes),'layer_counts':expected,
            'minimum_repair_exhaustive_candidates':candidate_factors,
            'minimal_basis':basis,'generating_source':str(n1),
            'split_divisors':[str(d) for d in hits1],'split_quotient':str(quotient),
            'full_primorial_plus_repair_uniform_multiplicity':'2^6528'}


def audit_hostiles_and_sources():
    hostiles=[]
    for q,a,b in [(5,[2,3,7],[2,3,11]),(7,[2,3,5,11],[2,3,5,17]),
                  (5,[11],[19]),(17,[2,43],[2,59])]:
        na,ka,pa,ha=signed_prime_data(a,q)
        nb,kb,pb,hb=signed_prime_data(b,q)
        need(ka==kb,('hostile identical kernels',q,a,b))
        need(not ha and len(hb)>0,('hostile different candidate outcomes',q,a,b))
        if a in ([11],[2,43]):
            need(pa==pb==1,('root hostile same positive phase',q))
        else:
            need(pa==-1 and pb==1 and min(ka)>0,('phase hostile full support',q))
        hostiles.append({'q':q,'sources':[na,nb],'kernel':ka,'phases':[pa,pb],
                         'ordered_hit_counts':[len(ha),len(hb)]})
    n,c,p,h=signed_prime_data([2],11)
    need(p==1 and not h, 'positive phase without support does not suffice')
    # Actual complete factorization of every unordered split: no residue support theorem.
    expected=[2,3,5,11,37,13,7,29,17,19]
    selected=[]
    rows=[]
    for target in expected:
        n=prod(selected)
        ds=divisors(selected)
        sums=sorted({d+n//d for d in ds})
        factors=[factor(v) for v in sums]
        chosen=min(p for row in factors for p,e in row)
        need(chosen==target,('literal split-factor greedy step',selected,chosen,target))
        rows.append({'selected':chosen,'source':n,'distinct_unordered_split_sums':len(sums)})
        selected.append(chosen)
    # A nonresidue-multiplier countermodel for the original N+1 process is not a trajectory.
    cumulative=[]
    n=1
    for r in [3,5]*8:
        need(pow(r,3,7)==6,'original-Mullin hostile multiplier is nonresidue')
        n=n*r%7
        cumulative.append(n)
    need(set(cumulative)=={1,3},'nonresidue tail does not force product -1')
    return {'arithmetic_hostiles':hostiles,'literal_greedy_factorization_steps':rows,
            'original_Mullin_residue_countermodel_products':sorted(set(cumulative))}


def main():
    out={'status':'PASS','scope':'Independent bounded controls, not a prime completeness proof',
         'cyclic':audit_cyclic(),'arithmetic':audit_arithmetic(),
         'hostiles_and_sources':audit_hostiles_and_sources()}
    out['always_active_checks']=CHECKS
    payload=json.dumps(out,sort_keys=True,indent=2)+'\n'
    path=Path(__file__).with_name('independent_receipt.json')
    path.write_text(payload,encoding='utf-8')
    print('PASS:', CHECKS, 'always-active checks')
    print('receipt sha256:',hashlib.sha256(payload.encode()).hexdigest())
    print('cyclic multisets:',out['cyclic']['all_multisets_h_1_to_24_k_0_to_4'])
    print('minimum repair candidates:',json.dumps(out['arithmetic']['minimum_repair_exhaustive_candidates']))


if __name__=='__main__':
    main()
