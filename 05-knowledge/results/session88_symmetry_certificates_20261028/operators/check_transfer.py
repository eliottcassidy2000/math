#!/usr/bin/env python3
"""Stdlib exact checks for TRANSFER.md. No input files or stochastic sampling."""
from fractions import Fraction as F
from itertools import combinations, product
from functools import reduce
from math import gcd, lcm, isqrt
import json


def require(value, detail):
    if not value:
        raise RuntimeError(detail)


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def ident(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def minus(a, b):
    return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]


def scale(c, a):
    return [[c*x for x in r] for r in a]


def partitions(n):
    if not n:
        yield ()
        return
    def rec(a):
        if len(a) == n:
            yield tuple(a)
        else:
            for x in range(max(a)+2):
                yield from rec(a+[x])
    yield from rec([0])


def join(a, b):
    parent = list(range(len(a)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(len(a)):
        for j in range(i):
            if a[i] == a[j] or b[i] == b[j]:
                parent[root(i)] = root(j)
    return tuple(root(i) for i in range(len(a)))


def operator(a, w, normalized=True):
    mass = {x:sum(v for t,v in zip(a,w) if t == x) for x in a}
    return [[(F(w[j],mass[a[i]]) if normalized else w[j]) if a[i] == a[j] else 0
             for j in range(len(a))] for i in range(len(a))]


def factors(a, b, w):
    c = join(a,b)
    for z in set(c):
        inds = [i for i,x in enumerate(c) if x == z]
        total = sum(w[i] for i in inds)
        for x in {a[i] for i in inds}:
            for y in {b[i] for i in inds}:
                row = sum(w[i] for i in inds if a[i] == x)
                col = sum(w[i] for i in inds if b[i] == y)
                cell = sum(w[i] for i in inds if a[i] == x and b[i] == y)
                if total*cell != row*col:
                    return False
    return True


def descends(t, b):
    # Direct strong-lumpability check, independent of projection multiplication.
    for i in range(len(b)):
        for j in range(i):
            if b[i] == b[j]:
                for label in set(b):
                    if sum(t[i][k] for k in range(len(b)) if b[k] == label) != sum(
                       t[j][k] for k in range(len(b)) if b[k] == label):
                        return False
    return True


def mass_constant_on_join(a,b,w):
    c = join(a,b)
    ma = [sum(w[j] for j in range(len(w)) if a[j] == a[i]) for i in range(len(w))]
    return all(c[i] != c[j] or ma[i] == ma[j] for i in range(len(w)) for j in range(i))


def det(a):
    if not a:
        return 1
    return sum((-1)**j*a[0][j]*det([r[:j]+r[j+1:] for r in a[1:]]) for j in range(len(a)))


def determinantal_divisors(a):
    n = len(a)
    return [reduce(gcd,(abs(det([[a[i][j] for j in cols] for i in rows]))
                       for rows in combinations(range(n),r)
                       for cols in combinations(range(n),r)),0)
            for r in range(1,n+1)]


def rank_mod(a,p):
    a = [[int(x)%p for x in r] for r in a]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank,len(a)) if a[r][col]),None)
        if pivot is None:
            continue
        a[rank],a[pivot] = a[pivot],a[rank]
        inv = pow(a[rank][col],-1,p)
        a[rank] = [x*inv%p for x in a[rank]]
        for r in range(len(a)):
            if r != rank:
                c = a[r][col]
                a[r] = [(x-c*y)%p for x,y in zip(a[r],a[rank])]
        rank += 1
    return rank


def normalized_push(f,w):
    return [[F(w[i],sum(w[j] for j in range(len(w)) if f[j] == y)) if f[i] == y else 0
             for i in range(len(w))] for y in range(max(f)+1)]


def pull(f):
    return [[int(y == f[i]) for y in range(max(f)+1)] for i in range(len(f))]


def main():
    count = 0
    commuting = 0
    for n in range(1,6):
        ps = list(partitions(n))
        weights = list(product((1,2),repeat=n)) if n <= 4 else [(1,)*n,tuple(range(1,n+1)),(1,2,1,2,1)]
        for w in weights:
            ops = {a:operator(a,w) for a in ps}
            raw = {a:operator(a,w,False) for a in ps}
            for a in ps:
                for b in ps:
                    pa,pb = ops[a],ops[b]
                    ab,ba = mm(pa,pb),mm(pb,pa)
                    cond = factors(a,b,w)
                    require((ab == ba) == cond, ('commutation',a,b,w))
                    require(descends(pa,b) == cond, ('normalized descent',a,b,w))
                    require(descends(raw[a],b) == (cond and mass_constant_on_join(a,b,w)),
                            ('unnormalized descent',a,b,w))
                    if cond:
                        require(ab == operator(join(a,b),w), ('join product',a,b,w))
                    count += 1
                    commuting += cond
    arithmetic = 0
    for k in range(1,5):
        for m in product(range(1,5),repeat=k):
            M,g = sum(m),reduce(gcd,m)
            q = [list(m) for _ in m]
            require(mm(q,q) == scale(M,q), ('quasi-idempotence',m))
            require(determinantal_divisors(q) == [g]+[0]*(k-1), ('Q Smith',m))
            lap = minus(scale(M,ident(k)),q)
            expected = [g*M**(r-1) for r in range(1,k)]+[0]
            require(determinantal_divisors(lap) == expected, ('L Smith',m))
            denom = reduce(lcm,(F(x,M).denominator for x in m),1)
            require(denom == M//g, ('denominator',m))
            for p in (2,3,5,7):
                r = rank_mod(q,p)
                require(r == int(any(x%p for x in m)), ('mod rank',m,p))
                if M%p == 0:
                    require(all(x%p == 0 for row in mm(q,q) for x in row),('square zero',m,p))
                    require(k-2*r >= 0, ('homology dimension',m,p))
            arithmetic += 1
    # Minimal noncommuting example and its rational principal-angle witness.
    a,b,w = (0,0,1),(0,1,1),(1,1,1)
    p,q = operator(a,w),operator(b,w)
    j = operator((0,0,0),w)
    u = [[1],[1],[-2]]
    require(mm(mm(mm(p,q),p),u) == scale(F(1,4),u), 'lambda=1/4')
    require(mm(j,u) == [[0],[0],[0]], 'join kills witness')
    s = mm(mm(p,q),p)
    require(scale(F(1,3),minus(scale(4,s),p)) == j, 'signed rational recovery')
    require(scale(F(1,3),minus(scale(4,mm(s,s)),s)) == j, 'whole-space spectral polynomial')
    for length in range(1,13):
        for first in (0,1):
            t = ident(3)
            for i in range(length):
                t = mm((p,q)[(i+first)%2],t)
            require(t != j, ('finite-word hostile',length,first))
    # Modular positivity hostile: all masses are units, but PRP=join and PR!=RP.
    def mod(a,p):
        return [[(x.numerator*pow(x.denominator,-1,p))%p if isinstance(x,F) else x%p for x in r] for r in a]
    ma,mb,mw = (0,0,1,0),(0,1,0,2),(1,1,1,2)
    mp,mr = mod(operator(ma,mw),3),mod(operator(mb,mw),3)
    mj = mod(operator((0,0,0,0),mw),3)
    require(mod(mm(mm(mp,mr),mp),3) == mj, 'mod3 three-letter exact completion')
    require(mod(mm(mp,mr),3) != mod(mm(mr,mp),3), 'mod3 noncommutation')
    # Literal divisor maps; full product law versus a deleted top corner.
    for prime,other in ((2,3),(2,7),(3,7),(5,7)):
        da,db,dw = (0,1,0,1),(0,0,1,1),(1,prime,other,prime*other)
        require(factors(da,db,dw), ('full divisor product',prime,other))
        da,db,dw = (0,1,0),(0,0,1),(1,prime,other)
        require(not factors(da,db,dw), ('proper divisor hostile',prime,other))
        dp,dr = operator(da,dw),operator(db,dw)
        du = [[1],[-F(1+other,prime)],[1]]
        lam = F(prime*other,(1+prime)*(1+other))
        require(mm(mm(mm(dp,dr),dp),du) == scale(lam,du), ('divisor lambda',prime,other))
        require(mm(operator((0,0,0),dw),du) == [[0],[0],[0]], ('divisor mean',prime,other))
    # Same row/column margins but opposite commutation verdicts: joint sidecar matters.
    # Use connected positive tables with row/column margins (6,6).
    aa,bb = (0,0,1,1),(0,1,0,1)
    require(factors(aa,bb,(3,3,3,3)), 'same-margin independent 2x2')
    require(not factors(aa,bb,(2,4,4,2)), 'same-margin correlated 2x2')
    # Tower law with actual pushed source masses, not uniform intermediate fibers.
    f,g = (0,0,1,1,2),(0,0,1)
    w = (1,2,3,4,5)
    wf = tuple(sum(w[i] for i in range(len(w)) if f[i] == y) for y in range(3))
    h = tuple(g[y] for y in f)
    require(mm(normalized_push(g,wf),normalized_push(f,w)) == normalized_push(h,w), 'weighted tower')
    pf = mm(pull(f),normalized_push(f,w))
    ph = mm(pull(h),normalized_push(h,w))
    require(mm(pf,ph) == ph == mm(ph,pf), 'nested expectations')
    # Cartesian base-change: duplicate selected target fibers.
    hbase = (2,0,2,1)
    pairs = [(yp,x) for yp,y in enumerate(hbase) for x,z in enumerate(f) if z == y]
    fprime = tuple(yp for yp,x in pairs)
    hprime = tuple(x for yp,x in pairs)
    wp = tuple(w[x] for yp,x in pairs)
    hpull = [[int(x == old) for old in range(len(f))] for x in hprime]
    hbpull = [[int(y == old) for old in range(3)] for y in hbase]
    require(mm(normalized_push(fprime,wp),hpull) == mm(hbpull,normalized_push(f,w)), 'base change')
    # gcd-cancellation lift distinguishes P mod 3 from Q mod 3.
    require(rank_mod([[3]],3) == 0 and F(3,3) == 1, 'lift-data hostile')
    # Variable-fiber-size raw transfer does not descend, normalized transfer does.
    va,vb,vw = (0,0,1),(0,0,0),(1,1,1)
    require(descends(operator(va,vw),vb) and not descends(operator(va,vw,False),vb), 'mass sidecar hostile')
    # All prime pairs through 100: direct prime supports versus proved criterion.
    def is_prime(n):
        return n >= 2 and all(n%d for d in range(2,isqrt(n)+1))
    def prime_support(n):
        result = set()
        d = 2
        while d*d <= n:
            if n%d == 0:
                result.add(d)
                while n%d == 0:
                    n //= d
            d += 1
        if n > 1:
            result.add(n)
        return result
    primes = [p for p in range(2,101) if is_prime(p)]
    no_new_prime = []
    prime_pairs = 0
    for prime,other in combinations(primes,2):
        total = 1+prime+other
        remainder = total
        while remainder%prime == 0:
            remainder //= prime
        predicted = (prime,other) == (2,3) or remainder == 1
        actual = prime_support(total) <= prime_support((1+prime)*(1+other))
        require(predicted == actual, ('prime denominator boundary',prime,other))
        if actual:
            no_new_prime.append([prime,other,total])
        prime_pairs += 1
    # Inherited excess-two family: independent divisor-pair sum of n.
    radical_bridge = []
    for b in range(2,11):
        prime = 2**(b+1)-3
        if not is_prime(prime):
            continue
        n = 2**b*prime
        sigma = sum(d+(n//d if d*d != n else 0) for d in range(1,isqrt(n)+1) if n%d == 0)
        require(sigma-2*n == 2, ('excess-two member',b,n))
        require(3+prime == 2**(b+1), ('radical image sum',b,n))
        require(prime_support(3+prime) <= prime_support(3*(1+prime)), ('radical denominator',b,n))
        radical_bridge.append([b,n,2*prime,3+prime])
    # All ordered unitary splits in the declared sample, full coordinate spaces.
    def divisors(n):
        return sorted({d for a in range(1,isqrt(n)+1) if n%a == 0 for d in (a,n//a)})
    def rank_q(matrix):
        rows = [[F(x) for x in r] for r in matrix]
        rank = 0
        for col in range(len(rows[0])):
            pivot = next((r for r in range(rank,len(rows)) if rows[r][col]),None)
            if pivot is None:
                continue
            rows[rank],rows[pivot] = rows[pivot],rows[rank]
            divisor = rows[rank][col]
            rows[rank] = [x/divisor for x in rows[rank]]
            for r in range(rank+1,len(rows)):
                factor = rows[r][col]
                rows[r] = [x-factor*y for x,y in zip(rows[r],rows[rank])]
            rank += 1
        return rank
    split_controls = []
    for n in (6,10,12,14,18,20,28,30,36,60,210,496):
        dn = divisors(n)
        sigma_n = sum(dn)
        proper_n = sigma_n-n
        splits = []
        for a in dn[1:-1]:
            b = n//a
            if gcd(a,b) != 1:
                continue
            da,db = divisors(a),divisors(b)
            pairs = [(u,v) for u in da for v in db if (u,v) != (a,b)]
            weights = tuple(u*v for u,v in pairs)
            require(sorted(weights) == dn[:-1], ('unitary multiplication',n,a,b))
            flags = tuple(1 if u == a else 2 if v == b else 0 for u,v in pairs)
            sx,sy = sum(da)-a,sum(db)-b
            fw = (sx*sy,a*sy,sx*b)
            require(tuple(sum(w for w,label in zip(weights,flags) if label == j) for j in range(3)) == fw,
                    ('flag masses',n,a,b))
            p3,r3 = operator((0,1,0),fw),operator((0,0,1),fw)
            e3 = operator((0,0,0),fw)
            s3 = mm(mm(p3,r3),p3)
            pzero,qzero = minus(p3,e3),minus(r3,e3)
            require(mm(pzero,pzero) == pzero and mm(qzero,qzero) == qzero, ('mean-zero idempotents',n,a,b))
            require(mm(mm(pzero,qzero),pzero) == scale(F(n,sigma_n),pzero), ('first finite TL relation',n,a,b))
            require(mm(mm(qzero,pzero),qzero) == scale(F(n,sigma_n),qzero), ('second finite TL relation',n,a,b))
            require(scale(F(1,proper_n),minus(scale(sigma_n,s3),scale(n,p3))) == e3,
                    ('three-flag projector',n,a,b))
            pa,pb = operator(tuple(u for u,v in pairs),weights),operator(tuple(v for u,v in pairs),weights)
            h = [[int(label == j) for j in range(3)] for label in flags]
            require(mm(pa,h) == mm(h,p3) and mm(pb,h) == mm(h,r3), ('flag intertwinings',n,a,b))
            ss = mm(mm(pa,pb),pa)
            ss2 = mm(ss,ss)
            ee = operator((0,)*len(weights),weights)
            require(scale(F(1,proper_n),minus(scale(sigma_n,ss2),scale(n,ss))) == ee,
                    ('full-coordinate projector',n,a,b))
            require(rank_q(ss) == 2, ('full rank two',n,a,b))
            require(sum(ss[i][i] for i in range(len(ss))) == 1+F(n,sigma_n), ('full trace',n,a,b))
            require(minus(scale(sigma_n,mm(ss2,ss)),minus(scale(sigma_n+n,ss2),scale(n,ss))) == scale(0,ss),
                    ('full spectral polynomial',n,a,b))
            lam = F(fw[1]*fw[2],(fw[0]+fw[1])*(fw[0]+fw[2]))
            require(lam == F(n,sigma_n), ('split-independent lambda',n,a,b))
            require((1-lam)/lam == F(proper_n,n), ('abundancy angle',n,a,b))
            require((lam == F(1,2)) == (proper_n == n), ('perfect angle',n,a,b))
            require((lam < F(1,2)) == (proper_n > n), ('abundant angle',n,a,b))
            if (n,a,b) == (12,4,3):
                hostile = [[{1:2,2:-1,4:0}[u]] for u,v in pairs]
                require(mm(pa,hostile) == hostile and mm(ss,hostile) == [[0] for _ in pairs], 'extra A zero mode')
                wrong = scale(F(1,proper_n),minus(scale(sigma_n,ss),scale(n,pa)))
                require(mm(wrong,hostile) == scale(F(-3,4),hostile) != mm(ee,hostile), 'wrong full formula hostile')
            splits.append([a,b])
        split_controls.append({'N':n,'sigma':sigma_n,'s':proper_n,'splits':splits,'lambda':str(F(n,sigma_n))})
    require(sum(divisors(6))**2 == 144 and sum(divisors(36)) == 91, 'noncoprime full mass')
    require(144-36 == 108 and 91-36 == 55, 'noncoprime proper mass')
    require(F(36*36,72*72) == F(1,4) != F(36,91), 'noncoprime angle hostile')
    result = {
        'status':'PASS', 'partition_weight_pairs':count, 'commuting_cases':commuting,
        'prime_pairs_through_100':prime_pairs, 'no_new_prime_cases':no_new_prime,
        'radical_bridge_b_n_radical_target':radical_bridge,
        'unitary_split_controls':split_controls,
        'arithmetic_blocks':arithmetic, 'primes':[2,3,5,7],
        'partition_universe':'all ordered pairs of set partitions for n=1..4, all weights in {1,2}^n; n=5 all partition pairs for three stated weight vectors',
        'arithmetic_universe':'all positive integer vectors of lengths 1..4 with entries 1..4; all square minors',
        'hostiles':['3-state noncommutation and eigenvalue 1/4','24 alternating words through length 12','same-margins different joint tables','gcd-cancellation lift','variable-mass raw descent','F3 isotropic exact completion','proper divisor missing corner'],
        'composition_controls':['weighted tower','nested expectations','Cartesian weighted base-change','signed rational recovery','full divisor product weights'],
        'dependencies':'Python standard library only; checks remain active under -O'}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
