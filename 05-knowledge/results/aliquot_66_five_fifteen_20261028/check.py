"""Exact checks for the typed bridges; no removable assertions or dependencies."""
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
from fractions import Fraction
from pathlib import Path
import json
import hashlib

HERE = Path(__file__).resolve().parent
CHECKS = 0

def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def sigma(n):
    return prod((p ** (a + 1) - 1) // (p - 1) for p, a in factor(n).items())

def sf(n):
    return prod(p for p, a in factor(n).items() if a % 2)

def square(n):
    return isqrt(n) ** 2 == n

def prime(n):
    return n >= 2 and all(n % p for p in range(2, isqrt(n) + 1))

def balance(n):
    f = factor(n)
    r = len(f)
    F = prod(a + 1 for a in f.values()) - 2
    S = 2 ** r - 1 - int(all(a == 1 for a in f.values()))
    U = r - int(len(f) == 1 and next(iter(f.values())) == 1)
    return F - S - U

def rank(signatures):
    primes = sorted({p for x in signatures for p in factor(x)})
    pivots = {}
    for x in signatures:
        v = sum((1 << i) for i, p in enumerate(primes) if x % p == 0)
        while v:
            b = v.bit_length()
            if b in pivots:
                v ^= pivots[b]
            else:
                pivots[b] = v
                break
    return len(pivots)

def psi(n):
    return prod(p ** (a-1) * (p+1) for p,a in factor(n).items())

def hnf(n):
    return [(a,b,n//a) for a in range(1,n+1) if n%a == 0 for b in range(a)]

def reduced_forms(D):
    out = []
    for a in range(1,isqrt((-D)//3)+1):
        for b in range(-a,a+1):
            if (b*b-D) % (4*a):
                continue
            c = (b*b-D)//(4*a)
            if a > c or gcd(gcd(a,abs(b)),c) != 1:
                continue
            if (abs(b)==a or a==c) and b < 0:
                continue
            out.append((a,b,c))
    return out

def genus(n):
    # Elliptic points counted as roots, independently of Legendre-product formulas.
    e2 = sum((x*x+1)%n == 0 for x in range(n))
    e3 = sum((x*x+x+1)%n == 0 for x in range(n))
    def phi(d):
        return sum(gcd(i,d)==1 for i in range(1,d+1))
    cusps = sum(phi(gcd(d,n//d)) for d in range(1,n+1) if n%d==0)
    g = 1 + Fraction(psi(n),12) - Fraction(e2,4) - Fraction(e3,3) - Fraction(cusps,2)
    check(g.denominator == 1 and g>=0, f"genus integer {n}")
    return {"index":psi(n),"cusps":cusps,"e2":e2,"e3":e3,"genus":int(g)}

def main():
    limit = 200000
    sums = [0] * (limit + 1)
    for d in range(1, limit + 1):
        for n in range(d, limit + 1, d):
            sums[n] += d
    for n in range(1, 5001):
        check(sums[n] == sigma(n), f"independent sigma {n}")
    squares = [n for n in range(1, limit + 1) if square(sums[n])]
    primes = [p for p in range(2, 98) if prime(p)]
    edges = []
    for n in squares:
        for p in primes:
            actual = square(sigma(n * p))
            check(actual == (p == 3 and n % 3 != 0), f"prime dilation {n},{p}")
            if actual:
                edges.append([n, n * p])
    cubes = {}
    for support in ((2, 3, 11), (2, 5, 7)):
        rows = []
        for mask in product((0, 1), repeat=3):
            n = prod(p ** (1 + bit) for p, bit in zip(support, mask))
            rows.append({"mask": list(mask), "n": n, "sigma": sigma(n),
                         "signature": sf(sigma(n)), "balance_defect": balance(n)})
            check((balance(n) == 0) == (sum(mask) == 1), f"FUS cube {n}")
            check(square(sigma(n)) == (sum(mask) == 0), f"square cube {n}")
        jumps = [sf((p + 1) * (p * p + p + 1)) for p in support]
        check(rank(jumps) == 3, f"independent jumps {support}")
        check(len({r['signature'] for r in rows}) == 8, f"distinct signatures {support}")
        cubes[str(support)] = {"jumps": jumps, "rank": rank(jumps), "rows": rows}
    # Literal divisors, independent of the exponent-box formulas.
    for n in range(2, 1501):
        ds = [d for d in range(2, n) if n % d == 0]
        S = sum(all(a == 1 for a in factor(d).values()) for d in ds)
        U = sum(prime(d) for d in ds)
        check(len(ds) - S - U == balance(n), f"literal balance {n}")
    mixed = []
    for n in range(2, limit + 1):
        if not square(sums[n]):
            continue
        f = factor(n)
        if sorted(f.values()) == [1, 1, 2]:
            mixed.append(n)
            p = next(p for p, a in f.items() if a == 2)
            q, r = sorted(p for p, a in f.items() if a == 1)
            d = sf((p*p + p + 1) * (q + 1))
            check((r + 1) % d == 0 and square((r + 1) // d), f"parametrization {n}")
    check(mixed[0] == 364, "first mixed balanced square sigma")
    check(sigma(364) == 784 and sigma(182) == 336, "repair hostile")
    for n in squares:
        fs=factor(n)
        if len(fs)!=3 or any(a!=1 for a in fs.values()):
            continue
        ps=tuple(fs)
        zero_subsets=[]
        for size in (1,2):
            zero_subsets += [subset for subset in combinations(ps,size) if square(sigma(prod(subset)))]
        check(bool(zero_subsets)==(3 in fs),f"three-column dependence type {n}")
    check(7*19*89==11837 and sigma(11837)==14400,
          "incoming Frobenius preserves square total")
    check(2*2-4*3*3946==-4*11837 and gcd(gcd(3,2),3946)==1 and 0<2<3<3946,
          "incoming Frobenius does not preserve idoneality")
    # Root-leg map: all primitive Euclid pairs in this finite parameter box.
    root_rows = []
    hostile = None
    for m in range(2, 61):
        for k in range(1, m):
            if gcd(m, k) != 1 or (m - k) % 2 != 1:
                continue
            a, b, c = m*m-k*k, 2*m*k, m*m+k*k
            check(b % 4 == 0 and a*a + b*b == c*c, "primitive marked leg")
            q = 3*(b//2)**2 - 1
            n = 2*q
            if prime(q):
                check(sigma(n) == (3*b//2)**2, "leg square-sigma map")
                check(sigma(3*n) == (3*b)**2, "leg adjoin-three lift")
                root_rows.append([a,b,c,n,sigma(n)-n,3*n,sigma(3*n)-3*n])
            elif hostile is None:
                hostile = {"triple":[a,b,c],"q":q,"q_factors":factor(q),
                           "actual_sigma":sigma(n),"unguarded_total":(3*b//2)**2}
    check(root_rows[0] == [3,4,5,22,14,66,78], "root map")
    check(hostile['actual_sigma'] != hostile['unguarded_total'], "prime guard hostile")
    # A5 matching bridge, retaining the marked sixth label.
    def matchings(xs):
        if not xs:
            yield ()
        else:
            a = xs[0]
            for b in xs[1:]:
                rest = tuple(x for x in xs if x not in (a,b))
                for tail in matchings(rest):
                    yield tuple(sorted(((a,b),) + tail))
    ms = list(matchings(tuple(range(6))))
    ds = list(combinations(range(5),2))
    images = set()
    for matching in ms:
        perm = list(range(5))
        for a,b in matching:
            if b != 5:
                perm[a],perm[b] = b,a
        images.add(tuple(perm))
    check(len(ms) == len(images) == 15, "marked A5 involutions")
    # K5 two-subsets with disjointness is Petersen; its odd cycle forbids a sign gauge.
    pedges = [(i,j) for i in range(10) for j in range(i+1,10) if not set(ds[i]) & set(ds[j])]
    check(len(pedges) == 15, "Petersen edges")
    bip_solutions = sum(all(s[i] != s[j] for i,j in pedges) for s in product((0,1),repeat=10))
    check(bip_solutions == 0, "Petersen cannot carry sum-to-difference vertex sign gauge")
    cycle = [(0,1),(2,3),(0,4),(1,2),(3,4)]
    check(all(not set(cycle[i]) & set(cycle[(i+1)%5]) for i in range(5)), "odd cycle witness")
    # Base-height normalization and prime partition.
    exps = {2:5,3:7,5:3,7:4,11:1,13:2}
    def height(b, es):
        return min(es.get(p,0)//a for p,a in factor(b).items())
    check([height(b,exps) for b in (4,2,5,3,7,13)] == [2,5,3,7,4,2], "height cycle")
    corrected = dict(exps, **{})
    corrected[2] -= 1
    check(height(2,corrected) == 4 and height(4,corrected) == 2, "centralizer normalization")
    # HNF enumeration is independent of product sigma and psi formulas.
    for n in range(1,501):
        mats = hnf(n)
        check(len(mats) == sigma(n), f"lattice mass {n}")
        primitive = sum(gcd(gcd(a,b),d)==1 for a,b,d in mats)
        check(primitive == psi(n), f"primitive lattice mass {n}")
        contents = {}
        for a,b,d in mats:
            t=gcd(gcd(a,b),d)
            contents[t]=contents.get(t,0)+1
        expected={d:psi(n//(d*d)) for d in range(1,isqrt(n)+1) if n%(d*d)==0}
        check(contents==expected, f"content strata {n}")
    check((8,0,1) in hnf(8) and (8,2,1) in hnf(8), "mod-two invisible shear")
    chain_fibres={}
    for p in (2,3,5,7):
        rows=[]
        for a,b,c in hnf(p*p):
            count=sum(a%A==0 and c%C==0 and (b-B*(c//C))%A==0 for A,B,C in hnf(p))
            expected=p+1 if (a,b,c)==(p,0,p) else 1
            check(count==expected,f"two-step fibre {p},{a},{b},{c}")
            rows.append([a,b,c,count])
        check(sum(r[3] for r in rows)==(p+1)**2,f"unmerged path mass {p}")
        chain_fibres[str(p)]=rows
    # Direct isotropic-span enumeration, independent of the double-count formula.
    lag_counts={}
    for p in (2,3):
        vs=list(product(range(p),repeat=4))[1:]
        spaces=set()
        for u in vs:
            for v in vs:
                pairing=(u[0]*v[2]+u[1]*v[3]-u[2]*v[0]-u[3]*v[1])%p
                if pairing:
                    continue
                span=tuple(sorted({tuple((a*u[i]+b*v[i])%p for i in range(4))
                                   for a in range(p) for b in range(p)}))
                if len(span)==p*p:
                    spaces.add(span)
        check(len(spaces)==sigma(p**3),f"Lagrangian and lattice count {p}")
        lag_counts[str(p)]=len(spaces)
    forms = {}
    for n, expected in ((5,(2,2)),(15,(2,2)),(21,(4,4)),(63,(4,2)),(66,(8,4)),(70,(4,4)),(1848,(16,16))):
        fs = reduced_forms(-4*n)
        amb = [f for f in fs if f[1]==0 or abs(f[1])==f[0] or f[0]==f[2]]
        check((len(fs),len(amb)) == expected, f"reduced form carrier {n}")
        forms[str(n)]={"discriminant":-4*n,"forms":fs,"self_inverse":len(amb),
                       "genus_kernel_size":len(fs)//len(amb)}
    totals = [n for n in range(1,144) if sigma(n)==144]
    check(totals==[66,70,94,115,119],"complete sigma144 fibre")
    genera={str(n):genus(n) for n in totals+[132,198,726,140,350,490,364]}
    check([genera[str(n)]['genus'] for n in totals]==[9,9,11,11,11],"two genera of same total")
    for n in (66,70):
        for p in factor(n):
            check(genera[str(n*p)]['genus']==10*p-1,f"promoted genus {n},{p}")
    ray_classes=[]
    check(reduced_forms(-8)==[(1,0,2)],"parent family small idoneal boundary")
    check((3,2,9) in reduced_forms(-104) and (3,-2,9) in reduced_forms(-104),
          "parent family t2")
    for t in range(3,1001):
        a,b,c=2*t-1,4,4*t+2
        check(b*b-4*a*c == -4*(8*t*t-6) and 0<b<a<c and gcd(gcd(a,b),c)==1,
              f"parent family nonambiguous form {t}")
    for k in range(1,1001):
        q=12*k*k-1
        a,b,c=6*k-1,4,12*k+2
        check(b*b-4*a*c == -24*q and 0<b<a<c and gcd(gcd(a,b),c)==1,
              f"uniform nonambiguous form {k}")
        if k<=30 and prime(q):
            fs=reduced_forms(-24*q)
            check((a,b,c) in fs and (a,-b,c) in fs, f"inverse pair on ray {k}")
            amb=sum(y==0 or abs(y)==x or x==z for x,y,z in fs)
            check(amb==4 and len(fs)>=8,f"nonzero genus kernel {k}")
            ray_classes.append({"k":k,"q":q,"class_number":len(fs),"genus_kernel":len(fs)//4})
    result = {"status":"FINITE-EXACT controls, not an independent proof audit",
              "checks":CHECKS,"sigma_source_limit":limit,"square_sigma_sources":len(squares),
              "tested_multiplier_primes":primes,"prime_edges":len(edges),
              "first_edges":edges[:12],"promotion_cubes":cubes,
              "mixed_balanced_square_sources":mixed,
              "marked_leg_parameter_limit":60,"marked_leg_prime_rows":len(root_rows),
              "root_rows_first":root_rows[:8],"marked_leg_hostile":hostile,
              "matching_count":len(ms),"petersen_edges":len(pedges),
              "class_forms":forms,"modular_curves":genera,
              "berggren_ray_class_groups":ray_classes,
              "two_step_lattice_fibres":chain_fibres,"lagrangian_counts":lag_counts,
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__ == "__main__":
    main()
