"""Exact cyclic-repair audit. All gates remain active under python -O.

Universe: integer mean-zero vectors with first n-1 coordinates in {-1,0,1},
n=3..7; all reflection-symmetric such vectors, including even-cycle controls.
Capacity tests: every nonnegative {0,1,2} edge-capacity vector at n=3,4 and
all source vectors from that universe. Analytic profile KKT certificates:
every odd p=5..31 and the explicitly listed rational capacity samples.
No canonical owner/word realization is assumed.
"""
from fractions import Fraction as F
from itertools import product
import hashlib
import json

GATES = 0


def check(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def lap(x):
    return [2*x[s]-x[s-1]-x[(s+1)%len(x)] for s in range(len(x))]


def div(f):
    return [f[s]-f[s-1] for s in range(len(f))]


def prefixes(a):
    total = 0
    out = []
    for t in a:
        total += t
        out.append(total)
    check(total == 0, "prefix domain")
    return out


def potential(a):
    # Solves Lap(x)=-a over Q with x[0]=0, using zero-circulation flux.
    p = len(a)
    A = prefixes(a)
    h = F(sum(A), p)
    f = [h-t for t in A]
    x = [F(0)]
    for t in f[:-1]:
        x.append(x[-1]-t)
    check(sum(f) == 0, "zero circulation")
    check(lap(x) == [-t for t in a], "potential reconstruction")
    return x, f


def gauss_potential(a):
    # Independent reduced-Laplacian solve, no cumulative-coordinate formula.
    n = len(a)
    basis = []
    for j in range(1, n):
        e = [F(0)] * n
        e[j] = F(1)
        basis.append(lap(e))
    M = [[basis[j-1][i] for j in range(1,n)] + [F(-a[i])]
         for i in range(1,n)]
    for j in range(n-1):
        pivot = next(i for i in range(j,n-1) if M[i][j])
        M[j], M[pivot] = M[pivot], M[j]
        q = M[j][j]
        M[j] = [v/q for v in M[j]]
        for i in range(n-1):
            if i != j:
                q = M[i][j]
                M[i] = [u-q*v for u,v in zip(M[i], M[j])]
    return [F(0)] + [M[i][-1] for i in range(n-1)]


def capacity_feasible(a, k):
    A = prefixes(a)
    return max(t-c for t,c in zip(A,k)) <= min(t+c for t,c in zip(A,k))


def arc_capacity_feasible(a, k):
    n = len(a)
    for start in range(n):
        total = 0
        for length in range(1,n):
            end = (start+length-1)%n
            total += a[end]
            if abs(total) > k[(start-1)%n] + k[end]:
                return False
    return True


def role(p):
    check(p >= 5 and p%2 == 1, "role domain")
    a0 = (p-1)**2
    a1 = (p*p-5*p+2)//2
    b = 2*p-1
    return [a0,a1] + [-b]*(p-3) + [a1]


def optimal_profile(p,k):
    a = role(p)
    A = (p-3)*(2*p-1)
    k0 = F(p*(p+1),4)
    ks = F(A,2)
    if k <= k0:
        c = [a[0]-2*k,a[1]] + [F(-(2*p-1))+2*k/F(p-3)]*(p-3) + [a[1]]
    elif k <= ks:
        c = [F(A-2*k,3)]*2 + [F(-(A-2*k),p-3)]*(p-3) + [F(A-2*k,3)]
    else:
        c = [F(0)]*p
    x,f = potential([u-v for u,v in zip(a,c)])
    check(all(abs(v) <= k for v in f), "optimizer capacity")
    check([u+v for u,v in zip(a,div(f))] == c, "optimizer profile")
    for s in range(p):
        delta = c[s]-c[(s+1)%p]
        check((delta == 0) or (delta > 0 and f[s] == -k)
              or (delta < 0 and f[s] == k), "box first-variation certificate")
    return c,f


def main():
    audit = {}
    torsion_cases = 0
    symmetric_counts = {}
    for n in range(3,8):
        syms = 0
        for first in product((-1,0,1),repeat=n-1):
            a = list(first)+[-sum(first)]
            x,f = potential(a)
            q = sum(s*a[s] for s in range(n))%n
            check(all(v.denominator == 1 for v in x) == (q == 0), "torsion iff")
            if n <= 5:
                check(x == gauss_potential(a), "independent reduced solve")
            if all(a[s] == a[-s%n] for s in range(n)):
                syms += 1
                check((2*q)%n == 0, "reflection reverses critical class")
                if n%2:
                    check(q == 0, "odd reflection kills torsion")
            torsion_cases += 1
        symmetric_counts[n] = syms
    audit["torsion_vectors"] = torsion_cases
    audit["symmetric_vectors_by_order"] = symmetric_counts
    check(sum(s*v for s,v in enumerate((1,0,-1,0)))%4 == 2,
          "even symmetric torsion hostile")
    check(sum(s*v for s,v in enumerate((1,-1,0)))%3 != 0,
          "nonsymmetric odd torsion hostile")
    capacity_cases = 0
    for n in (3,4):
        for first in product((-1,0,1),repeat=n-1):
            a = list(first)+[-sum(first)]
            for k in product((0,1,2),repeat=n):
                check(capacity_feasible(a,k) == arc_capacity_feasible(a,k),
                      "all arcs versus interval intersection")
                capacity_cases += 1
    audit["capacity_instances"] = capacity_cases
    analytic_cases = 0
    for p in range(5,32,2):
        k0 = F(p*(p+1),4)
        ks = F((p-3)*(2*p-1),2)
        samples = sorted({F(0),F(1,3),k0/2,k0,k0+F(1,4),(k0+ks)/2,
                          ks-F(1,3),ks,ks+1,2*ks})
        for k in samples:
            c,f = optimal_profile(p,k)
            if k >= k0 and k <= ks:
                check(sum(v*v for v in c)/p == (2*ks-2*k)**2/F(3*(p-3)),
                      "sharp arc-energy formula")
            analytic_cases += 1
        a = role(p)
        x,f = potential(a)
        check(all(v.denominator == 1 for v in x), "role integral repair")
        check(max(abs(v) for v in f) == ks, "canonical minimum flux")
        check(max(x)-min(x) == F((p-1)*(2*p*p-3*p-1),8),
              "canonical potential range")
    audit["analytic_capacity_certificates"] = analytic_cases
    p=13
    a=role(p)
    x,f=potential(a)
    check(a == [144,53]+[-25]*10+[53], "p13 role")
    check(x == list(map(F,[0,72,197,297,372,422,447,447,422,372,297,197,72])),
          "p13 exact potential")
    c,_=optimal_profile(p,F(100))
    D=1183
    audit["p13"]={"a":a,"D":D,"torsion_class":sum(s*a[s] for s in range(p))%p,
                   "potential":list(map(str,x)),"cancel_flux":list(map(str,f)),
                   "k0":"91/2","k_cancel":"125","potential_range":"447",
                   "variance_at_k100":str(sum(t*t for t in c)/F(p*D*D)),
                   "variance_at_zero_capacity":str(sum(F(t*t) for t in a)/F(p*D*D))}
    # Independent integer-flux exhaustive check for a small member, where
    # the analytic optimum happens to be integral at the chosen capacities.
    brute = 0
    a5=role(5)
    for k in (0,1,2):
        best=None
        for f0 in product(range(-k,k+1),repeat=5):
            c0=[u+v for u,v in zip(a5,div(f0))]
            value=sum(v*v for v in c0)
            best=value if best is None else min(best,value)
            brute += 1
        c5,_=optimal_profile(5,F(k))
        check(best == sum(v*v for v in c5), "independent integer flux optimum")
    audit["brute_flux_vectors"] = brute
    audit["gates"] = GATES
    body=json.dumps(audit,sort_keys=True,indent=2)
    print("CYCLIC REPAIR: PROVED SYMBOLIC + FINITE-EXACT AUDIT; LRC realization OPEN")
    print(body)
    print("semantic_sha256="+hashlib.sha256(body.encode()).hexdigest())
    print("PASS: explicit gates remain active with -O")


if __name__ == "__main__":
    main()
