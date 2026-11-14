"""Small strong-core census and literal graph replays for Session 330."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import canon, ends, extend, glue as fermat_glue, rooted, strong


def slots(A, top):
    h = [sum(row) for row in ends(A)]
    h[0] = 1
    prev = h
    out = [0, h[-1]]
    S = len(h)
    for z in range(2, top + 1):
        f = [0] * S
        for s in range(S):
            t = s
            while True:
                f[s] += h[t] * prev[s ^ t]
                if not t:
                    break
                t = (t - 1) & s
        prev = f
        out.append(f[-1])
    return out


def ordinal(A, B):
    n = len(A)
    return (tuple(x | (((1 << len(B)) - 1) << n) for x in A)
            + tuple(x << n for x in B))


def triangle(A, B):
    n, m = len(A), len(B)
    return (tuple(x | (((1 << m) - 1) << n) for x in A)
            + tuple((x << n) | (1 << (n + m)) for x in B)
            + ((1 << n) - 1,))


def scalar(a, b):
    h, q = a
    k, r = b
    return (q - h) * (r - k) + 2 * h * k


def verify(flag, why):
    if not flag:
        raise RuntimeError(why)


def main():
    ps, F = {}, {}
    T = {(0,)}
    bounds = (1, 1, 2, 4, 12, 56, 456, 6880)
    for n, count in enumerate(bounds, 1):
        if n > 1:
            T = {canon(extend(A, c)) for A in T for c in range(1 << (n-1))}
        verify(len(T) == count, ("class count", n))
        ps[n], F[n] = set(), set()
        for A in T:
            if not strong(A):
                continue
            if n <= 7:
                ps[n].add(tuple(slots(A, 2)[1:]))
            if n <= 7:
                F[n].update(h + (l << j) for h,l,q,_,_,root in rooted(A) for j in (1,2,4))
        if n >= 7:
            print("isomorphism classes", n, len(T), flush=True)

    A4, A7 = (2,4,9,3), (6,36,104,19,7,89,27)
    B5, B6 = (6,12,24,17,3), (6,12,24,49,35,7)
    pairs = [[(a,b) for a in ps[u] for b in ps[v]] for u,v in [(4,7),(5,6)]]
    products = [{(a[0]*b[0], a[1]*b[1]) for a,b in P} for P in pairs]
    verify(products[0] & products[1] == {(615,49920)}, "stacked singleton alias")
    verify(not ({scalar(a,b) for a,b in pairs[0]} & {scalar(a,b) for a,b in pairs[1]}),
           "crossed-cyclic overlap")
    # Monic polynomials at equal order 11: the difference has degree at most 10.
    # Eleven sample evaluations plus parity are sufficient, and we check 12.
    a, b, c, d = [slots(R,12) for R in [A4,A7,B5,B6]]
    verify(all(a[z]*b[z]-c[z]*d[z] == 2*z**3*(z*z-1)*(z*z-4)*(z*z-2)
               for z in range(13)), "complete polynomial difference")
    L, R = ordinal(A4,A7), ordinal(B5,B6)
    for H in [(615,615),(50535,50535),(2632275,2602035)]:
        verify((sum(ends(L)[-1]),sum(ends(R)[-1])) == H, ("literal depth", H))
        print("literal wrapper observation", len(L), H, flush=True)
        if len(L)<15:
            L, R = triangle(L,(0,)), triangle(R,(0,))

    # Same added vertex, different cyclic chamber, not the ordinal alias.
    left  = ((2,4,9,3), (2,12,9,81,39,15,55))
    right = ((16,5,9,3,14),(2,4,25,19,35,15))
    padded = [(extend(Z,0), W) for Z,W in [left,right]]
    for (Z,W), expected in zip(padded, [(50739,105203,104075),(50739,104523,105411)]):
        tests = [(Z,W),(extend(Z,0),W),(Z,extend(W,0))]
        verify(tuple(sum(ends(triangle(J,K))[-1]) for J,K in tests)==expected,
               "cyclic ambiguity and probes")
        print("three literal probe readings",expected,flush=True)

    physical_profile=(133,773)
    source=(2,8,67,149,39,15,59,119)
    verify(list(rooted(source))[7][:2]==physical_profile
           and sum(ends(fermat_glue(source,7,2))[-1])==3225,
           "same-block eleven-vertex 3225")
    for factor in (59,97):
        verify(any(h==factor for h,q in ps[7])
               and all(h!=factor for n in range(1,7) for h,q in ps[n]),
               ("genuine ordinal factor grade",factor))
    verify(75**2+2*7**2==45**2+2*43**2==5723, "norm of 5723")
    verify(55**2+2*10**2==5**2+2*40**2==3225, "norm of 3225")

    def blowup(A,m):
        return tuple(sum(1<<j for j in range(m*len(A))
                         if (j//m==i//m and i<j) or (A[i//m]>>(j//m)&1))
                     for i in range(m*len(A)))
    verify(sum(ends(blowup((2,4,1),2))[-1])==45,
           "sharp first-level inflation")
    A5={(0,)}
    for n in range(2,6):
        A5={canon(extend(A,c)) for A in A5 for c in range(1<<(n-1))}
    C,D=(2,8,3,21,7),(2,4,9,19,7)
    verify(slots(C,5)==slots(D,5)==[z**5+6*z**3+2*z for z in range(6)],
           "entire old slot-polynomial alias")
    verify([sum(ends(blowup(A,2))[-1]) for A in (C,D)]==[3225,2785],
           "old twin-inflation examples")
    doubled={sum(ends(blowup(A,2))[-1]) for A in A5}
    verify(all(H%4==1 for H in doubled), "level-one congruence")
    for X,Y,a,b in [
        ((2,12,9,49,7,23),(6,12,56,81,35,75,23),25,129),
        ((6,12,48,21,35,11),(2,12,9,113,39,71,23),43,75)]:
        verify(strong(X) and strong(Y)
               and (sum(ends(X)[-1]),sum(ends(Y)[-1]))==(a,b)
               and sum(ends(ordinal(X,Y))[-1])==3225,
               "two same-size factorizations of 3225")

    for N,A,v,j,h,l in [(3225,(2,20,65,7,45,15,59),6,4,73,197),
                         (5723,(6,40,82,21,35,77,155,63),7,2,283,1360)]:
        verify(all(N not in F[n] for n in range(1,len(A))), "minimum strong Fermat source")
        verify(list(rooted(A))[v][:2]==(h,l) and h+(l<<j)==N, "Fermat source profile")
        verify(sum(ends(fermat_glue(A,v,j))[-1])==N, "literal Fermat target")
        print("Fermat target",N,len(A),j,h,l,flush=True)


if __name__ == "__main__":
    main()
