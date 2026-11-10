# Exact checks for the tilted (biased) tournament ensemble:
#   P_p(T) = p^(forward arcs) q^(back arcs), q = 1-p, reference order 1<...<n.
#   Claim: E_p[H_n] = sum_j gamma_{n,j} (pq)^j   (gamma-vector of the Eulerian polynomial).
from fractions import Fraction as Fr
from math import comb, factorial
import sys
def eulerian(n):
    A = [1]
    for m in range(2, n+1):
        B = [0]*m
        for k in range(m):
            B[k] = (k+1)*(A[k] if k < len(A) else 0) + (m-k)*(A[k-1] if k >= 1 else 0)
        A = B
    return A
def gamma_vec(n):
    A = eulerian(n)[:]  # coefficients of t^k, degree n-1
    g = []
    d = n-1
    for j in range(d//2 + 1):
        c = A[j]; g.append(c)
        # subtract c * t^j (1+t)^(d-2j)
        for i in range(d-2*j+1):
            A[j+i] -= c*comb(d-2*j, i)
    assert all(a == 0 for a in A), A
    return g
def load(n):
    rows = []
    for line in open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "out", f"lab_{n}.txt")):
        b, h, c = map(int, line.split()); rows.append((b, h, c))
    return rows
for n in range(3, 8):
    rows = load(n); m = n*(n-1)//2
    g = gamma_vec(n)
    ok = True
    for p in (Fr(1,2), Fr(3,5), Fr(2,3), Fr(9,10)):
        q = 1-p
        Z = sum(c * q**b * p**(m-b) for b, h, c in rows)
        E1 = sum(c*h * q**b * p**(m-b) for b, h, c in rows)
        E2 = sum(c*h*h * q**b * p**(m-b) for b, h, c in rows)
        rhs = sum(gj*(p*q)**j for j, gj in enumerate(g))
        ok &= (Z == 1 and E1 == rhs)
    print(f"n={n}: gamma={g}  mean formula exact: {ok}")
# coefficient-of-variation and residue structure at n = 7
n = 7; rows = load(n); m = 21
print("\nn=7 tilted ensemble: p, mean, sd, cv, P(H=1 mod 4), mode prob, 2/(sd*sqrt(2pi))")
import math
for p in (0.5, 0.55, 0.6, 0.65, 0.7, 0.8):
    q = 1-p
    pmf = {}
    for b, h, c in rows:
        pmf[h] = pmf.get(h, 0.0) + c * q**b * p**(m-b)
    mean = sum(h*x for h, x in pmf.items()); var = sum(h*h*x for h, x in pmf.items()) - mean**2
    sd = math.sqrt(var)
    p1 = sum(x for h, x in pmf.items() if h % 4 == 1)
    near = max(x for h, x in pmf.items() if abs(h-mean) <= sd)
    low = min(pmf.get(h, 0.0) for h in range(1, 190, 2) if abs(h-mean) <= sd)
    print(f"p={p:.2f} mean={mean:8.3f} sd={sd:7.3f} cv={sd/mean:.3f} P(1 mod 4)={p1:.4f} max/min pmf in mean±sd = {near:.4f}/{low:.5f}  gauss={2/(sd*math.sqrt(2*math.pi)):.4f}")
# uniform measure: residue classes mod 4 and mod 8 for n = 3..7
print("\nuniform measure: P(H = r mod 8), r = 1,3,5,7")
for n in range(3, 8):
    rows = load(n); m = n*(n-1)//2; tot = 2**m
    d = [0]*8
    for b, h, c in rows: d[h % 8] += c
    print(n, [round(d[r]/tot, 4) for r in (1,3,5,7)])
# exact labelled mod-4 discrepancy |P(H = 1 mod 4) - P(H = 3 mod 4)| under the uniform measure
print("\nuniform measure: exact mod-4 discrepancy by order")
for n in range(3, 8):
    rows = load(n); m = n*(n-1)//2
    d = sum(c if h % 4 == 1 else -c for b, h, c in rows)
    print(n, Fr(abs(d), 2**m), float(Fr(abs(d), 2**m)))
