from fractions import Fraction as Fr
from math import factorial
import math
def load(n):
    return [tuple(map(int, l.split())) for l in open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "out", f"lab_{n}.txt"))]
def polyfit_w(n, k):
    """coefficients c_j with E_p[H^k] = sum_j c_j w^j, w = p(1-p), by exact interpolation"""
    rows = load(n); m = n*(n-1)//2; deg = m//2
    pts = []
    for i in range(deg+1):
        p = Fr(1,2) + Fr(i+1, 2*(deg+2))      # distinct w values
        q = 1-p
        val = sum(c * h**k * q**b * p**(m-b) for b, h, c in rows)
        pts.append((p*q, val))
    # Newton / Lagrange exact
    xs = [x for x, _ in pts]; ys = [y for _, y in pts]
    coef = [Fr(0)]*(deg+1)
    for i in range(deg+1):
        # Lagrange basis polynomial
        num = [Fr(1)]; den = Fr(1)
        for j in range(deg+1):
            if j == i: continue
            num = [Fr(0)] + num
            for t in range(len(num)-1):
                num[t] -= xs[j]*num[t+1]
            den *= (xs[i]-xs[j])
        for t in range(len(num)):
            coef[t] += ys[i]*num[t]/den
    return coef
for n in range(3, 8):
    c2 = polyfit_w(n, 2)
    while c2 and c2[-1] == 0: c2.pop()
    print(f"n={n}: E_p[H^2] in powers of w:", [str(x) for x in c2])
# tilt needed to bring the mean of order n down to the uniform mean of order n-1, and cv there
def gamma_vec(n):
    from math import comb
    A = [1]
    for mm in range(2, n+1):
        B = [0]*mm
        for k in range(mm):
            B[k] = (k+1)*(A[k] if k < len(A) else 0) + (mm-k)*(A[k-1] if k >= 1 else 0)
        A = B
    g = []; d = n-1
    for j in range(d//2+1):
        c = A[j]; g.append(c)
        for i in range(d-2*j+1): A[j+i] -= c*comb(d-2*j, i)
    return g
print("\n n : w_n (mean_n(w_n) = mu_{n-1}),  p_n,  (p_n-1/2)*sqrt(n/ln n)")
for n in list(range(4, 13)) + [16, 24, 40, 80, 160]:
    g = gamma_vec(n); target = factorial(n-1)/2**(n-2)
    lo, hi = 0.0, 0.25
    for _ in range(200):
        mid = (lo+hi)/2
        if sum(gj*mid**j for j, gj in enumerate(g)) > target: hi = mid
        else: lo = mid
    w = lo; p = 0.5 + math.sqrt(0.25 - w)
    print(n, round(w, 6), round(p, 5), round((p-0.5)*math.sqrt(n/math.log(n)), 4))
