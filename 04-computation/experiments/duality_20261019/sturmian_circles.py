"""General slopes a (irrational, a != log_3 2): balanced rational 3x+1 cycles along the convergents p/q -> a converge to
 subcritical (a < log_3 2):  x_a(theta) = (1/2) sum_{r>=0} 3^r 2^(-floor(r/a + theta)),  theta in [0, 1/a),  positive;
 supercritical (a > log_3 2): x_a(theta) = - sum_{m>=1} 2^(ceil(m/a - theta) - 1) / 3^m,  negative;
 rotation theta -> theta + 1 mod 1/a, letter 1 iff theta >= 1/a - 1 (subcritical form; the supercritical form uses the same word)."""
from mpmath import mp, mpf, floor, ceil, sqrt, log
from math import gcd
mp.dps = 30
def mech(p, q, t): return [(((k+1)*p + t)//q) - ((k*p + t)//q) for k in range(q)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
def convergents(x, n):
    h0,h1,k0,k1 = 0,1,1,0; y=x; out=[]
    for _ in range(n):
        ai=int(floor(y)); h0,h1 = h1, ai*h1+h0; k0,k1 = k1, ai*k1+k0; out.append((h1,k1))
        if y-ai < mpf(10)**-20: break
        y = 1/(y-ai)
    return out
theta_c = log(2)/log(3)
def x_sub(a, th, R=3000): return sum(mpf(3)**r / mpf(2)**int(floor(r/a + th)) for r in range(R)) / 2
def x_sup(a, th, M=3000): return -sum(mpf(2)**(int(ceil(m/a - th)) - 1) / mpf(3)**m for m in range(1, M))
for name, a in (("sqrt2-1 = 0.4142", sqrt(2)-1), ("1/phi = 0.6180", (sqrt(5)-1)/2), ("pi-3+0.5 = 0.6416", mpf(3.14159265358979)-3+mpf('0.5')), ("1/sqrt2 = 0.7071", 1/sqrt(2))):
    sub = a < theta_c
    print(f"slope a = {name} ({'subcritical' if sub else 'supercritical'}); predicted min/max:",
          mp.nstr(x_sub(a, 1/a - mpf(10)**-12), 12) if sub else mp.nstr(x_sup(a, mpf(0)), 12), "/",
          mp.nstr(x_sub(a, mpf(0)), 12) if sub else mp.nstr(x_sup(a, 1/a - mpf(10)**-12), 12))
    for (p, q) in convergents(a, 16):
        if q < 5 or q > 300 or p >= q: continue
        D = 2**q - 3**p; el = [carry(mech(p,q,t)) for t in range(q)]
        print(f"   {p}/{q}: x_min = {mp.nstr(mpf(min(el))/D, 10)}, x_max = {mp.nstr(mpf(max(el))/D, 10)}")
    # conjugacy check on the limit circle
    bad = 0
    import random; random.seed(1)
    for _ in range(40):
        th = mpf(random.random()) / a
        s = 1 if th >= 1/a - 1 else 0
        th2 = (th + 1) % (1/a)
        f = x_sub if sub else x_sup
        y = (3*f(a, th) + 1)/2 if s else f(a, th)/2
        if abs(y - f(a, th2)) > mpf(10)**-20: bad += 1
    print(f"   rotation conjugacy on the limit circle: failures {bad}/40")
