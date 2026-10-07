# Verify #090 node sets, sine-product coefficients, node constants; list artificial nodes.
import numpy as np, mpmath as mp
from math import gcd
mp.mp.dps = 30
LIM = 5000
losch = set()
for j in range(-80, 81):
    for l in range(-80, 81):
        v = j*j + j*l + l*l
        if v <= LIM: losch.add(v)
def is_losch(n):
    # n = x^2+xy+y^2 iff every prime p = 2 mod 3 has even exponent
    m = n; p = 2
    if n == 0: return True
    while p*p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0: m//=p; e+=1
            if p % 3 == 2 and e % 2: return False
        p += 1
    if m > 1 and m % 3 == 2: return False
    return True
assert all(is_losch(n) == (n in losch) for n in range(LIM+1))
for mod in (12, 36, 4*9*25, 36*7):
    res = sorted({n % mod for n in losch})
    print(f"Loeschian residues mod {mod}: count {len(res)} density {len(res)/mod:.4f}", res if mod <= 36 else "")
L = [0,1,3,4,7,9]
I = [0,1,3,4,7,9,12,13,16,19,21,25,27,28,31]
assert sorted({n % 12 for n in losch}) == L
assert sorted({n % 36 for n in losch}) == I
# artificial nodes (in superset, not Loeschian)
art12 = [n for n in range(1, 400) if n % 12 in L and not is_losch(n)]
art36 = [n for n in range(1, 400) if n % 36 in I and not is_losch(n)]
print("artificial nodes mod-12 superset (<400):", art12[:25])
print("artificial nodes mod-36 superset (<400):", art36[:25])
# fraction of artificial among superset nodes up to X
for X in (168, 1000, 5000):
    s12 = [n for n in range(1, X+1) if n % 12 in L]; a12 = [n for n in s12 if not is_losch(n)]
    s36 = [n for n in range(1, X+1) if n % 36 in I]; a36 = [n for n in s36 if not is_losch(n)]
    print(f"X={X}: superset12 {len(s12)} artificial {len(a12)} | superset36 {len(s36)} artificial {len(a36)}")
# sine product coefficients P_l, P(s)=prod (2 sin(pi(s-a)/12))^2 = sum_{l=-6}^6 P_l e^{i pi l s/6}
w = np.poly1d([1.0])
for a in L:
    r = np.exp(1j*np.pi*a/6)
    w = w * np.poly1d([1.0, -r]) * np.poly1d([1.0, -r])
coef = w.coeffs[::-1]  # coefficient of w^k, k=0..12 ; P = w^{-6} * prod
P = {k-6: coef[k] for k in range(13)}
b = np.sqrt(3)/2
claimed = [5, -1+2j*b, 1+2j*b, -2, -0.5+1j*b, -1-2j*b, 1]
print("P_l (l=0..6):", [complex(round(P[l].real,12), round(P[l].imag,12)) for l in range(7)])
print("matches paper:", all(abs(P[l]-claimed[l]) < 1e-9 for l in range(7)))
# Q_a, D_a
for a in L:
    Q = (mp.pi/6)**2 * mp.fprod([(2*mp.sin(mp.pi*(a-d)/12))**2 for d in L if d != a])
    D = mp.pi/6 * mp.fsum([mp.cot(mp.pi*(a-d)/12) for d in L if d != a])
    print(f"a={a}: Q={mp.nstr(Q,12)}  D={mp.nstr(D,12)}")
print("pi^2/3 =", mp.nstr(mp.pi**2/3,12), " (2pi^2/3)(2-sqrt3) =", mp.nstr(2*mp.pi**2/3*(2-mp.sqrt(3)),12),
      " (2pi^2/3)(2+sqrt3) =", mp.nstr(2*mp.pi**2/3*(2+mp.sqrt(3)),12))
print("7pi sqrt3/18 =", mp.nstr(7*mp.pi*mp.sqrt(3)/18,12), " pi(3+sqrt3)/18 =", mp.nstr(mp.pi*(3+mp.sqrt(3))/18,12), " pi(3-sqrt3)/18 =", mp.nstr(mp.pi*(3-mp.sqrt(3))/18,12))
# r(D) for triangular lattice: 6 * sum_{d|D} chi_{-3}(d)
def chi3(d): return 0 if d % 3 == 0 else (1 if d % 3 == 1 else -1)
def r(D): return 6*sum(chi3(d) for d in range(1, D+1) if D % d == 0)
print("r(D) D=1..31:", [(D, r(D)) for D in range(1, 32) if r(D)])
print("r(183) =", r(183), " 183 mod 12 =", 183 % 12, " mod 36 =", 183 % 36, " 7,21,13 in L-superset:", [x % 12 in L for x in (7,21,13)])
