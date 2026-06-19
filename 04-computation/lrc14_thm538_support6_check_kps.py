import sys, itertools, cmath
from cmath import exp, pi, sin
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
# chat_T(n) = Fourier coeff of 1_{[0,1) minus union_{j in T} sector_j}, sector_j=[j/7,(j+1)/7).
# n=0: 1 - |T|/7. n!=0: -sum_{j in T} e^{-2pi i n (j+0.5)/7} * sin(pi n/7)/(pi n).
def chat(T, n):
    if n == 0:
        return complex(1 - len(T)/7.0, 0)
    s = sin(pi*n/7)/(pi*n)
    tot = 0+0j
    for j in T:
        tot += cmath.exp(-2j*pi*n*(j+0.5)/7) * s
    return -tot
def K(nvec):
    tot = 0+0j
    for r in range(0,7):
        for T in itertools.combinations(range(1,7), r):
            prod = 1+0j
            for nj in nvec:
                prod *= chat(set(T), nj)
            tot += ((-1)**r) * prod
    return tot
print("=== THM-538 support-6 floor check: is K(n)=0 for support < 6 ? ===")
tests = [
  ("support2 (1,-1,0,0,0,0,0)", (1,-1,0,0,0,0,0)),
  ("support2 (1,1,0,0,0,0,0)",  (1,1,0,0,0,0,0)),
  ("support3 (1,1,-1,0,0,0,0)", (1,1,-1,0,0,0,0)),
  ("support3 (1,-1,1,0,0,0,0)", (1,-1,1,0,0,0,0)),
  ("support4 (1,1,-1,-1,0,0,0)",(1,1,-1,-1,0,0,0)),
  ("support5 (1,1,1,-1,-1,0,0)",(1,1,1,-1,-1,0,0)),
  ("support6 (1,1,1,-1,-1,-1,0)",(1,1,1,-1,-1,-1,0)),
  ("support7 (1,1,1,-1,-1,-1,1)",(1,1,1,-1,-1,-1,1)),
]
for name, n in tests:
    k = K(n)
    print(f"  {name}: K = {k.real:+.6f} {k.imag:+.6f}i   |K|={abs(k):.6f}  {'~0 (floor holds)' if abs(k)<1e-9 else '*** NONZERO -> floor FAILS for K'}")
