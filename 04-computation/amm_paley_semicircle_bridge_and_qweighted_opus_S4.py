"""
The AMM<->Paley semicircle bridge, and the new q-weighted Paley tournament (exact closed forms).

BRIDGE: AMM golden floor = Catalan GF at branch point w=-1 (THM-3009); Paley signed cluster GF (THM-438 (**))
is F(x)=sum(-1)^k C_k x^k = C(-x), loop equation x F^2 + F - 1 = 0. So golden = F(1) = C(-1) = 1/phi is the
EDGE value of the same free-semicircle object whose exponential resummation is the path-ratio R -> e.

NEW: q-weighted Paley  R_q(p)=E_sigma[prod_k (1+q chi(d_k))]:
   path-ratio   R_q(p) -> e^{q^2}                       (only the cherry L=2 cluster survives; a_2=1)
   cluster edge F_q(1) = (sqrt(1+4q^2)-1)/(2q^2)         (F_q(x)=F(q^2 x); q^2 x F^2 + F - 1 = 0)
q=1 recovers R->e and F_1(1)=1/phi (AMM golden).
"""
import mpmath as mp
from itertools import permutations
mp.mp.dps = 20
phi = (1+mp.sqrt(5))/2

# --- bridge: golden solves the Paley loop equation x F^2 + F - 1 = 0 at x=1 ---
F1 = 1/phi
print("BRIDGE: Paley loop eq  1*F^2 + F - 1  at F = 1/phi = C(-1):", mp.nstr(F1*F1+F1-1, 4), "(=0)")
print("        => AMM golden edge C(-1)=F(1)=1/phi is the edge of the Paley semicircle GF.\n")

# --- q-weighted Paley path-ratio, brute force over Paley tournament on F_p ---
def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p-1)//2, p) == 1 else -1)
def Rq(p, q):
    tot = mp.mpf(0); cnt = 0
    for perm in permutations(range(1, p)):
        seq = (0,)+perm; prod = mp.mpf(1)
        for k in range(p-1):
            prod *= (1+q*legendre(seq[k+1]-seq[k], p))
        tot += prod; cnt += 1
    return tot/cnt
print("q-weighted path-ratio R_q(7) vs limit e^{q^2}:")
for q in [mp.mpf('0.5'), mp.mpf('1')]:
    print("  q=%.2f: R_q(7)=%s  ->  e^{q^2}=%s" % (float(q), mp.nstr(Rq(7, q), 7), mp.nstr(mp.e**(q*q), 7)))

# --- cluster edge closed form ---
def Fq1(q): return (mp.sqrt(1+4*q*q)-1)/(2*q*q)
print("\ncluster edge F_q(1)=(sqrt(1+4q^2)-1)/(2q^2)  [solves q^2 F^2+F-1=0]:")
for q in [mp.mpf('0.5'), mp.mpf('1'), mp.mpf('2')]:
    f = Fq1(q)
    print("  q=%.2f: F_q(1)=%s   loop residual=%s" % (float(q), mp.nstr(f, 12), mp.nstr(q*q*f*f+f-1, 3)))
print("  q=1 => F_1(1)=%s == 1/phi (AMM golden)" % mp.nstr(Fq1(mp.mpf(1)), 15))
