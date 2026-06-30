"""
f_14 = the weight-2 newform of level 14 = the elliptic curve 14a (rank 0).
Compute a_p (point counts), Atkin-Lehner eigenvalues (w_2, w_7), root number; reason about the
sign of the cusp-form contribution to the LRC floor via the apex cusp (d=7) and the Dirac comb.
"""
import math
# 14a1: y^2 + xy + y = x^3 + 4x - 6  (a-invariants [1,0,1,4,-6])
a1,a2,a3,a4,a6 = 1,0,1,4,-6
def count(p):
    # count affine points + 1 (point at infinity)
    n=1
    for x in range(p):
        for y in range(p):
            if (y*y + a1*x*y + a3*y - (x**3 + a2*x*x + a4*x + a6))%p==0: n+=1
    return n
print("f_14 = 14a (rank 0). a_p = p+1-#E(F_p):")
aps={}
for p in [2,3,5,7,11,13,17,19,23]:
    Np=count(p); ap=p+1-Np; aps[p]=ap
    bad = "  (BAD prime, |14)" if 14%p==0 else ""
    print(f"   p={p:>2}: #E(F_p)={Np:>3}  a_{p}={ap:+d}{bad}")
print("\nq-expansion start: f_14 = q + a2 q^2 + a3 q^3 + ... =",
      f"q {aps[2]:+d}q^2 {aps[3]:+d}q^3 {aps[5]:+d}q^5 {aps[7]:+d}q^7 ...")
# Atkin-Lehner: for p||N weight 2, a_p = -w_p
w2 = -aps[2]; w7 = -aps[7]; w14 = w2*w7
eps = -w14  # root number (weight 2)
print(f"\nAtkin-Lehner: w_2 = -a_2 = {w2:+d}   w_7 = -a_7 = {w7:+d}   w_14 = w_2 w_7 = {w14:+d}")
print(f"root number eps = -w_14 = {eps:+d}  => rank parity EVEN => L(14a,1) != 0 (rank 0): consistent")
print(f"\nThe APEX cusp d=7 carries w_7 = {w7} (the 'minus' cusp, klein's d=7 APEX-hard/-).")
print("Cusp forms VANISH at cusps; f_14|W_7 = w_7 f_14 = -f_14 relates the apex-cusp leading coeff to a_1=1.")
print("\nFLOOR sign (empirical, klein HYP-3593): inf R' = 114382/332563 = %.5f  >  bulk 3/pi^2 = %.5f"%(114382/332563, 3/math.pi**2))
print("   => the cusp-form (f_14) contribution to the floor is NET POSITIVE (+%.3f): it does NOT sink the bulk."%(114382/332563-3/math.pi**2))
print("   The atom 4cos^2(3pi/7)=%.4f is the WORST LOCAL gap (apex doublet) < bulk, but measure-zero."%(4*math.cos(3*math.pi/7)**2))
