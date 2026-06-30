"""
The Z_7 cyclotomic SOS floor (the user's hint, made concrete).
THM-515: kernel h(0)=6/7, h(t)=-sin(pi t/7)/(pi t); positive-definite because hat-h(theta)=1_safe(theta)>=0.
TRANSITIVE-GROUP MECHANISM: a G-invariant 2nd moment = a single G-average; positivity in the G-Fourier
basis is an SOS (Bochner/Fejer). On Z_7 (the apex-7 scale) the eigenvalues are hat-h sampled at the 7th
roots theta=j/7: rho_j = 1_safe(j/7). SET-INDEPENDENT (the 7-scale sample ignores the speeds).
"""
import numpy as np
n=14
def safe(theta):  # 1 if ||theta||>=1/14 else 0
    t=theta%1.0; d=min(t,1-t); return 1.0 if d>=1/n-1e-12 else 0.0
print("Z_7 cyclotomic Gram eigenvalues rho_j = 1_safe(j/7), j=0..6 (Bochner positivity on the apex-7 scale):")
rho=[safe(j/7) for j in range(7)]
print(f"   rho = {rho}")
print(f"   => PSD (all rho_j >= 0); floor c = min over the 6 SAFE modes = {min(rho[1:])}; the only zero is")
print(f"      j=0 (DC = the danger point at the origin). SET-INDEPENDENT: depends only on 1/14 & the 7-scale.")
print(f"      ||j/7|| for j=1..6 = {[round(min((j/7)%1,1-(j/7)%1),3) for j in range(1,7)]} all >= 1/7 > 1/14 = safe.")
print()
# kernel h positive-definite: Toeplitz [h(i-j)] PSD
import math
def h(t): return 6/7 if t==0 else -math.sin(math.pi*t/7)/(math.pi*t)
N=20
T=np.array([[h(i-j) for j in range(N)] for i in range(N)])
ev=np.linalg.eigvalsh(T)
print(f"kernel h is positive-definite: Toeplitz [h(i-j)]_{{{N}x{N}}} min eigenvalue = {ev.min():.5f} (>=0 PSD): {ev.min()>=-1e-9}")
print(f"   (Bochner: hat-h(theta)=1_safe(theta)>=0; Fejer-Riesz: 1_safe=|q|^2 => SOS)")
print()
# the DOUBLING-2: danger band is TWO-SIDED (+-1/14), width 2/14=1/7. The '2' = the Z_2 mirror.
print("The doubling-2: danger band (-1/14,1/14) is two-sided, width 2/14 = 1/7. 14 = 2 (mirror Z_2) * 7 (apex Z_7).")
print("CV^2(H) ~ 2/n (metagraph, S_n-transitive group-average) carries the '2' (even-overlap doubling) --")
print("the PROVEN FINITE REHEARSAL of the same transitive-group -> group-average -> Fourier-positivity = SOS.")
print()
# set-independence demo: the j!=0 cyclotomic modes are safe for ANY integer speed (k*j/7, k coprime to 7)
print("set-independence: for any speed v with 7 not| v, ||v*(j/7)|| = ||(vj mod 7)/7|| in {1/7,2/7,3/7} >= 1/7 SAFE")
for v in [1,2,3,5,11,13]:
    vals=[round(min((v*j/7)%1,1-(v*j/7)%1),3) for j in range(1,7)]
    print(f"   v={v:>2}: ||v*j/7|| (j=1..6) = {vals}  all-safe={all(x>=1/7-1e-9 for x in vals)}")
