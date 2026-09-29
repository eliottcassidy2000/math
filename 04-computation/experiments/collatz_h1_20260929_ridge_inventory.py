"""Ridge inventory of the negative family by the 3-adic depth of the seed coincidence.
A ridge of the negative family {2^-m mod 3^n} is born where 2^-Q = -+u mod 3^k with u a small odd multiplier:
depth k(u, sign, Q) = v_3(u 2^Q -+ 1).  For fixed (u, sign) the exponents Q with depth >= k form one residue class
mod L_k = 2 3^(k-1) (2 is a primitive root mod 3^k), so the seeds of depth >= k in a window of length W number
about W/L_k per (u, sign); LTE gives v_3(2^Q - 1) = 1 + v_3(Q) (Q even), v_3(2^Q + 1) = 1 + v_3(Q) (Q odd) for u = 1.
Predicted birth amplitude (S22 multiplier law): M(k) u^(-0.44) at level k, i.e. normalised v ~ 3^(k/2) M(k) u^(-0.44),
and the ridge then lives on the line m = Q - resonant exponent(k, u) - 1.3 (n - k).
Output: all seeds with depth >= kmin, u <= umax odd (3 not dividing u), Q in [Qmin, Qmax]; counts by depth."""
import sys, math
from collections import Counter
umax = int(sys.argv[1]) if len(sys.argv) > 1 else 200
Qmax = int(sys.argv[2]) if len(sys.argv) > 2 else 700
kmin = int(sys.argv[3]) if len(sys.argv) > 3 else 5
LOG23 = math.log2(3)
# M(k) from the S22 profile (FFT to 18, recursion beyond): use the no-descent proxy 0.46 P_k? we only need order of magnitude
# use the level values M(h) h<=18 from the note
Mh = {1:.5774,2:.3779,3:.2522,4:.1770,5:.1293,6:.0961,7:.0759,8:.0609,9:.0480,10:.0383,11:.0319,12:.0265,13:.0221,14:.0191,15:.0163,16:.0144,17:.0125,18:.0112,20:8.88e-3,30:3.29e-3,40:1.38e-3}
def v3(x):
    if x == 0: return 99
    c = 0
    while x % 3 == 0: x //= 3; c += 1
    return c
seeds = []
for u in range(1, umax+1, 2):
    if u % 3 == 0: continue
    p = 1  # 2^Q
    for Q in range(0, Qmax+1):
        for sign, lab in ((-1, '-'), (1, '+')):     # u 2^Q + sign*1 ... depth of u 2^Q -+ 1: sign -1: u2^Q - 1 ; +1: u2^Q + 1
            d = v3(u*p + sign)
            if d >= kmin: seeds.append((d, u, lab if sign == 1 else '-', Q))
        p <<= 1
seeds.sort(reverse=True)
print(f"seeds 2^-Q = -+u mod 3^k with depth k = v_3(u 2^Q -+ 1) >= {kmin}, u <= {umax} (odd, 3 not | u), 0 <= Q <= {Qmax}: {len(seeds)} seeds")
print("depth  u  sign  Q   [meaning: u*2^Q + (sign)1 = 0 mod 3^depth, i.e. 2^-Q = -(sign) u mod 3^depth]   predicted birth v ~ 3^(k/2) M(k) u^-0.44")
for (d, u, lab, Q) in seeds[:60]:
    Mk = Mh.get(d, 0.46*0.9465**d*d**-1.3*9)  # rough
    print(f"  {d:2d}  {u:4d}   {lab}  {Q:4d}    v_birth ~ {3**(d/2)*Mk*u**-0.44:.2f}")
cnt = Counter(d for (d, *_ ) in seeds)
print("count by depth:", sorted(cnt.items(), reverse=True))
print("expected count per depth k over all (u,sign) pairs: (#pairs) * (Qmax+1) / (2*3^(k-1)) * (1 - 1/3) [exact depth]:")
npairs = 2*sum(1 for u in range(1, umax+1, 2) if u % 3)
for k in range(kmin, 17):
    print(f"   k={k}: expected {npairs*(Qmax+1)/(2*3**(k-1))*(2/3):.2f}, observed {cnt.get(k,0)}")
print("the seeds named by S22: (u=55,+,Q=423) depth", v3(55*2**423+1), "; (u=13,-,Q=154) depth", v3(13*2**154-1), "; (u=1,-,Q=162) depth", v3(2**162-1), "; (u=1,-,Q=486) depth", v3(2**486-1))
# top seeds by predicted birth amplitude within the S22 window Q <= 600
print("\ntop 15 seeds by predicted birth amplitude, Q <= 600, u <= 200:")
amp = []
for (d, u, lab, Q) in seeds:
    if Q <= 600:
        Mk = Mh.get(d, 0.46*0.9465**d*d**-1.3*9)
        amp.append((3**(d/2)*Mk*u**-0.44, d, u, lab, Q))
amp.sort(reverse=True)
for a in amp[:15]: print(f"  v~{a[0]:6.2f}  depth {a[1]:2d}  u={a[2]:4d} {a[3]}  Q={a[4]}")
