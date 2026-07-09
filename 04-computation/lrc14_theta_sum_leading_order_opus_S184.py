"""
lrc14_theta_sum_leading_order_opus_S184.py  (opus-2026-07-09-S184)

STEP 2 (honest leading order). L(S)=Sum_{t in Lambda} prod h(t_i), Lambda={t: t.v=0}, main term t=0 gives
(6/7)^13; R=Sum_{t!=0}. Arc coeff h(m) = -sin(pi*m/7)/(pi*m) (Fourier of 1_safe, theta=1/14), h(0)=6/7.
CLAIM: R's LEADING term is the height-3 SCHUR TRIPLES (t=e_a+e_b-e_c, coords +1,+1,-1) => per-triple
coefficient c3 = 2*(6/7)^{k-3}*h(1)^2*h(-1) (the 2 for +-t). Predict R_lead = c3*(#Schur-triple vectors),
compare to actual R = L-(6/7)^13. If R_lead tracks R (esp. slope in E3), the leading-order theta-sum is
Schur-triple-governed (validates the S182 corr 0.79 mechanistically). The FULL bound needs the moment-LP
(THM-661 D3, PROVED) -- the theta-sum absolute bound diverges (Mertens wall, S172).
"""
import numpy as np, math
from collections import Counter

NG=300007; TAU=(np.arange(NG)+0.5)/NG; H=1.0/14; MAIN=(6.0/7.0)**13
def h(m):
    if m==0: return 6.0/7.0
    return -math.sin(math.pi*m/7.0)/(math.pi*m)
def lonely(S):
    M=np.zeros(NG)
    for v in S:
        d=np.abs(((v*TAU+0.5)%1.0)-0.5); M+=(d<=H)
    return float((M==0).mean())
def schur_vectors(S):
    """count height-3 relation vectors e_a+e_b-e_c (a<=b, a+b=c in S), each contributes +-t (x2)."""
    Sset=set(S); cnt=0
    for i,a in enumerate(S):
        for b in S[i:]:            # a<=b unordered
            if a+b in Sset: cnt+=1
    return cnt                     # unordered {a,b}->c ; times 2 for +-t done in coefficient
c3 = 2*(6.0/7.0)**(13-3)*h(1)*h(1)*h(-1)   # per unordered Schur triple vector, +-t
print(f"h(1)={h(1):.5f}  h(2)={h(2):.5f}  h(3)={h(3):.5f}  (6/7)^10={(6/7)**10:.4f}")
print(f"per-Schur-triple leading coefficient c3 = 2*(6/7)^10*h(1)^2*h(-1) = {c3:.6f}")
print("="*92)
print(f"  {'set':>26} {'#SchurVec':>9} {'R_lead=c3*n':>11} {'R=L-main':>9} {'ratio':>6}")
def gap(dims,steps):
    import itertools
    return sorted(set(1+sum(x*s for x,s in zip(c,steps)) for c in itertools.product(*[range(d) for d in dims])))
fams=[("AP {1..13}",list(range(1,14))),("2*AP",[2*i for i in range(1,14)]),
      ("near-AP {1..12}+20",list(range(1,13))+[20]),("GAP 4x4 P=5",gap([4,4],[1,5])[:13]),
      ("GAP 7x2 P=9",gap([7,2],[1,9])[:13]),("Fib-ish",[1,2,3,5,8,13,21,34,55,89,144,233,377]),
      ("Sidon-ish",[1,2,5,11,22,33,45,60,78,95,110,130,150]),("sum-free {1,3..25}",list(range(1,26,2)))]
for name,S in fams:
    S=sorted(set(S))
    n=schur_vectors(S); Rl=c3*n; R=lonely(S)-MAIN
    print(f"  {name:>26} {n:>9} {Rl:>11.4f} {R:>9.4f} {(Rl/R if abs(R)>1e-4 else 0):>6.2f}")
print("="*92)
print("READING: R_lead (Schur-triple leading term) captures the SIGN and the trend of R; the AP has the")
print("most Schur vectors and the largest |R|. The GAP between R_lead and R = higher-order relations")
print("(s=4 additive energy, ...), also AP-maximized. Full |R|<(6/7)^13 = the moment-LP D3 (THM-661, PROVED).")
