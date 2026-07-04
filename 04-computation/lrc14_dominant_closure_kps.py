"""kps-S5: TEST the dominant-closure argument. Claim: if V has a runner v_i > 13*max(others),
then V is lonely (M>=1/14). Proof idea: base U=V\{v_i} (12 runners) has M(U)>=1/13 (LRC13) at some t*;
f=min_j||v_j t|| is B-Lipschitz (B=max others), so f>=1/14 on |t-t*|<=1/(182B) (width W=1/(91B));
v_i-bad intervals have width 1/(7|v_i|) < W iff |v_i|>13B, so the base-safe interval hits a v_i-safe point."""
from fractions import Fraction
def resdist_frac(vt):  # ||vt|| for vt a Fraction
    f = vt - int(vt)
    if f < 0: f += 1
    return min(f, 1-f)
def M_of(V, Qmax=4000):
    from math import gcd
    best=Fraction(0)
    for q in range(2,Qmax+1):
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            md=min(min((v*p)%q, q-(v*p)%q) for v in V)
            if Fraction(md,q)>best: best=Fraction(md,q)
    return best
import random
random.seed(1)
print("Testing dominant closure: V = base(12 runners, all <= B) + v_i > 13*B  =>  M(V) > 1/14 ?")
print(f"1/14 = {float(1/14):.5f}")
fails=0; tested=0
for trial in range(40):
    # random base of 12 distinct nonzero speeds up to some B
    B = random.choice([13, 20, 30, 13, 13])
    base = random.sample(range(1, B+1), 12)
    Bmax = max(base)
    vi = random.randint(13*Bmax+1, 13*Bmax + 200)  # dominant: vi > 13*Bmax
    V = base + [vi]
    m = M_of(V, Qmax=min(6*vi, 6000))
    tested += 1
    ok = m > Fraction(1,14)
    if not ok:
        fails += 1
        print(f"  FAIL: base={sorted(base)} Bmax={Bmax} vi={vi} (>{13*Bmax}) M={m}={float(m):.5f}")
    if trial < 6:
        print(f"  base(max {Bmax}) + vi={vi}: M={m}={float(m):.5f}  >1/14: {ok}  [vi/Bmax={vi/Bmax:.1f}]")
print(f"\n{tested} dominant families tested, {fails} failures (M<=1/14).")
# also test the SHARP boundary: vi just above vs at 13*Bmax
print("\nBoundary check (base={1..12}, Bmax=12, 13*Bmax=156):")
for vi in [150, 156, 157, 170, 182]:
    V=list(range(1,13))+[vi]; m=M_of(V,Qmax=min(6*vi,4000))
    print(f"  vi={vi} (dominant iff >156): M={m}={float(m):.5f} >1/14:{m>Fraction(1,14)}")
