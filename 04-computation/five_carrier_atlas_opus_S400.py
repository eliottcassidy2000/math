# opus-2026-07-20-S400 -- HYP-8210: THE FIVE-CARRIER ATLAS from THM-1155.
#
# codex THM-1169 closed the three- and four-carrier r=6 rows.  The engine is my
# THM-1155 density bound  L*(1 - 2k*lam) <= 2*lam*sum(1/w_i), lam = 1/14.
# The first carrier x is safe throughout the owner window [1/(14x), 13/(14x)],
# of length  12/(14x) = 6/(7x).  The LATER carriers must cover it.
#
#   FOUR carriers  -> three later combs, k=3:  L <= (1/4)*sum(1/w_i) < 3/(4x)
#                     window 6/(7x) = 0.857/x  >  0.75/x  => CLOSED immediately.
#   FIVE carriers  -> four later combs,  k=4:  L <= (1/3)*sum(1/w_i) < 4/(3x)
#                     window 0.857/x  <  1.333/x  => NOT closed by the crude bound.
# So the five-carrier row is exactly where the uniform argument first fails, and
# the residual is governed by how much the ACTUAL later carriers (distinct
# integers > x) fall short of the crude 4/x.  This computes that atlas exactly.
from fractions import Fraction as F
from itertools import combinations
LAM=F(1,14)
def budget(k, ws):
    """THM-1155: max coverable length with k combs of moduli ws."""
    return (2*LAM*sum(F(1,w) for w in ws))/(1-2*k*LAM)
print("(1) THE CRUDE COMPARISON, carrier count c (c-1 later combs)")
print("    window = 6/(7x);  crude budget with all later carriers > x")
for c in range(3,8):
    k=c-1
    coef=1-2*k*LAM
    if coef<=0:
        print(f"    c={c}: k={k}, 1-2k*lam = {coef} <= 0 -- density bound DEAD (THM-1155 ceiling)")
        continue
    crude=(2*LAM*k)/coef              # sum(1/w_i) < k/x, so budget < crude/x
    print(f"    c={c}: k={k}  budget < {crude}/x   vs window {F(6,7)}/x   "
          f"=> {'CLOSED' if crude<F(6,7) else 'NOT closed by crude bound'}")

print()
print("(2) THE FIVE-CARRIER ATLAS: exact residual over first carrier x")
print("    later carriers are distinct integers > x; best case for the adversary")
print("    is the four smallest, x+1..x+4.")
print("    x   sum(1/w) max    budget        window 6/(7x)    closed?")
surv=[]
for x in range(1,40):
    ws=[x+1,x+2,x+3,x+4]
    b=budget(4,ws); w=F(6,7*x)
    cl = b < w
    if not cl: surv.append(x)
    if x<=14 or not cl:
        print(f"    {x:3d}  {str(sum(F(1,v) for v in ws)):12s} {float(b):.6f}   {float(w):.6f}      "
              f"{'CLOSED' if cl else 'survives'}")
print()
print(f"    surviving x (density bound alone insufficient): {surv[:20]}{'...' if len(surv)>20 else ''}")
print(f"    total surviving in [1,40): {len(surv)}")

print()
print("(3) WHAT TRUNCATION WOULD CLOSE THE SURVIVORS")
print("    (three-carrier case used truncation at the owner threshold 1/(14p);")
print("     here we ask what window length W would be closed at each x)")
for x in surv[:8]:
    ws=[x+1,x+2,x+3,x+4]; b=budget(4,ws); w=F(6,7*x)
    frac=b/w
    print(f"    x={x:3d}: need window <= {float(b):.6f} = {float(frac):.4f} x the full window"
          f"  -> truncate to {float(frac)*100:.1f}%")
