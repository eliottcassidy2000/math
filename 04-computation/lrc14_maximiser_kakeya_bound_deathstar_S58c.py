import sys
from fractions import Fraction as F
from math import gcd
exec(open('/home/claude/math/04-computation/lrc14_maximiser_sojourn_bound_deathstar_S58b.py').read().split('maxsoj = F(0)')[0])
N=int(sys.argv[1]) if len(sys.argv)>1 else 20
# Check: (a) run <= 1/(7D) exactly; (b) the needed inequality 3*#runs <= 2D; (c) total <= 2/21.
bad_run=0; bad_ineq=[]; maxsoj=F(0); worstineq=None; wr=None
for d2 in range(1,N+1):
 for d3 in range(1,N+1):
  for d4 in range(1,N+1):
   if gcd(gcd(d2,d3),d4)!=1: continue
   ivs=bad_intervals([d2,d3,d4])
   if not ivs: continue
   D=max(d2,d3,d4); k=len(ivs)
   tot=sum(hi-lo for lo,hi in ivs)
   if tot>maxsoj: maxsoj=tot
   for lo,hi in ivs:
     if hi-lo>F(1,7*D): bad_run+=1
   # needed: 3k <= 2D
   if 3*k>2*D:
     if worstineq is None or (k,-D)>worstineq[0]: worstineq=((k,-D),(d2,d3,d4),k,D,float(tot))
   # ratio k/D
   if wr is None or F(k,D)>wr[0]: wr=(F(k,D),(d2,d3,d4),k,D)
print("primitive d, components 1..%d:"%N)
print("  runs violating run<=1/(7D): %d  (should be 0 -> the Kakeya box-exit bound is exact)"%bad_run)
print("  MAX ratio #runs/D = %s at %s (#runs=%d, D=%d)"%(wr[0],wr[1],wr[2],wr[3]))
print("  needed inequality 3*#runs <= 2D:  %s"%("HOLDS everywhere" if worstineq is None else "FAILS e.g. %s (#runs=%d,D=%d,total=%.4f)"%(worstineq[1],worstineq[2],worstineq[3],worstineq[4])))
print("  overall max sojourn = %s = %.6f (2/21=%.6f)"%(maxsoj,float(maxsoj),2/21))
