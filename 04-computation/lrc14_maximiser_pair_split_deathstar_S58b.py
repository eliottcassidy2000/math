import sys
from fractions import Fraction as F
from math import gcd
exec(open('/home/claude/math/04-computation/lrc14_maximiser_sojourn_bound_deathstar_S58b.py').read().split('maxsoj = F(0)')[0])
N=int(sys.argv[1]) if len(sys.argv)>1 else 20
# split by #mirror-pairs: 2 runs (1 pair) vs >=4 runs (>=2 pairs). Track max sojourn in each, and
# confirm: (a) 1/21-run occurs ONLY for (1,2,3)-orbit, (b) multi-pair max sojourn < 2/21.
max1=F(0);d1=None; maxm=F(0);dm=None; run21=[]
for d2 in range(1,N+1):
 for d3 in range(1,N+1):
  for d4 in range(1,N+1):
   if gcd(gcd(d2,d3),d4)!=1: continue
   ivs=bad_intervals([d2,d3,d4]); 
   if not ivs: continue
   soj=sum(hi-lo for lo,hi in ivs)
   if any(hi-lo==F(1,21) for lo,hi in ivs): run21.append((d2,d3,d4))
   if len(ivs)<=2:
     if soj>max1: max1=soj;d1=(d2,d3,d4)
   else:
     if soj>maxm: maxm=soj;dm=(d2,d3,d4)
def is123orbit(d):
    s=sorted(d); return s[1]==2*s[0] and s[2]==3*s[0]
print("components 1..%d, primitive:"%N)
print("  1-mirror-pair (<=2 runs): MAX sojourn = %s=%.6f at %s  (=2/21 & equality is 123-orbit: %s)"%(
    max1,float(max1),d1,max1==F(2,21) and is123orbit(d1)))
print("  >=2-mirror-pairs (>=4 runs): MAX sojourn = %s=%.6f at %s   (< 2/21: %s, gap %.4f)"%(
    maxm,float(maxm),dm,maxm<F(2,21),float(F(2,21)-maxm)))
print("  directions with a 1/21-length run: %d ; ALL are (1,2,3)-orbit: %s"%(
    len(run21), all(is123orbit(d) for d in run21)))
