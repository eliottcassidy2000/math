#!/usr/bin/env python3
"""Independent all-factor-pairs K4,4 exclusion; bounded positive integer universe."""
import math,itertools,json,sys,time
from checker import require,common_factor

def run(N):
 require(type(N) is int and N>=1,'invalid N')
 start=time.monotonic();seen=set();eligible=0;factor_capable_deltas=0
 for d in range(1,N):
  bases=[]
  for p in range(1,math.isqrt(d)+1):
   if d%p:continue
   q=d//p
   if p<q and (q-p)%2==0:bases.append(((q-p)//2)**2)
  if len(bases)<4:continue
  factor_capable_deltas+=1
  # Every a is covered unless the factor bound already makes four neighbors impossible.
  for a in range(1,N-d+1):
   b=a+d;C=sorted(t-a for t in bases if 1<=t-a<=N and t-a not in (a,b))
   if len(C)<4:continue
   eligible+=1
   for R in itertools.combinations(C,4):
    if R in seen:continue
    seen.add(R)
    L=[x for x in common_factor(R[0],R[1],N) if x not in R and math.isqrt(x+R[2])**2==x+R[2] and math.isqrt(x+R[3])**2==x+R[3]]
    require(len(L)<4,f'K44 found: L={L[:4]}, R={R}')
 return {'status':'PASS','N':N,'claim':'No ordinary K4,4 in G_N','factor_capable_deltas':factor_capable_deltas,'eligible_seed_pairs':eligible,'distinct_right_quadruples':len(seen),'seconds':time.monotonic()-start}
if __name__=='__main__':
 try:print(json.dumps(run(int(sys.argv[1]) if len(sys.argv)>1 else 10000),indent=2))
 except Exception as e:print(json.dumps({'status':'FAIL','type':type(e).__name__,'message':str(e)}),file=sys.stderr);sys.exit(1)
