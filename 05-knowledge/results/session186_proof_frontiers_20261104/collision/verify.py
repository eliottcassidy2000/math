"""Literal independent Hamilton-DP attachment witnesses for S185 perturbations."""
import json
# No inherited counting routines are used.
SEEDS=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3))
def need(x,msg):
 if not x:raise RuntimeError(msg)
def hp(a):
 n=len(a);N=1<<n;d=[[0]*n for _ in range(N)]
 for m in range(1,N):
  for v in range(n):
   if not m>>v&1:continue
   prev=m^(1<<v)
   if not prev:d[m][v]=1
   else:d[m][v]=sum(d[prev][u] for u in range(n) if prev>>u&1 and a[u]>>v&1)
 return sum(d[-1])
def inflate(seed,r):
 a=[0]*(5+r)
 for i in range(5):
  a[i]=seed[i]&31
  for j in range(5,5+r):
   if seed[i]>>5&1:a[i]|=1<<j
   else:a[j]|=1<<i
 for i in range(5,5+r):
  for j in range(i+1,5+r):a[i]|=1<<j
 return a
def attach(t,a,c):
 n=len(t);b=list(t)+[1<<(n+1),0]
 for v in range(n):
  for p,mask in ((n,a),(n+1,c)):
   if mask>>v&1:b[v]|=1<<p
   else:b[p]|=1<<v
 return b
def pairs(t):
 return [(i,j) for i in range(len(t)) for j in range(i+1,len(t)) if not ((t[i]^t[j])&~((1<<i)|(1<<j)))]
def modules(t):
 n=len(t)
 return [m for m in range(1,(1<<n)-1) if m.bit_count()>1 and all((t[w]&m) in (0,m) for w in range(n) if not m>>w&1)]
def strong(t):
 n=len(t)
 for v in range(n):
  seen=1<<v
  while True:
   more=seen
   for u in range(n):
    if seen>>u&1:more|=t[u]
   if more==seen:break
   seen=more
  if seen!=(1<<n)-1:return False
 return True
def run():
 base=inflate(SEEDS[0],6)
 controls=[(592,5,2017,2023,2030,1987,1666,1795,1603,1026,3),(656,5,2017,2023,2030,1987,1667,1794,1603,1026,3),(784,5,2017,2023,2030,1987,1667,1795,1602,1026,3)]
 out=[]
 for t in controls:
  need(not modules(t),'prime');need(strong(t),'strong')
  diffs=[(i,j) for i in range(11) for j in range(i+1,11) if (base[i]^t[i])>>j&1]
  need(len(diffs)==3,'three flips');need(hp(t)==hp(base)==3489,'scalar collision')
  found=None
  for a,c in [(1<<i,0) for i in range(11)]+[(0,1<<i) for i in range(11)]:
   x,y=hp(attach(base,a,c)),hp(attach(t,a,c))
   if x!=y:found={'a_mask':a,'c_mask':c,'base_response':x,'other_response':y};break
  need(found is not None,'literal separating response')
  out.append({'rows':t,'reversed_pairs':diffs,'H':3489,'prime':True,'strong':True,'query':found})
 return {'base_rows':base,'base_homogeneous_pairs':pairs(base),'scalar_hostiles':out,'scope':'Three scalar matches in prescribed finite search, all separated by an actual full-response query.'}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
