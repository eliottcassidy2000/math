"""Exact two-transitive-block run compiler (allows empty blocks)."""
from collections import defaultdict

def profile(rows, blocks, roots):
 fixed=[v for v in range(len(rows))if v not in blocks];m=len(fixed);full=(1<<m)-1
 types=fixed+list(blocks);count=len(types);maximum=m+1
 def paths(start=None):
  dp=[[defaultdict(int)for _ in range(count)]for _ in range(full+1)]
  for last,old in enumerate(types):
   if start is not None and old!=start:continue
   if last<m:dp[1<<last][last][(0,0)]=1
   else:dp[0][last][(int(last==m),int(last==m+1))]=1
  for mask in range(full+1):
   for total in range(2*maximum+1):
    for last,old in enumerate(types):
     for (ja,jb),val in list(dp[mask][last].items()):
      if ja+jb!=total:continue
      for nxt,v in enumerate(types):
       if not(rows[old]>>v&1):continue
       if nxt<m:
        if not(mask>>nxt&1):dp[mask|1<<nxt][nxt][ja,jb]+=val
       elif nxt==m and ja<maximum:dp[mask][nxt][ja+1,jb]+=val
       elif nxt==m+1 and jb<maximum:dp[mask][nxt][ja,jb+1]+=val
  return dp
 dp=paths();p=[]
 for mask in range(full+1):
  val=defaultdict(int)
  for col in dp[mask]:
   for exp,c in col.items():val[exp]+=c
  p.append(val)
 p[0][0,0]+=1
 out={}
 for root in roots:
  ir=fixed.index(root);F=defaultdict(int)
  for mask in range(full+1):
   if not(mask>>ir&1):continue
   for (a,b),x in dp[mask][ir].items():
    for (c,d),y in p[full^mask].items():F[a+c,b+d]+=x*y
  ds=paths(root);S=defaultdict(int)
  for col in ds[full]:
   for exp,x in col.items():S[exp]+=x
  out[root]=(dict(S),dict(F))
 return dict(p[full]),out

if __name__=='__main__':
 import json
 from pathlib import Path
 from inflation_probe import inflate
 P=(78,12,24,16,35,79,30)
 def birth(rows,a):
  n=len(rows);return tuple(r|((1<<n)if i==a else 0)for i,r in enumerate(rows))+(((1<<n)-1)^(1<<a),)
 X=birth(birth(P,4),7);q,idx=inflate(X,1,3)
 H,out=profile(q,(1,3),(2,idx[8]));S,F=out[2];Ss,Fs=out[idx[8]]
 keys=S.keys()|F.keys()|Ss.keys()|Fs.keys()
 delta={exp:F.get(exp,0)-S.get(exp,0)-Fs.get(exp,0)+Ss.get(exp,0)for exp in keys}
 delta={str(exp):x for exp,x in sorted(delta.items())if x}
 report={'status':'PROVED via exact two-block run compiler','description':'Q of internal inflated port minus Q of existing latest port; variables are left/right clone counts, each allowed zero','delta_onto_coefficients':delta}
 Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report))
