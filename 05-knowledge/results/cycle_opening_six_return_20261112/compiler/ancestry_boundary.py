"""Exact structural controls for the restricted-birth ancestry boundary."""
from pathlib import Path
import json
import gzip
from collections import Counter,defaultdict
from hashlib import sha256
P=(78,12,24,16,35,79,30)
Q=(250,508,505,496,224,192,128,256,625,255)
def need(ok,why):
 if not ok:raise RuntimeError(why)
def indegrees(rows):return [sum(r>>j&1 for r in rows)for j in range(len(rows))]
def birth(rows,p):
 n=len(rows);return tuple(r|((1<<n)if i==p else 0)for i,r in enumerate(rows))+(((1<<n)-1)^(1<<p),)
def substitute(rows,e,module):
 owner=[i for i in range(len(rows))for _ in range(len(module)if i==e else 1)]
 local=[j for i in range(len(rows))for j in range(len(module)if i==e else 1)]
 out=tuple(sum(((module[local[x]]>>local[y]&1)if u==v else(rows[u]>>v&1))<<y for y,v in enumerate(owner)if x!=y)for x,u in enumerate(owner))
 return out,{i:owner.index(i)for i in range(len(rows))}
def strong(rows):
 for start in range(len(rows)):
  seen=1<<start
  while True:
   nxt=seen
   for i,r in enumerate(rows):
    if seen>>i&1:nxt|=r
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(rows))-1:return False
 return True
def remove(rows,v):
 old=[i for i in range(len(rows))if i!=v]
 return tuple(sum((rows[i]>>j&1)<<k for k,j in enumerate(old))for i in old)
def main():
 report=[]
 modules=[('TT1',(0,)),('TT2',(2,0)),('C3',(2,4,1)),('C4',(2,4,9,3)),('H9_five_core',(22,4,8,3,14)),('weak_C3_then_singleton',(10,12,9,0))]
 for name,M in modules:
  m=len(M);p,idx=substitute(P,1,M);s=birth(p,idx[4]);x=birth(s,len(s)-1);raw=birth(x,len(x)-1);target,qi=substitute(Q,3,M)
  ds,dr,dq=map(indegrees,(s,raw,target));z=len(s)-1;xr=len(raw)-2;yr=len(raw)-1
  need([i for i,d in enumerate(ds)if d==1]==[z]and ds[idx[4]]==3,'source unique root and degree-three predecessor')
  need([i for i,d in enumerate(dr)if d==1]==[xr,yr]and raw[xr]>>yr&1,'raw ordered two-root ancestry')
  need([i for i,d in enumerate(dq)if d==1]==[qi[9]]and dq[qi[8]]==m+3,'returned unique root and predecessor degree')
  expectedS={0:3,2:m+4,3:m+5,4:3,5:2,6:3}
  need(all(ds[idx[v]]==d for v,d in expectedS.items()),'complete source exterior degrees')
  need(all(ds[idx[1]+j]==5+d for j,d in enumerate(indegrees(M))),'source module degrees')
  expectedQ={0:3,1:2,2:2,4:m+5,5:m+6,6:m+7,7:m+7,8:m+3,9:1}
  need(all(dq[qi[v]]==d for v,d in expectedQ.items()),'complete returned exterior degrees')
  need(all(dq[qi[3]+j]==4+d for j,d in enumerate(indegrees(M))),'returned module degrees')
  need(all(strong(T)for T in(p,s,x,raw,target,remove(target,qi[9]))),'strong outputs and adaptive parents')
  high=birth(p,idx[0]);dh=indegrees(high)
  need([i for i,d in enumerate(dh)if d==1]==[len(high)-1]and dh[idx[0]]==2,'high side also strictly primitive')
  report.append(dict(module=name,size=m,source_indegrees=ds,raw_indegrees=dr,returned_indegrees=dq,source_root=z,source_port=idx[4],returned_root=qi[9],returned_port=qi[8]))
 bank=Path(__file__).resolve().parents[2]/'metagraph_aspect_sequences_20261111/records.jsonl.gz'
 census=defaultdict(Counter);ordinary_rows=0;inverse_checks=0
 with gzip.open(bank,'rt',encoding='utf-8')as fh:
  for line in fh:
   rec=json.loads(line);ordinary_rows+=1
   if rec['n']<3 or not rec['strong']:continue
   rows=tuple(rec['rows']);ones=[v for v,d in enumerate(indegrees(rows))if d==1];census[len(rows)][len(ones)]+=1
   need(len(ones)<=2 or len(rows)==3,'at most two roots outside triangle')
   if len(ones)==2:
    a,b=ones
    if not(rows[a]>>b&1):a,b=b,a
    parent=remove(rows,b);port=a-int(a>b)
    need(strong(parent)and indegrees(parent)[port]==1,'strong two-root strict inverse')
    inverse_checks+=1
 need(census[4][2]==1,'triangle first strict descendant')
 for n in range(5,max(census)+1):need(census[n][2]==census[n-1][1]+2*census[n-1][2],('ordinary strict recurrence',n))
 out={'status':'VERIFIED exact structural controls; universal statements proved in ANCESTRY_BOUNDARY.md','controls':report,'bank_sha256':sha256(bank.read_bytes()).hexdigest(),'ordinary_bank_rows':ordinary_rows,'bank_scope':'complete ordinary tournament bank through order8, strong members only for recurrence','root_count_census':{str(n):dict(sorted(v.items()))for n,v in sorted(census.items())},'strict_inverse_checks':inverse_checks}
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print('PASS six module degree profiles and complete order8 strong bank: strict inverse and ordinary recurrence; '+str(inverse_checks)+' inverse controls')
if __name__=='__main__':main()
