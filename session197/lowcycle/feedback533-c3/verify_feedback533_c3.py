from pathlib import Path
from itertools import combinations
import hashlib,json
P=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
a=[2,268,905,881,999,967,903,779,513,3];n=len(a);full=(1<<n)-1
need(all(0<=r<=full and not(r>>i&1)for i,r in enumerate(a)),'row range')
need(all(((a[i]>>j&1)+(a[j]>>i&1))==1 for i in range(n)for j in range(i+1,n)),'tournament')
d=[[0]*n for _ in range(1<<n)];c=[[0]*n for _ in range(1<<n)]
for v in range(n):d[1<<v][v]=1
c[1][0]=1
for S in range(1,1<<n):
 for u in range(n):
  q=a[u]&~S
  while q:
   b=q&-q;q-=b;v=b.bit_length()-1;d[S|b][v]+=d[S][u];c[S|b][v]+=c[S][u]
H=sum(d[full]);C=sum(c[full][v]for v in range(n)if a[v]&1)
paths=0;cycles=[]
def rec(w,S):
 global paths
 if S==full:
  paths+=1
  if w[0]==0 and a[w[-1]]&1:cycles.append(w)
  return
 q=a[w[-1]]&~S
 while q:
  b=q&-q;q-=b;rec(w+[b.bit_length()-1],S|b)
for v in range(n):rec([v],1<<v)
expected=[[0,1,2,3,4,5,6,7,8,9],[0,1,2,7,3,4,5,6,8,9],[0,1,3,4,5,6,2,7,8,9]]
need(H==paths==533 and C==len(cycles)==3 and cycles==expected,'independent exact H/C')
ports=[1,3];tail=[4,5,6,2,7,8,9,0];masks=[1,1,1,2,3,0,1,1]
need(a[1]>>3&1,'initial edge');birth=ports+tail;arcs=1
for k,(v,mask)in enumerate(zip(tail,masks),2):
 for u in birth[:k]:
  desired=u in ports and bool(mask>>ports.index(u)&1)
  need(bool(a[v]>>u&1)==desired,'fixed port arc');arcs+=1
need(arcs==45 and all(a[u]>>v&1 for i,u in enumerate(tail)for v in tail[i+1:]),'upper feedback witness')
triangles=[]
for deleted in range(n):
 found=None
 for x,y,z in combinations([v for v in range(n)if v!=deleted],3):
  for t in ([x,y,z],[x,z,y]):
   if all(a[t[i]]>>t[(i+1)%3]&1 for i in range(3)):found=t;break
  if found:break
 need(found is not None,'no feedback singleton');triangles.append({'deleted':deleted,'cycle':found})
negative='24c37bda3f05d45bb98f685f570046313ae1a464573932b145b7d7a58f83a3cb'
need(hashlib.sha256((P/'FEEDBACK533.md').read_bytes()).hexdigest()==negative,'negative theorem binding')
out={'status':'PASS','rows':a,'H':H,'C':C,'cycles':cycles,'feedback_number':2,'fixed_ports':ports,'tail':tail,'masks':masks,'checked_arcs':arcs,'singleton_deletion_triangles':triangles,'negative_theorem_sha256':negative,'minimum_C_for_strong_H533_feedback_at_most_two':3}
s=json.dumps(out,sort_keys=True,indent=2)+'\n';(P/'feedback533_c3_certificate.json').write_text(s);print(s,end='')
