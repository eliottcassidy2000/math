import json,hashlib
from pathlib import Path
from functools import cache
checks=0
def need(ok,label):
 global checks;checks+=1
 if not ok:raise RuntimeError(label)
def reverse(T):
 full=(1<<len(T))-1
 return tuple(full^(1<<i)^row for i,row in enumerate(T))
def induced(T,ids):return tuple(sum(1<<j for j,w in enumerate(ids) if T[v]>>w&1) for v in ids)
def strong(T):
 for v in range(len(T)):
  seen={v};todo=[v]
  while todo:
   z=todo.pop()
   for j in range(len(T)):
    if T[z]>>j&1 and j not in seen:seen.add(j);todo.append(j)
  if len(seen)!=len(T):return False
 return True
@cache
def table(T):
 n=len(T);full=(1<<n)-1;E=[[0]*n for _ in range(full+1)];H=[0]*(full+1);H[0]=1
 pred=[sum(1<<j for j in range(n) if T[j]>>i&1) for i in range(n)]
 for mask in range(1,full+1):
  for v in range(n):
   if not mask>>v&1:continue
   old=mask^(1<<v)
   if not old:E[mask][v]=1
   else:
    q=old&pred[v]
    while q:
     b=q&-q;q-=b;E[mask][v]+=E[old][b.bit_length()-1]
  H[mask]=sum(E[mask])
 F=[sum(E[mask][v]*H[full^mask] for mask in range(full+1)) for v in range(n)]
 return H[full],E[full],F,sum(H[m]*H[full^m] for m in range(full+1))
p=Path('05-knowledge/results/nonzero_preparation_20261103/surgery/targets.json')
receipt=json.loads(p.read_text());results=[]
for rec in receipt['records']:
 T=tuple(rec['rows']);r=rec['root'];a=rec['aux'];n=len(T)
 need(all(not(row>>i&1) and not(row>>n) for i,row in enumerate(T)) and all(((T[i]>>j&1)+(T[j]>>i&1))==1 for i in range(n) for j in range(i+1,n)),'literal tournament')
 need(strong(T),'actual strong target');need(T[a]==1<<r,'sole outgoing root')
 need(n==rec['order'],'declared witness order')
 h,_,F,_=table(T);_,_,FS,_=table(reverse(T))
 ids=[i for i in range(n) if i!=a];_,_,FL,_=table(induced(T,ids));L=FL[ids.index(r)]
 ids2=[i for i in ids if i!=r];K=table(induced(T,ids2))[3]
 actual={'H':h,'F':F[r],'S':FS[a],'L':L,'K':K,'G':L+K,'root_defect':F[r]-h,'sigma':FS[a]-h,'tau':L+K-F[r]}
 need(actual==rec['observers'],'all nine independent observers')
 need(actual['sigma']==2 and actual['tau']==0,'actual joint guard')
 core=tuple(rec['core']);J=set(rec['ports']);word=rec['word'];m=len(core);t=len(word);rr=m+t
 need(strong(core) and 7<=m<=10,'literal strong core scope')
 rows=[row|(((1<<t)-1)<<m)|(0 if i in J else 1<<rr) for i,row in enumerate(core)]
 rows += [sum(1<<(m+j) for j in range(i+1,t))|((1<<rr) if bit else 0) for i,bit in enumerate(word)]
 rows += [sum(1<<v for v in J)+sum(1<<(m+i) for i,bit in enumerate(word) if not bit)]
 nn=len(rows);rows=tuple(row|((1<<nn) if i==rr else 0) for i,row in enumerate(rows))+(((1<<nn)-1)^(1<<rr),)
 need(reverse(rows)==T and rr==r and nn==a,'independent funnel/birth/converse construction')
 need(sum(table(reverse(core))[1][j] for j in J)==2,'actual core entry budget')
 results.append({'H':h,'order':n,'root':r,'aux':a,'observers':actual})
expected=sorted(27*1429+9*a+3*b+c for a in (2,4) for b in (2,4) for c in (2,4))
need(sorted(r['H'] for r in results)==expected,'exact eight third-generation scalar values')
print(json.dumps({'checks':checks,'receipt_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'records':results},sort_keys=True,indent=2))
