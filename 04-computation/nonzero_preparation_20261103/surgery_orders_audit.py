import sys,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path('04-computation/nonzero_preparation_20261103').resolve()))
import surgery_orders as producer
checks=0
def need(ok,label):
 global checks;checks+=1
 if not ok:raise RuntimeError(label)
def rotate(T,k):
 n=len(T);mask=(1<<n)-1
 return tuple(((T[(i+k)%n]>>k)|(T[(i+k)%n]<<(n-k)))&mask for i in range(n))
def canon(T):return min(rotate(T,k) for k in range(len(T)))
def insert(T,p):
 n=len(T);R=tuple(row|((1<<n) if i==p else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<p),)
 order=list(range(p+1))+[n]+list(range(p+1,n));pos={v:i for i,v in enumerate(order)}
 return canon(tuple(sum(1<<pos[j] for j in range(n+1) if R[v]>>j&1) for v in order))
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
 return H[full],E[full],F
previous={canon((2,4,1))}; records=[];order11=None
for n in range(3,15):
 if n>3:
  fresh=set()
  for T in previous:
   if len(T)==3:parents=[0]
   else:
    high=max(v.bit_count() for v in T);M=[i for i,row in enumerate(T) if row.bit_count()==high]
    need(len(M) in (1,2),'classical one/two maximum ports')
    if len(M)==1:parents=[(M[0]-1)%len(T),M[0]]
    else:
     a,b=M
     if not T[a]>>b&1:a,b=b,a
     need((a+1)%len(T)==b,'maximum vertices consecutive on the directed cycle')
     parents=[(a-1)%len(T),a,b]
   for p in parents:
    Q=insert(T,p)
    need(Q not in fresh,'classical canonical newest inverse gives distinct children')
    fresh.add(Q)
  previous=fresh
 constructed={canon(producer.unrank_one(n,i)) for i in range(producer.U[n])}
 need(previous==constructed,'independent max-degree insertion equals rooted-unranking inventory')
 if n==10:
  from collections import Counter
  hist=Counter();marked=set()
  for T in previous:
   h,ends,F=table(T);full=(1<<n)-1
   reverse=tuple(full^(1<<i)^row for i,row in enumerate(T))
   hr,starts,_=table(reverse);need(hr==h,'independent converse H')
   low={i for i,row in enumerate(T) if row.bit_count()==1}
   entries=low|{(i+1)%n for i in low}
   need(entries=={i for i,v in enumerate(starts) if v==1},'all actual unique-start roots')
   for j in entries:hist[h+F[j]-1]+=1;marked.add(rotate(T,j))
  need(len(marked)==sum(hist.values())==987,'entire independent987 marked universe')
  raw=(json.dumps(sorted(hist.items()),separators=(',',':'))+'\n').encode()
  supplied=json.loads(Path('05-knowledge/results/nonzero_preparation_20261103/surgery/order11.json').read_text())
  need(hashlib.sha256(raw).hexdigest()==supplied['histogram_sha256'],'entire independent order11 histogram')
  need(hist[1429]==0,'actual empty1429 fibre')
  order11={'classes':987,'histogram_sha256':hashlib.sha256(raw).hexdigest(),'minimum':min(hist),'maximum':max(hist),'target_hits':hist[1429]}
 paircount=0;pairrows=[]
 for T in sorted(previous):
  low={i for i,row in enumerate(T) if row.bit_count()==1};entries=sorted(low|{(i+1)%n for i in low})
  need(len(entries)==(3 if n==3 else len(low)+1),'entry inventory')
  if n==3:J=[(0,1)]
  else:J=[(entries[i],entries[j]) for i in range(len(entries)) for j in range(i+1,len(entries))]
  paircount+=len(J);pairrows.extend((T,j) for j in J)
 need(paircount==producer.P[n],'independent complete pair count')
 raw=json.dumps(pairrows,separators=(',',':')).encode()
 records.append({'order':n,'ordinary_classes':len(previous),'pair_classes':paircount,'pair_sha256':hashlib.sha256(raw).hexdigest()})
 print('completed',n,len(previous),paircount,file=sys.stderr,flush=True)
print(json.dumps({'checks':checks,'construction':'independent maximum-outdegree insertion, compared with producer root-deletion unranking','orders':records,'independent_order11':order11},indent=2,sort_keys=True))
