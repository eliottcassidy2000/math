import json,hashlib,sys
from pathlib import Path
from array import array
D=Path(__file__).resolve().parent
C=D.parent/'closure'/('candidate'+sys.argv[1]+'.json')
c=json.loads(C.read_text()); s=json.loads((D/'source.json').read_text())
def req(x,m):
 if not x: raise RuntimeError(m)
def tables(a):
 n=len(a); size=1<<n; E=array('Q',[0])*(size*n); H=array('Q',[0])*size; H[0]=1
 for m in range(1,size):
  bs=m;total=0
  while bs:
   q=bs&-bs;bs-=q;v=q.bit_length()-1;sub=m^q
   if not sub:val=1
   else:
    pred=sub & (((1<<n)-1)^a[v]^(1<<v));val=0
    while pred:
     p=pred&-pred;pred-=p;val+=E[sub*n+p.bit_length()-1]
   E[m*n+v]=val;total+=val
  H[m]=total
 return E,H
names=c['names'];a=c['rows'];n=len(a);V=(1<<n)-1
req(names[:len(s['names'])]==s['names'] and len(set(names))==len(names),'old names preserved')
for i in range(n):
 req(not a[i]>>i&1,'loop')
 for j in range(i):req(((a[i]>>j)^(a[j]>>i))&1,'tournament')
changes=[[names[i],names[j]] for i in range(len(s['names'])) for j in range(i+1,len(s['names'])) if (a[i]^s['rows'][i])>>j&1]
req({frozenset(x) for x in changes}=={frozenset(x) for x in c['all_changed_pairs']},'flip list')
for i in range(n):
 seen=1<<i
 while True:
  new=seen
  for j in range(n):
   if seen>>j&1:new|=a[j]
  if new==seen:break
  seen=new
 req(seen==V,'strong')
r=names.index(c['root']);P=[names.index(x) for x in c['prefix']];B=[names.index(x) for x in c['suffix_names']];W=[names.index(x) for x in c['in_ports']]
req([names[i] for i in range(n) if a[i]>>r&1]==c['in_ports'],'physical inports')
req(all(a[r]>>i&1 for i in P),'root beats prefix')
req(all(a[i]>>j&1 for i in P for j in B),'prefix beats suffix')
req(all(a[P[i]]>>P[j]&1 for i in range(len(P)) for j in range(i+1,len(P))),'prefix transitive')
E,H=tables(a);req(H[V]==c['target_H'],'full H')
root_E2=sum(E[m*n+r]*H[V^m] for m in range(1,V+1))
b=[sum(((a[i]>>j)&1)<<k for k,j in enumerate(B)) for i in B];bn=len(b);bv=(1<<bn)-1;w=[B.index(x) for x in W]

for i in range(bn):
 seen=1<<i
 while True:
  new=seen
  for j in range(bn):
   if seen>>j&1:new|=b[j]
  if new==seen:break
  seen=new
 req(seen==bv,'suffix strong')
BE,BH=tables(b);rev=[bv^x^(1<<i) for i,x in enumerate(b)];BS,_=tables(rev)
F=J=0
for m in range(1,bv+1):
 end=sum(BE[m*bn+i] for i in w);F+=end*BH[bv^m];J+=end*sum(BS[(bv^m)*bn+i] for i in w)
h=BH[bv]
req((h,F,J)==(c['h'],c['F'],c['J']),'joint fields')
req(h-J+2**len(P)*F==H[V],'attachment formula')
req(root_E2==c['root_E2'] and root_E2-H[V]==J,'root defect')
_,SH=tables(s['rows']);req(SH[-1]==c['source_H'],'source H')
report={'status':'PASS','source_H':int(SH[-1]),'target_H':int(H[V]),'h':int(h),'F':F,'J':J,'root_E2':root_E2,'root_defect':root_E2-int(H[V]),'strong':True,'physical_inports':c['in_ports'],'old_names_preserved':True,'changed_pairs':changes,'candidate_sha256':hashlib.sha256(C.read_bytes()).hexdigest(),'native_compilation':'independently replayed in native_audit.py'}
(D/('audit'+sys.argv[1]+'.json')).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
