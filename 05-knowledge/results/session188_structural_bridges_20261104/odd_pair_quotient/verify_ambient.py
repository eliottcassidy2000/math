"""Focused literal checks for the separate arbitrary-ambient corollary."""
from pathlib import Path
from itertools import combinations,product
import json,hashlib,importlib.util
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('odd_pair_base',HERE/'verify.py');base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
need,hp,ext=base.need,base.hp,base.ext

def root(K,r,A):
 outside=[v for v in range(len(K)) if v!=r];R=list(A)+[0]*len(outside);n=len(A)
 for i,v in enumerate(outside):
  for j,w in enumerate(outside):
   if K[v]>>w&1:R[n+i]|=1<<(n+j)
  if K[v]>>r&1:R[n+i]|=(1<<n)-1
  else:
   for u in range(n):R[u]|=1<<(n+i)
 return tuple(R)
def pairbits(R,bits):
 R=list(R)
 for i in range(3):
  if bits>>i&1:R[2*i]^=1<<(2*i+1);R[2*i+1]^=1<<(2*i)
 return tuple(R)
def main():
 Q=(2,4,1);q=3;A,B=base.inflate(Q),base.inflate(base.op(Q));coarse=[]
 for j in range(q):
  for a in [0]+[1<<i for i in range(q)]+[(1<<i)|(1<<j) for i in range(q) if Q[j]>>i&1]:coarse.append((a,1<<j))
 coarse += [(a,0) for a in [0]+[1<<i for i in range(q)]+[(1<<i)|(1<<j) for i,j in combinations(range(q),2)]]
 need(len(coarse)==len(set(coarse))==2*q*q+q+1,'root battery count')
 contexts=[((0,),0),((2,0),0),((2,0),1),((2,4,1),0),((6,4,0),1),((24,17,3,6,12),0)];records=[]
 for K,r in contexts:
  T,U=root(K,r,A),root(K,r,B);outside=sum(1<<u for u in range(6,len(T)) if T[u]&1)
  queries=[(outside|sum(3<<(2*i) for i in range(q) if a>>i&1),outside|sum(3<<(2*i) for i in range(q) if c>>i&1)) for a,c in coarse]
  values=[(hp(ext(T,a,c)),hp(ext(U,a,c))) for a,c in queries];d={m:f-g for m,(f,g) in zip(coarse,values)};L=[]
  for j in range(q):
   av=1<<j;s=0
   for i in range(q):
    if i==j:continue
    au=1<<i;s+=d[au,av]-d[au,0]-d[0,av]+d[0,0]
    if Q[j]>>i&1:s+=d[au|av,av]-d[au,av]-d[av,av]+d[0,av]-d[au|av,0]+d[au,0]+d[av,0]-d[0,0]
   L.append(s)
  factor=hp(root(K,r,(2,0)));need(factor%2==1 and sum(L)%8==4,'ambient trace and odd unit')
  for (a,c),(f,g) in zip(coarse,values):
   aa=sum(3<<(2*i) for i in range(q) if a>>i&1);cc=sum(3<<(2*i) for i in range(q) if c>>i&1)
   need((f-g-factor*(hp(ext(A,aa,cc))-hp(ext(B,aa,cc))))%8==0,'query-independent unit transfer')
  records.append({'outer_quotient':K,'root':r,'source':T,'target':U,'coarse_queries':coarse,'physical_queries':queries,'source_target_values':values,'L':L,'unit':factor})
 # Complete physical two-level cube in the genuinely cyclic outer context.
 rec=records[3];T0,T1=tuple(rec['source']),tuple(rec['target']);queries=[tuple(x) for x in rec['physical_queries']];decoders=[]
 for theta,T in enumerate((T0,T1)):
  for i in range(q):
   mask=3<<(2*i);incoming=sum(1<<u for u in range(len(T)) if not mask>>u&1 and T[u]>>(2*i)&1);U=pairbits(T,1<<i);choice=None
   for a,c in product(range(4),repeat=2):
    query=(incoming|(a<<(2*i)),incoming|(c<<(2*i)));f,g=hp(ext(T,*query)),hp(ext(U,*query))
    if (f-g)%4==2:choice=query;break
   need(choice is not None,'ambient pair isolation');decoders.append({'theta':theta,'child':i,'query_index':len(queries),'values':[f,g]});queries.append(choice)
 states=[];seen=set()
 for theta,T in enumerate((T0,T1)):
  for bits in range(8):
   R=pairbits(T,bits);values=[hp(ext(R,*v)) for v in queries];need(values[:len(coarse)]==[v[theta] for v in rec['source_target_values']],'ambient root independence')
   for d in decoders:
    if d['theta']==theta:need(values[d['query_index']]==d['values'][(bits>>d['child'])&1],'ambient child independence')
   need(tuple(values) not in seen,'ambient full cube signatures');seen.add(tuple(values));states.append({'theta':theta,'bits':bits,'physical_rows':R,'query_values':values})
 need(len(queries)==2*q*q+3*q+1,'fixed query bound')
 names=[HERE/'AMBIENT.md',HERE/'verify_ambient.py',HERE/'PROOF.md',HERE/'verify.py',HERE.parent/'mod8_transfer'/'PROOF.md'];pins={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in names}
 result={'status':'PASS','contexts':records,'ambient_cube':{'states':states,'queries':queries,'decoders':decoders,'adaptive_bound':2*q*q+2*q+1,'fixed_count':len(queries)},'source_sha256':pins}
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
