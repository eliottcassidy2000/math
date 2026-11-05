#!/usr/bin/env python3
"""Standalone exact controls for endpoint-one recursive classification."""
import hashlib
import itertools
import json
from pathlib import Path
from enumerate_states import census,worddata
from crosscheck_join import data,tournaments

HERE=Path(__file__).resolve().parent

def require(test,detail):
 if not test:raise RuntimeError(detail)

def literal(rows):
 n=len(rows);h=0;e=[0]*n;s=[0]*n;l=[0]*n;r=[0]*n
 for p in itertools.permutations(range(n)):
  ok=[bool(rows[p[i]]>>p[i+1]&1) for i in range(n-1)]
  if all(ok):h+=1;s[p[0]]+=1;e[p[-1]]+=1
  for k in range(1,n+1):
   if all(ok[:k-1]) and all(ok[k:]):l[p[k-1]]+=1
  for k in range(n):
   if all(ok[:max(0,k-1)]) and all(ok[k:]):r[p[k]]+=1
 return h,e,s,l,r

def induced(rows,vs):
 return [sum(1<<j for j,b in enumerate(vs) if rows[a]>>b&1) for a in vs]

def reachability(rows):
 n=len(rows);reach=[rows[i]|1<<i for i in range(n)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return reach

def strong(rows):
 return all(x==(1<<len(rows))-1 for x in reachability(rows))

def build(rows,u,v):
 # v=0 denotes empty word; otherwise binary expansion begins in 1.
 w=[] if v==0 else list(map(int,bin(v)[2:]));n=len(rows);m=len(w);r=n+m
 out=rows[:]+[0]*(m+1)
 for j,b in enumerate(w):
  z=n+j;out[z]|=(1<<n)-1
  for k in range(j+1,m):out[z]|=1<<(n+k)
  if b:out[r]|=1<<z
  else:out[z]|=1<<r
 for a in range(n):
  if a==u:out[a]|=1<<r
  else:out[r]|=1<<a
 return out,r

def transfer(st,v):
 h,s,l=st
 if v==0:return h+l-s,h-s,h+l
 g,a,c,q,e=worddata(v)
 return g*h-e*s+(a+1)*l,h,(1+c)*h+q*l

def binary(st,v):
 h,s,l=st;u=h;x=y=c=ph=ps=0;qh=h;ql=l
 for b in map(int,bin(v)[2:]):
  u,x,y,c,ph,ps,qh,ql=(u+(1-b)*(x+ph),2*x+b*h,2*y+b*l,c+(1-b)*qh,b*h,b*s,2*qh,2*ql)
 return u-s+ps+y+l,h,h+c+ql

def decode(rows,v):
 rem=[i for i in range(len(rows)) if i!=v];base=induced(rows,rem);reach=reachability(base)
 # Every vertex reaches every vertex in the sink SCC.
 ci=[j for j in range(len(rem)) if all(reach[i]>>j&1 for i in range(len(rem)))]
 C=[rem[i] for i in ci];Z=[i for i in rem if i not in C]
 Z0=tuple(Z)
 Z.sort(key=lambda z: -sum(bool(rows[z]>>q&1) for q in Z0))
 core=induced(rows,C);inn=[a for a in C if rows[a]>>v&1]
 require(len(inn)==1,('unique_core_inneighbor',rows,v,C,inn))
 u=C.index(inn[0]);w=[int(bool(rows[v]>>z&1)) for z in Z]
 require(not w or w[0]==1,('first_bit',rows,v,w))
 value=0
 for b in w:value=2*value+b
 require(w or len(C)>1,('base_case',rows,v))
 order=C+Z+[v];rebuilt,r=build(core,u,value)
 require(rebuilt==induced(rows,order),('reconstruction',rows,v,order,rebuilt))
 dc=data(core);st=(dc[0][-1],dc[2][u],dc[3][u])
 require(dc[1][u]==1 and strong(core),('core_endpoint',rows,v))
 require(transfer(st,value)==(data(rows)[0][-1],data(rows)[2][v],data(rows)[3][v]),('transfer',rows,v))
 return st,value

def witness(F,st):
 if F[st] is None:return [0],0
 prev,v=F[st];rows,r=witness(F,prev)
 return build(rows,r,v)

def main():
 graphchecks=rootchecks=decomp=0;smallcores=[]
 for n in range(1,6):
  for rows in tournaments(n):
   D=data(rows);L=literal(rows)
   require(L==(D[0][-1],D[1],D[2],D[3],D[4]),('literal_dp',rows,L,D))
   graphchecks+=1;rootchecks+=n
   if strong(rows):
    for r in range(n):
     if D[1][r]!=1:continue
     if n>1:decode(rows,r);decomp+=1
     if n<=4:smallcores.append((rows,r,(D[0][-1],D[2][r],D[3][r])))
 macrochecks=0
 for rows,r,st in smallcores:
  for v in range(16):
   if v==0 and len(rows)==1:continue
   child,root=build(rows,r,v);D=data(child)
   require(strong(child),('macro_strong',rows,r,v))
   require(D[1][root]==1,('macro_endpoint',rows,r,v))
   actual=(D[0][-1],D[2][root],D[3][root]);require(actual==transfer(st,v),('macro_state',rows,r,v,actual))
   if v:require(binary(st,v)==actual,('binary_state',st,v))
   macrochecks+=1
 F13=census(13);F=census(1000)
 expected13=[(1,1,1),(3,1,3),(5,1,5),(5,1,7),(5,2,6),(9,1,9),(9,1,11),(9,1,15),(9,3,9),(9,3,11),(9,4,10),(11,1,13),(11,4,12)]
 require(sorted(F13)==expected13,'states13')
 values=sorted({s[0] for s in F});missing=sorted(set(range(1,1001,2))-set(values))
 require(len(F)==9293 and len(values)==475,'census1000_cardinality')
 expectedmissing=[7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,277,299,421,449,503,533,569,583,637,821]
 require(missing==expectedmissing,'census1000_missing')
 selected=set(F13)
 for h in [17,19,23,35,49,63,127,343,511,959,997]:selected.add(next(s for s in sorted(F) if s[0]==h))
 witnesses=[]
 for st in sorted(selected):
  rows,r=witness(F,st);D=data(rows)
  require(strong(rows) and D[1][r]==1,'witness_shape')
  require((D[0][-1],D[2][r],D[3][r])==st,'witness_state')
  if st in F13:require(literal(rows)==(D[0][-1],D[1],D[2],D[3],D[4]),'literal_witness13')
  witnesses.append({'state':st,'rows':rows,'root':r,'order':len(rows)})
  for v in range(1,256):require(binary(st,v)==transfer(st,v),'binary_extended')
 hostile=[2,12,9,17,7];DH=literal(hostile)
 require(DH[0]==13 and DH[1]==[4,3,2,2,2],'H13_hostile')
 M=[[1,-1,1],[1,-1,0],[1,0,1]]
 def mul(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
 M2=mul(M,M);M3=mul(M2,M)
 require(all(M3[i][j]-M2[i][j]-M[i][j]-(i==j)==0 for i in range(3) for j in range(3)),'tribonacci_matrix')
 seq=[];st=(3,1,3)
 for j in range(12):seq.append(st[0]);st=transfer(st,0)
 require(all(seq[j+3]==sum(seq[j:j+3]) for j in range(9)),'tribonacci_values')
 Fbig=census(10000)
 valuesbig=sorted({t[0] for t in Fbig});missingbig=sorted(set(range(1,10001,2))-set(valuesbig))
 require(len(Fbig)==284844 and len(valuesbig)==4971,'census10000_cardinality')
 require(missingbig==expectedmissing+[1067,1339,1861,2647],'census10000_missing')
 high={'bound':10000,'states':len(Fbig),'values':len(valuesbig),'missing':missingbig,'state_sha256':hashlib.sha256(json.dumps(sorted(Fbig),separators=(',',':')).encode()).hexdigest()}
 statebytes=json.dumps(sorted(F),separators=(',',':')).encode()
 receipt={'census10000':high,'status':'PASS','labelled_graphs_orders_1_to_5':graphchecks,'root_marks':rootchecks,'strong_endpoint_one_decompositions':decomp,'graph_macro_checks':macrochecks,'states_through13':expected13,'census_bound':1000,'state_count':len(F),'value_count':len(values),'missing':missing,'state_sha256':hashlib.sha256(statebytes).hexdigest(),'tribonacci':seq,'hostile13':{'rows':hostile,'terminal_counts':DH[1]},'witnesses':witnesses}
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
 print(json.dumps({k:v for k,v in receipt.items() if k not in ('witnesses','states_through13')},sort_keys=True))

if __name__=='__main__':main()
