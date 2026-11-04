"""Literal controls for even cyclic-core ambient synthesis."""
from itertools import combinations
from functools import lru_cache
from pathlib import Path
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def rows(n,code):
 a=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(j,i) if code>>bit&1 else(i,j);a[u]|=1<<v
 return tuple(a)

def induced(a,mask):
 names=[i for i in range(len(a)) if mask>>i&1]
 return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)

def op(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))

def module(a,M):return all(r&M in(0,M) for i,r in enumerate(a) if not M>>i&1)

def flip(a,M):return tuple(r^(M^(1<<i)) if M>>i&1 else r for i,r in enumerate(a))

@lru_cache(None)
def hp(a):
 n=len(a);full=(1<<n)-1
 @lru_cache(None)
 def finish(mask,last):
  if mask==full:return 1
  out=0;rem=a[last]&~mask
  while rem:
   b=rem&-rem;rem-=b;out+=finish(mask|b,b.bit_length()-1)
  return out
 return sum(finish(1<<i,i) for i in range(n)) if n else 1

def blow(Q,blocks):
 offsets=[0]
 for b in blocks:offsets.append(offsets[-1]+len(b))
 a=[0]*offsets[-1]
 for i,b in enumerate(blocks):
  for u,r in enumerate(b):
   a[offsets[i]+u]=r<<offsets[i]
   for j,c in enumerate(blocks):
    if Q[i]>>j&1:a[offsets[i]+u]|=((1<<len(c))-1)<<offsets[j]
 return tuple(a)

def quotient(a,blocks):
 reps=[(M&-M).bit_length()-1 for M in blocks]
 return tuple(sum(1<<j for j,v in enumerate(reps) if a[u]>>v&1) for u in reps)

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
  if c>>v&1:b[v]|=1<<(n+1)
  else:b[n+1]|=1<<v
 return tuple(b)

def tree(a):
 mods=[M for M in range(1,1<<len(a)) if module(a,M)]
 strong=[M for M in mods if all(M&N in(0,M,N) for N in mods)]
 nodes={}
 for M in strong:
  if M.bit_count()<2:continue
  inside=[N for N in strong if N!=M and N&M==N]
  children=[N for N in inside if not any(N!=K and N&K==N for K in inside)]
  nodes[M]=tuple(children)
 return mods,nodes

def highest(T,U,nodes):
 changed=[M for M,ch in nodes.items() if quotient(T,ch)!=quotient(U,ch)]
 return [M for M in changed if not any(M!=N and M&N==M for N in changed)]


def put(Q,r,A):return blow(Q,tuple(A if i==r else(0,) for i in range(len(Q))))

def context_case(C,Q,r):
 M=blow((2,4,1),(rows(2,0),rows(2,0),C));m=len(M);full=(1<<m)-1
 mods=[S for S in range(1,full) if S.bit_count()>1 and module(M,S)]
 T=put(Q,r,M);Mmask=full<<r;rep=r
 need(module(T,Mmask),'core not ambient module')
 for S in mods:need(module(T,S<<r),'local original module not global')
 E=sum(1<<v for v in range(len(T)) if not Mmask>>v&1 and T[v]>>rep&1)
 a,c=E|(3<<r),E|(12<<r);baseline=[hp(T),hp(extend(T,a,c))]
 need(hp(extend(T,0,0))==baseline[0],'global scalar query')
 values=[];equal=strict=0
 for bits in range(1<<len(mods)):
  local=flip(M,full)
  for j,S in enumerate(mods):
   if bits>>j&1:local=flip(local,S)
  U=put(Q,r,local);hu=hp(U);fu=hp(extend(U,a,c))
  need(hu>=baseline[0],'global scalar decrease')
  if hu==baseline[0]:
   equal+=1;need((baseline[1]-fu)%8==4,'even core lifted residue')
  else:strict+=1
  need([hu,fu]!=baseline,'two-query family collision');values.append([hu,fu])
 return {'C':C,'Q':Q,'root':r,'local_original_modules':mods,'query':[a,c],'baseline':baseline,
         'families':len(values),'scalar_equal':equal,'scalar_strict':strict,'all_values_by_family_bits':values}

def highest_control():
 M=blow((2,4,1),(rows(2,0),rows(2,0),rows(2,0)));T=blow((2,4,1),(M,(2,4,1),(0,)))
 core=63;other=7<<6;V=flip(T,core)
 for child in(3,12,48):V=flip(V,child)
 U=flip(V,other)
 for selected in(core,other,3,12,48):need(module(T,selected),'not original selected module')
 _,nodes=tree(T);tops=highest(T,U,nodes)
 need(set(tops)=={core,other},'not the two claimed highest nodes')
 need(len(quotient(T,nodes[core]))==3,'core quotient is not C3')
 E=sum(1<<v for v in range(len(T)) if not core>>v&1 and T[v]&1)
 a,c=E|3,E|12;ft,fv,fu=[hp(extend(X,a,c)) for X in(T,V,U)]
 need(hp(T)==hp(U),'highest-node control not scalar equal')
 need(fu==fv,'other highest-node change did not preserve physical query')
 need((ft-fu)%8==4,'highest prime even-core residue lost')
 return {'T':T,'V':V,'U':U,'highest_nodes':tops,'query':[a,c],'values':[ft,fv,fu],'scalar':hp(T)}

def main():
 contexts=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):
   for r in range(n):contexts.append((rows(n,code),r))
 for code in(0,17,45):
  for r in range(4):contexts.append((rows(4,code),r))
 cases=[context_case(rows(2,0),Q,r) for Q,r in contexts]
 cases.append(context_case(rows(4,0),(2,4,1),0))
 Q5=tuple(sum(1<<((i+j)%5) for j in(1,2)) for i in range(5))
 cases.append(context_case((8,1,3,6),Q5,2))
 return {'status':'PASS','contexts':len(cases),'families':sum(x['families'] for x in cases),
         'scalar_equal':sum(x['scalar_equal'] for x in cases),'scalar_strict':sum(x['scalar_strict'] for x in cases),
         'cases':cases,'multiple_highest_prime_control':highest_control(),
         'scope':'Even C3[TT2,TT2,C] ambient root-family separation and highest-node exclusion; not arbitrary prime quotient reconstruction.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
