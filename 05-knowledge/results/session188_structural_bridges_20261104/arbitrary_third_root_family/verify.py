"""Arbitrary root-containing original-module XOR controls."""
from itertools import combinations
from functools import lru_cache
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

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
  if c>>v&1:b[v]|=1<<(n+1)
  else:b[n+1]|=1<<v
 return tuple(b)


def case(m,C,label):
 T=blow((2,4,1),(rows(m,0),rows(m,0),C));n=len(T);full=(1<<n)-1
 blocks=[(1<<m)-1,((1<<m)-1)<<m,((1<<len(C))-1)<<(2*m)]
 mods=[M for M in range(1,full) if M.bit_count()>1 and module(T,M)]
 need(all(any(M&B==M for B in blocks) for M in mods),'proper module crosses children')
 root_pairs=[(blocks[0],blocks[1])] if m==2 else [(blocks[1],blocks[2]),(blocks[2],blocks[0])]
 battery=[(0,0)]+root_pairs;base=[hp(extend(T,a,c)) for a,c in battery]
 need(base[0]==hp(T),'scalar physical query')
 equal=strict=0;records=[]
 for bits in range(1<<len(mods)):
  U=flip(T,full);selected=[full]
  for j,M in enumerate(mods):
   if bits>>j&1:U=flip(U,M);selected.append(M)
  values=[hp(extend(U,a,c)) for a,c in battery]
  need(values[0]>=base[0],'scalar decrease under original XOR')
  need(values!=base,'fixed root family battery collision')
  if values[0]>base[0]:strict+=1
  else:
   equal+=1;gap=sum(x-y for x,y in zip(base[1:],values[1:]))
   need(gap!=0,'root contrast cancelled after scalar equality')
   if m==2 and len(C)>=2:need(gap%8==4,'root mod8 residue')
   if m==3:need((gap<0)==(len(C)>=2),'TT3 root-sum sign')
  records.append({'family_bits':bits,'modules':selected,'values':values})
 return {'label':label,'m':m,'C':C,'T':T,'proper_modules':mods,'battery':battery,'baseline_values':base,
         'families':len(records),'scalar_equal':equal,'scalar_strict':strict,'records':records}

if __name__=='__main__':
 C5=tuple(sum(1<<((i+j)%5) for j in(1,2)) for i in range(5))
 deep=blow((2,4,1),(rows(2,0),(0,),(0,)))
 cases=[case(2,rows(4,0),'crossing third child'),case(3,(2,4,1),'crossing first two children'),
        case(2,deep,'deep cyclic third child'),case(2,C5,'prime five third child'),
        case(2,(0,),'singleton third child'),case(3,(0,),'TT3 singleton third child')]
 print(json.dumps({'status':'PASS','cases':cases,'scope':'Original-versus-root-containing-family separation, not all-pairs orbit injectivity.'},indent=2,sort_keys=True))
