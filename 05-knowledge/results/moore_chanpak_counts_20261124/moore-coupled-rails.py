from pathlib import Path
from runpy import run_path
core=run_path(str(Path(__file__).with_name("moore-compact-rail-final.py")))
sc,heights,flip,comb=(core[k] for k in ("sc","heights","flip","comb"))
from itertools import combinations
from collections import Counter
from functools import cache

@cache
def sig2(r,p,q):
 n=len(r);ext=[i for i in range(n) if i not in (p,q)];full=sum(1<<i for i in ext)
 @cache
 def walk(s,a,b,last):
  v=Counter()
  if s==full and a and b:v[a,b]=1
  for i in range(n):
   if last!=-1 and not(r[last]>>i&1):continue
   u=(s,a+1,b,i) if i==p and a<=len(ext) else ((s,a,b+1,i) if i==q and b<=len(ext) else ((s|1<<i,a,b,i) if i in ext and not(s>>i&1) else None))
   if u is not None:v.update(walk(*u))
  return v
 return walk(0,0,0,-1)

def expo(B):
 D=Counter()
 for (a,b),c in B.items():
  for i in range(1,a+1):
   for j in range(1,b+1):
    D[i*j]+=(-1)**(a+b-i-j)*comb(a,i)*comb(b,j)*c
 return +D if not any(v<0 for v in D.values()) else {k:v for k,v in D.items() if v}

def scan():
 R=(106,100,25,82,35,76,20);Z=(114,100,89,19,34,76,24)
 for src in (R,Z):
  h=heights(src)[0];E=list(combinations(range(7),2));tot=[]
  for prep in [None,*E]:
   r=src if prep is None else flip(src,*prep)
   hr=heights(r)[0]
   if not sc(r):continue
   for final in E:
    t=flip(r,*final);ds=heights(t)[0]-hr
    if abs(ds)!=2 or not sc(t):continue
    for J in E:
     if set(J)&set(final):continue
     a=sig2(r,*J);b=sig2(t,*J)
     d={x:b.get(x,0)-a.get(x,0) for x in set(a)|set(b)}
     if expo(d)!={1:ds}:continue
     assert prep is not None
     lo,up=(a,b) if ds==2 else (b,a)
     loexp, upexp=expo(lo),expo(up)
     if set(loexp)=={1}:continue
     raw={k:v for k,v in d.items() if v}
     print(h,prep,J,final,ds,loexp,'diff=',raw)
     tot.append((tuple(sorted(loexp.items())),len(raw)>1))
  expected=8 if h==147 else 4
  assert len(tot)==expected
  print('total',h,len(tot),Counter(tot))
scan()

M=15015
f=lambda k:(59*2**k-9,32*3**k+30*2**k-23,68*3**k-18*2**k-35,16*9**k-20*6**k+33*4**k-29*2**k+11,66*9**k-180*6**k+187*4**k+64*3**k-196*2**k+71)
S=[{f(k)[j]%M for k in range(1,61)} for j in range(5)]
L=set.union(*S); U={(x+2)%M for x in L}
assert list(map(len,S))==[60,48,44,48,48] and (len(L),len(L&U),len(L|U))==(243,2,484)
print('five rails: 243 lower residues, 484 endpoint residues')
