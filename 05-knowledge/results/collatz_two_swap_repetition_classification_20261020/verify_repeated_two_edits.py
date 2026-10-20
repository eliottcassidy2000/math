#!/usr/bin/env python3
"""Exact full repetition-factor syndrome for two actual adjacent swaps."""
import itertools,json,math
from collections import Counter
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ArithmeticError(msg)
def carry(w):
 b=0
 for i,t in enumerate(w):b=(3 if t else 1)*b+(t<<i)
 return b
def edits(w):
 for i in range(len(w)):
  j=(i+1)%len(w)
  if w[i]!=w[j]:
   v=list(w);v[i],v[j]=v[j],v[i];yield i,tuple(v)
def rot(w,i):return w[i:]+w[:i]
def necklace(w):return min(rot(w,i) for i in range(len(w)))
def run(v,d):
 global sources,paths,bases,rootpaths
 q0=len(v);p0=sum(v);A=2**q0;C=3**p0
 H=sum(A**(d-1-j)*C**j for j in range(d));c=v*d;q=len(c);D=2**q-C**d
 need(H>=7 and math.gcd(H,6)==1,'bad repetition factor')
 need(D==H*(A-C) and carry(c)==H*carry(v),'repetition law')
 sources+=1;mechanical=math.gcd(q0,p0)==1 and v==tuple(((i+1)*p0)//q0-(i*p0)//q0 for i in range(q0))
 for i,u in edits(c):
  for j,w in edits(u):
   paths+=1
   if q0<=7 and d<=4:rootpaths+=1
   gH=math.gcd(H,carry(w));pred=(i==j) or(d==2 and (j-i)%q==q0)
   need((gH==H)==pred,'H classification')
   if i==j:
    need(w==c and gH==H,'undo');cl='undo'
   elif (j-i)%q in [1,q-1]:
    start=i if (j-i)%q==1 else j
    cr=rot(c,start);wr=rot(w,start)
    delta=carry(wr[:3])-carry(cr[:3]);need(abs(delta) in[3,5],'overlap type')
    need(carry(wr)-carry(cr)==delta*3**sum(cr[3:]),'overlap carry')
    need(gH==math.gcd(H,abs(delta))<H,'overlap gcd depth')
    cl='overlap_smooth3' if abs(delta)==3 else 'overlap_five'
   else:
    start,other=(i,j) if(j-i)%q<=q//2 else(j,i)
    ell=(other-start)%q;cr=rot(c,start);wr=rot(w,start)
    m=sum(cr[2:ell+2]);suffix=sum(cr[ell+2:]);eps0=cr[0]-cr[1];eps1=cr[ell]-cr[ell+1]
    U=eps0*3**m+eps1*2**ell
    need(2<=ell<=q//2 and m>=1 and U!=0,'shorter arc')
    need(carry(wr)-carry(cr)==3**suffix*U,'disjoint exact perturbation')
    need(gH==math.gcd(H,U),'disjoint prime-power gcd')
    if d>=3:
     need(ell<(d-1)*q0 and m<=(d-1)*p0,'periodic interval bound')
     need(abs(U)<A**(d-1)+C**(d-1)<=H,'strict repetition bound')
    else:
     need(ell<=q0 and m<=p0 and abs(U)<=H,'double period bound')
     need((abs(U)==H)==(ell==q0 and m==p0 and eps0==eps1),'equality boundary')
    cl='disjoint_mirrored' if pred else 'disjoint_not_mirrored'
    if pred:
     need(d==2 and w[:q0]==w[q0:],'square output')
     need(carry(w)==H*carry(w[:q0]),'square carry')
     need(any(z==w[:q0] for e,z in edits(v)),'edited base')
   classes[cl]+=1
   if gH==H:divisible[cl]+=1
   if i==q-1 or j==q-1:crossings[cl]+=1
   if mechanical:
    integral=carry(w)%D==0
    expected=((p0,q0) in[(1,2),(2,3)]) and pred
    need(integral==expected,'mechanical integer classification')
    if integral:need(necklace(w)==necklace(c),'integral result not balanced return')
counts=Counter();classes=Counter();divisible=Counter();crossings=Counter();sources=paths=bases=rootpaths=0
for q0 in range(2,9):
 for v in itertools.product((0,1),repeat=q0):
  if sum(v) in[0,q0]:continue
  bases+=1
  for d in range(2,6):run(v,d)
# Explicit controls of proper prime-power retention and the d=2 exception.
hostiles=[]
for name,v,d,edges in [('overlap retains one5','0011',2,[3,2]),('arbitrary base restores minus17','11101111000',2,[3,14])]:
 base=tuple(map(int,v));c=base*d;w=c
 for edge in edges:w=next(z for e,z in edits(w) if e==edge)
 q0=len(base);p0=sum(base);H=sum(2**(q0*(d-1-j))*3**(p0*j) for j in range(d));D=2**len(c)-3**sum(c);bw=carry(w)
 dist=0 if necklace(w)==necklace(c) else 1 if necklace(w) in{necklace(z) for e,z in edits(c)} else 2
 if name.startswith('overlap'):need(H==25 and math.gcd(H,bw)==5 and bw==680 and D==175,'primepower hostile')
 else:need(H==4235 and bw==-17*D and dist==2,'true repeated integer hostile')
 hostiles.append(dict(name=name,base=v,repetitions=d,edges=edges,source=''.join(map(str,c)),output=''.join(map(str,w)),H=H,B=bw,D=D,H_gcd=math.gcd(H,bw),full_content=math.gcd(abs(D),bw),necklace_distance=dist,fixed_point=bw//D if bw%D==0 else None))
need(rootpaths==81592,'independent root subset count')
print(json.dumps(dict(status='PASS_EXACT',checks=checks,universe='all nonconstant binary bases length2..8, d2..5, all ordered actual cyclic adjacent swaps',bases=bases,sources=sources,paths=paths,classes=dict(classes),H_divisible=dict(divisible),crossing_cut_classes=dict(crossings),independent_root_subset_paths=rootpaths,hostiles=hostiles),indent=2))
