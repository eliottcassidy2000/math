#!/usr/bin/env python3
"""Two-prime compactification of the inherited physical reset. Standard library.
core.py and transfer.py are verbatim Session286 inherited literal definitions.
"""
from fractions import Fraction as Q
from itertools import product
from random import Random
import json
import core as c
import transfer as tr
q=Q(104,471); u=Q(13,471)
A=Q(-34104960,367); B=Q(3945600,367); C=Q(-140160,367)
D=Q(-27648000,471**2); E=Q(-829440,471)
CUTS={'A':9,'B':3,'E':18,'C':24}
TS={s:c.attach(tr.D,m) for s,m in CUTS.items()}
PS={s:c.packet(T) for s,T in TS.items()}
TABS={s:tr.table(T) for s,T in TS.items()}
FINAL=(tr.table(TS['A'],1),tr.table(TS['A'],2))
def ck(x,s):
 if not x: raise RuntimeError(s)
def vp(x,p):
 x=Q(x)
 if not x:return 1000000
 a,b=abs(x.numerator),x.denominator; v=0
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v
def mod(x,m):
 x=Q(x);return x.numerator*pow(x.denominator,-1,m)%m
def word(s):
 z=[0]+list(s)+[0,0]
 return ''.join({(0,0):'A',(0,1):'B',(1,0):'C',(1,1):'E'}[a,b] for a,b in zip(z,z[1:]))
def physical(s):
 W=word(s); st={(0,0,0):1}
 for a in W[:-1]:st=tr.step(st,TABS[a])
 V=sum(tr.finish(st,t) for t in FINAL)
 H=c.word([PS[a] for a in W])['H']
 return H,H+2*V,V

def finite(s):
 n=len(s)
 g=sum((A*q**i+B*q**(n-1-i)+C*q**(n-1))*x for i,x in enumerate(s))
 for i in range(n):
  for j in range(i+1,n):
   g+=(E*q**(n-2) if j==i+1 else D*q**(n-1-j+i))*s[i]*s[j]
 j=sum(s)-9*sum(s[i]*s[i+1] for i in range(n-1))
 return g,g+13*u**(n-1)*j

def surface(x,y):
 X=sum(q**i*a for i,a in enumerate(x));Y=sum(q**i*a for i,a in enumerate(y))
 return A*X+B*Y+D*X*Y

def decode_surface(residue,r,p):
 # Recover only r address pairs; input is z mod13^r or z/128 mod8^r.
 base=p if p==13 else 8; scale=1 if p==13 else 128
 x=[];y=[]
 for t in range(r):
  m=base**(t+1)
  cand=[]
  for a,b in product((0,1),repeat=2):
   if mod(surface(x+[a],y+[b])/scale,m)==residue%m:cand.append((a,b))
  if len(cand)!=1:return None
  a,b=cand[0];x.append(a);y.append(b)
 return x,y

STATES=list(product((0,1),repeat=2))
D13=[{mod((A+D*y)*a+(B+D*x)*b,13) for a,b in STATES} for x,y in STATES]
def root_triples():
 out={}
 for ids in product(range(4),repeat=3):
  if len({sum(v)%13 for v in product(*(D13[i] for i in ids))})==13:
   root=sum(mod(surface([STATES[i][0]],[STATES[i][1]]),13) for i in ids)%13
   out.setdefault(root,ids)
 return out
ROOTS=root_triples()
def sum_decode(target,r,p):
 base=13 if p==13 else 8;scale=1 if p==13 else 128;copies=3 if p==13 else 2
 xs=[[] for _ in range(copies)];ys=[[] for _ in range(copies)]
 for t in range(r):
  options=[tuple(STATES[i] for i in ROOTS[target%13])] if t==0 and p==13 else product(STATES,repeat=copies)
  found=None
  for bits in options:
   z=sum(surface(xs[i]+[bits[i][0]],ys[i]+[bits[i][1]]) for i in range(copies))/scale
   if mod(z,base**(t+1))==target%(base**(t+1)):
    found=bits;break
  ck(found is not None,'additive lift exists')
  for i,(a,b) in enumerate(found):xs[i].append(a);ys[i].append(b)
 return list(zip(xs,ys))

def main():
 ck(vp(q,13)==1 and vp(q,2)==3,'two contractions')
 ck([mod(z,13) for z in (A,B,D)]==[11,3,4],'13 leading')
 ck([mod(z/128,8) for z in (A,B,D)]==[5,7,0],'8 leading')
 ck(vp(E-D,13)==1 and vp((E-D)/128,2)==4,'adjacent correction')
 status13={}
 for x,y in product((0,1),repeat=2):
  vals=[mod((A+D*y)*a+(B+D*x)*b,13) for a,b in product((0,1),repeat=2)]
  ck(len(set(vals))==4,'all boundary Jacobian status')
  status13[f'{x}{y}']=vals
 records=[]; total=0
 for n in range(2,9):
  baseline=physical([0]*n); residues13=set();residues8=set()
  r=n//2
  for s in product((0,1),repeat=n):
   H,F,V=physical(s); g,f=finite(s)
   ck(g==Q(V-baseline[2],471**(n-1)),'literal strand response')
   ck(f==Q(F-baseline[1],2*471**(n-1)),'reset response')
   left=s[:n-r];right=s[n-r:][::-1]; z=surface(left,right)
   ck(vp(f-z,13)>=r,'uniform physical13 approximation')
   ck(vp((g-z)/128,2)>=3*r,'uniform physical2 approximation')
   a=mod(f,13**r); b=mod(g/128,8**r)
   ck(decode_surface(a,r,13)==(list(s[:r]),list(s[-r:][::-1])),'physical13 decode')
   ck(decode_surface(b,r,2)==(list(s[:r]),list(s[-r:][::-1])),'physical2 decode')
   residues13.add(a);residues8.add(b);total+=1
  ck(len(residues13)==4**r and len(residues8)==4**r,'exact depth counts')
  records.append(dict(n=n,k=n+2,vertices=5*(n+2)+3,depth=r,residues13=len(residues13),residues8=len(residues8)))
 # Exhaustive address tree, without constructing all enormous tournaments.
 trees=[]
 for r in range(1,7):
  R13=set();R8=set();paired=set()
  for bits in product((0,1),repeat=2*r):
   x=list(bits[:r]);y=list(bits[r:]);z=surface(x,y)
   a=mod(z,13**r);b=mod(z/128,8**r)
   ck(decode_surface(a,r,13)==(x,y),'surface13 inverse')
   ck(decode_surface(b,r,2)==(x,y),'surface2 inverse')
   R13.add(a);R8.add(b);paired.add((a,b))
  ck(len(R13)==len(R8)==len(paired)==4**r,'correlated paired support')
  trees.append(dict(depth=r,count=len(paired),ambient13=13**r,ambient8=8**r,ambient_product=104**r))
 rng=Random(28713);metric_checks=0
 for r in range(1,15):
  for _ in range(80):
   x=[rng.randrange(2) for _ in range(r)];y=[rng.randrange(2) for _ in range(r)]
   xx=[rng.randrange(2) for _ in range(r)];yy=[rng.randrange(2) for _ in range(r)]
   ts=[i for i in range(r) if (x[i],y[i])!=(xx[i],yy[i])]
   if not ts:continue
   t=min(ts);delta=surface(x,y)-surface(xx,yy)
   ck(vp(delta,13)==t,'exact13 metric')
   digit=5*(x[t]-xx[t])+7*(y[t]-yy[t])
   ck(vp(delta/128,2)==3*t+vp(digit,2),'exact2 leading digit metric')
   metric_checks+=1
 ck(len(ROOTS)==13,'good derivative root triples cover all13')
 ck(len({(a+b)%8 for a in (0,4,5,7) for b in (0,4,5,7)})==8,'twofold base8 sum')
 ck(len({(a+b)%13 for a in (0,3,11,5) for b in (0,3,11,5)})==9,'twofold13 not full')
 sumchecks=0; sumphysical=[]
 for p in (2,13):
  base=8 if p==2 else 13;scale=128 if p==2 else 1
  for r in range(1,8):
   targets=range(base**r) if r<=2 else [rng.randrange(base**r) for _ in range(12)]
   for z in targets:
    pairs=sum_decode(z,r,p)
    ck(mod(sum(surface(x,y) for x,y in pairs)/scale,base**r)==z,'additive target inverse')
    sumchecks+=1
   if r<=4:
    z=base**r-1;pairs=sum_decode(z,r,p);v=0
    baseline=physical([0]*(2*r))
    for x,y in pairs:
     H,F,V=physical(x+y[::-1])
     v+=Q(F-baseline[1],2*471**(2*r-1)) if p==13 else Q(V-baseline[2],128*471**(2*r-1))
    ck(mod(v,base**r)==z,'physical additive residue lift')
    sumphysical.append(dict(prime=p,depth=r,target=z,copies=len(pairs)))
 print(json.dumps(dict(status='PASS',inherited_physical_response_cases=total,leading13=[11,3,4],leading8=[5,7,0],boundary_status13=status13,physical_layers=records,cylinder_counts=trees,metric_checks=metric_checks,additive_basis_orders={'2':2,'13':3},good_root_triples=ROOTS,additive_inverse_checks=sumchecks,physical_sum_controls=sumphysical,limits=['normalization and central subtraction are retained','paired prime residues share one address','no Chan-Pak coverage']),indent=2,sort_keys=True))
if __name__=='__main__':main()
