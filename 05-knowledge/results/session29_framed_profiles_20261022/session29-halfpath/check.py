#!/usr/bin/env python3
"""Exact half-path DP with independent permutation/cycle/orbit controls."""
import itertools as it, json, math
from fractions import Fraction

def req(x,s):
 if not x: raise RuntimeError(s)
def arc(T,i,j):return bool(T[i]>>j&1)
def anti(T,p):return all(arc(T,i,j)==arc(T,p[j],p[i]) for i in range(len(T)) for j in range(i+1,len(T)))
def auto(T,p):return all(arc(T,i,j)==arc(T,p[i],p[j]) for i in range(len(T)) for j in range(i+1,len(T)))
def dp(T,s):
 n=len(T); pairs=[(i,s[i]) for i in range(n) if i<s[i]];fixed=[i for i in range(n) if s[i]==i];m=len(pairs)
 req(len(fixed)==n%2 and anti(T,s) and all(s[s[i]]==i for i in range(n)),'valid antiinvolution')
 if not m:return {'fixed_hp':1,'fixed_odd_cycles':0,'cycle_by_length':{},'states':0}
 vertices=[v for pair in pairs for v in pair];pows=[3**i for i in range(m)];size=3**m
 D=[[0]*(2*m) for _ in range(size)];V=[[0]*(2*m) for _ in range(size)] if fixed else None
 for j,x in enumerate(vertices):
  code=(1+j%2)*pows[j//2];D[code][j]=1
  if fixed:V[code][j]=int(arc(T,fixed[0],x))
 f=0;cyc={};nonzero=0
 for code in range(1,size):
  digits=[code//p%3 for p in pows];chosen=[2*i+d-1 for i,d in enumerate(digits) if d];q=len(chosen)
  for j in chosen:
   x=vertices[j];count=D[code][j]
   if count:nonzero+=1
   if q==m and arc(T,x,fixed[0] if fixed else s[x]):f+=count
   if fixed and arc(T,x,s[x]):cyc[2*q+1]=cyc.get(2*q+1,0)+V[code][j]
   for k,d in enumerate(digits):
    if d:continue
    for bit in (0,1):
     l=2*k+bit;y=vertices[l]
     if arc(T,x,y):
      nc=code+(bit+1)*pows[k];D[nc][l]+=count
      if fixed:V[nc][l]+=V[code][j]
 return {'fixed_hp':f,'fixed_odd_cycles':sum(cyc.values()),'cycle_by_length':cyc,'states':nonzero}
def normalize(T,p):
 return tuple(sum(1<<j for j in range(len(T)) if arc(T,p[i],p[j])) for i in range(len(T)))
def brute(T,s,with_cycles=True):
 n=len(T);A=[];hp=[];tiles=set()
 for p in it.permutations(range(n)):
  if auto(T,p):A.append(p)
  if all(arc(T,p[i],p[i+1]) for i in range(n-1)):
   hp.append(p);tiles.add(normalize(T,p))
 f=sum(tuple(s[x] for x in reversed(p))==p for p in hp)
 central=sum(all(p[s[i]]==s[p[i]] for i in range(n)) for p in A)
 rev=tuple(reversed(range(n)));blue=sum(anti(U,rev) for U in tiles)
 cyc={}
 if with_cycles:
  for k in range(3,n+1,2):
   for subset in it.combinations(range(n),k):
    for tail in it.permutations(subset[1:]):
     c=(subset[0],)+tail
     if not all(arc(T,c[i],c[(i+1)%k]) for i in range(k)):continue
     b=tuple(s[x] for x in reversed(c));a=b.index(min(b));b=b[a:]+b[:a]
     if b==c:cyc[k]=cyc.get(k,0)+1
 return {'canonical_spine':min(tiles),'H':len(hp),'aut':len(A),'centralizer':central,'fixed_hp':f,'blue_tiles':blue,'tiles':len(tiles),'fixed_odd_cycles':sum(cyc.values()),'cycle_by_length':cyc}
def find_sigma(T):
 for p in it.permutations(range(len(T))):
  if all(p[p[i]]==i for i in range(len(T))) and anti(T,p):return p
 raise RuntimeError('not self-converse')
def check(name,T,s=None,cycles=True):
 n=len(T)
 s=find_sigma(T) if s is None else s;d=dp(T,s);b=brute(T,s,cycles)
 req(d['fixed_hp']==b['fixed_hp'],name+' fixed HP')
 if cycles:req(d['fixed_odd_cycles']==b['fixed_odd_cycles'],name+' fixed cycles')
 req(b['fixed_hp']==b['centralizer']*b['blue_tiles'],name+' quotient')
 req(b['H']==b['aut']*b['tiles'],name+' free action')
 req(b['fixed_hp']%2==1,name+' odd fixed HP')
 if n>=2:req(b['fixed_hp']<=2**(n//2-1)*math.factorial(n//2),name+' half-path upper bound')
 return {'name':name,'outmasks':T,'sigma':s,**b,'dp_nonzero_states':d['states']}
def circulant(n,diffs):return [sum(1<<j for j in range(n) if (j-i)%n in diffs) for i in range(n)]
def rev_tournament(n,bits):
 s=tuple(reversed(range(n)));left={(i,j) for i in range(n) for j in range(i+1,n)};orbs=[]
 while left:
  a=min(left);orb={a,(s[a[1]],s[a[0]])};left-=orb;orbs.append(orb)
 T=[0]*n
 for z,orb in enumerate(orbs):
  for i,j in orb:
   a,b=(i,j) if bits>>z&1 else (j,i);T[a]|=1<<b
 return T,s,len(orbs)
def main():
 cases=[check('C3',circulant(3,{1})),check('C5{1,2}',circulant(5,{1,2})),check('Paley P7',circulant(7,{1,2,4})),check('seed A',[16,9,3,5,14]),check('seed B',[18,8,3,5,14])]
 controls=0;ensemble=[]
 for n in range(1,7):
  _,_,k=rev_tournament(n,0);total_f=0;classes={}
  for bits in range(1<<k):
   T,s,_=rev_tournament(n,bits);row=check(f'n{n}-{bits}',T,s,cycles=n%2==1);controls+=1;total_f+=row['fixed_hp'];classes[row['canonical_spine']]=row
  req(total_f==(1<<k)*math.factorial(n//2),f'ensemble expectation n{n}')
  mass=sum(Fraction(1,row['centralizer']) for row in classes.values());blue=sum(row['blue_tiles'] for row in classes.values())
  req(mass==Fraction(2**(k-n//2),math.factorial(n//2)),f'centralizer mass n{n}')
  req(blue==2**(k-n//2),f'blue mass n{n}')
  ensemble.append({'self_converse_classes':len(classes),'centralizer_reciprocal_mass':str(mass),'global_blue_tiles':blue,'n':n,'orientations':1<<k,'total_fixed_hp':total_f,'mean_fixed_hp':math.factorial(n//2)})
 for n,bits in [(7,374),(7,977),(8,15217),(8,53773)]:
  T,s,_=rev_tournament(n,bits);cases.append(check(f'sample n{n}-{bits}',T,s))
 print(json.dumps({'status':'PASS','exhaustive_fixed_involution_orientations_n1_to_6':controls,'fixed_involution_ensemble':ensemble,'named_and_larger_cases':cases,'verification':'Independent full permutations, normalized spine-tile sets, automorphisms, centralizers, and directed odd-cycle enumeration.'},indent=2))
if __name__=='__main__':main()
