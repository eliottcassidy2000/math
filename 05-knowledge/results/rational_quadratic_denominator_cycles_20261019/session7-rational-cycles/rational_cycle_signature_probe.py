#!/usr/bin/env python3
"""Exact denominator-first cycle search; no cutoff on numerator C of c=C/b^2.
A rigorous theorem bounds period by 2^omega(b)+2 (Eliahou--Fares).
This script uses a proved superset of the three-edge local signature table.
"""
from itertools import product
from math import gcd
from fractions import Fraction
import json, argparse, time

def factor(n):
 out=[];p=2
 while p*p<=n:
  if n%p==0:
   e=0
   while n%p==0:n//=p;e+=1
   out.append((p,e))
  p+=1
 if n>1:out.append((n,1))
 return out

def signatures(N):
 s={(0,0,0)}
 for k in range(2,N):s.update(((0,1,k),(0,k,1),(k,k-1,0)))
 for k in range(1,N):s.add((1,0,k))
 return sorted(s)

def canonical(a):return min(tuple(a[i:]+a[:i]) for i in range(len(a)))

def search(b,N=None):
 fs=factor(b);assert b%4==0
 N=N or 2**len(fs)+2
 sig=signatures(N);h=b//2
 local=[[(p**((e-1 if p==2 else e)*i),p**((e-1 if p==2 else e)*j),p**((e-1 if p==2 else e)*k)) for i,j,k in sig] for p,e in fs]
 triangles=set();found=set();count=0;primitive_triangles=0
 for parts in product(*local):
  A=B=D=1
  for a,bb,d in parts:A*=a;B*=bb;D*=d
  count+=1
  for u in (A,-A):
   for v in (B,-B):
    t=u*v-h*u+h*v
    if t!=u*D and t!=-u*D:continue
    w=t//u
    a0=u+v-w;a1=u+w-v;a2=v+w-u
    if len({a0,a1,a2})<3:continue
    if any(gcd(a,b)>1 for a in (a0,a1,a2)):continue
    C=b*a1-a0*a0
    if gcd(C,b)>1:continue
    assert a1*a1+C==b*a2
    triangle=(a0,a1,a2,C)
    if triangle in triangles:continue
    triangles.add(triangle);primitive_triangles+=1
    orbit=[];x=a0
    for step in range(N+1):
     if x in orbit:
      if x==a0 and len(orbit)>=3:
       found.add((C,canonical(orbit)))
      break
     if step==N:break
     orbit.append(x)
     numerator=x*x+C
     if numerator%b:break
     x=numerator//b
     if gcd(x,b)>1:break
 return {'b':b,'omega':len(fs),'period_bound':N,'signature_count':len(sig),'unsigned_triangle_count':count,'primitive_triangle_candidates':primitive_triangles,'cycles':[{'C':c,'c':str(Fraction(c,b*b)),'numerators':list(a),'period':len(a)} for c,a in sorted(found)]}

def main():
 p=argparse.ArgumentParser();p.add_argument('--max-b',type=int,default=400);p.add_argument('--b',type=int,nargs='*');p.add_argument('--max-omega',type=int,default=3);p.add_argument('--output');args=p.parse_args()
 bs=args.b or [b for b in range(4,args.max_b+1,4) if 3<=len(factor(b))<=args.max_omega]
 results=[]
 for b in bs:
  if len(factor(b))>args.max_omega:continue
  r=search(b);results.append(r)
  print(f"b={b}: {r['unsigned_triangle_count']} magnitudes, {r['primitive_triangle_candidates']} primitive triangles, cycles={[x['period'] for x in r['cycles']]}",flush=True)
 out={'method':'exact signature reconstruction with no parameter-numerator cutoff','results':results,'max_period_found':max((c['period'] for r in results for c in r['cycles']),default=0)}
 if args.output:open(args.output,'w').write(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'denominator_count':len(results),'max_period_found':out['max_period_found'],'total_cycles':sum(len(r['cycles']) for r in results)},sort_keys=True))
if __name__=='__main__':main()
