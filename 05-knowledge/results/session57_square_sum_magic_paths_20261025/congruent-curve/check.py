from fractions import Fraction as F
from math import isqrt,lcm
from itertools import combinations
import json
from pathlib import Path

def square(q):
 q=F(q)
 return q>=0 and isqrt(q.numerator)**2==q.numerator and isqrt(q.denominator)**2==q.denominator

def root(q):
 assert square(q)
 q=F(q);return F(isqrt(q.numerator),isqrt(q.denominator))

def add(P,Q,N):
 if P is None:return Q
 if Q is None:return P
 x,y=P;u,v=Q
 if x==u and y==-v:return None
 m=(3*x*x-N*N)/(2*y) if P==Q else (v-y)/(u-x)
 z=m*m-x-u
 return z,m*(x-z)-y

def triangle(P,N):
 x,y=P;assert y and y*y==x**3-N*N*x
 return tuple(abs(z/y) for z in (x*x-N*N,2*N*x,x*x+N*N))

def point(a,b,h,N):return (N*(a+h)/b,2*N*N*(a+h)/(b*b))

def integer_centers(d):
 out=[]
 if d%2:return out
 for s in range(1,isqrt(d//2)+1):
  if d%(2*s):continue
  r=d//(2*s)
  if r>s and square(r*r+s*s):out.append((r*r+s*s,r,s))
 return sorted(out)

def aps(vals):
 vals=sorted(set(vals));S=set(vals)
 return [(a,b,2*b-a) for a,b in combinations(vals,2) if 2*b-a in S]

def magic(C,d,e):
 return [[C+d,C-d-e,C+e],[C-d+e,C,C+d-e],[C-e,C+d+e,C-d]]

def guarded(C,d,e):return C>d+e and d>0 and e>0 and len({C+i*d+j*e for i in (-1,0,1) for j in (-1,0,1)})==9

N=F(6);d=4*N;P=(F(12),F(36));Q=None;rows=[];centers=[]
for n in range(1,25):
 Q=add(Q,P,N);a,b,h=triangle(Q,N);c=h*h
 assert a*a+b*b==h*h and a*b/2==N
 assert square(c-d) and square(c+d)
 assert c==4*add(Q,Q,N)[0]
 u,w=root(c-d),root(c+d);aa,bb=(w+u)/2,(w-u)/2
 assert aa*bb/2==N and aa*aa+bb*bb==c
 assert triangle(point(aa,bb,h,N),N)==(aa,bb,h)
 assert add(point(aa,bb,h,N),point(aa,bb,h,N),N)==(c/4,h*u*w/8)
 rows.append({'n':n,'P':[str(t) for t in Q],'triangle':[str(t) for t in (a,b,h)],'center':str(c)})
 centers.append(c)
assert centers[0]==25
assert centers[1]==F(1442401,4900)
assert centers[0]+centers[2]!=2*centers[1]
# Translation by 2-torsion does not create a distinct doubled x coordinate.
for T in [(F(0),F(0)),(N,F(0)),(-N,F(0))]:
 assert add(add(P,T,N),add(P,T,N),N)==add(P,P,N)
assert integer_centers(24)==[(25,4,3)]
# Complete integer-center catalogs for each fixed integer d <= 20000.
cat={d:integer_centers(d) for d in range(1,20001)}
cat={d:rows for d,rows in cat.items() if rows}
integer_ap=[]
for d,rs in cat.items():
 for a,b,c in aps([r[0] for r in rs]):
  if guarded(b,d,b-a):integer_ap.append((d,a,b,c))
# Independently check factor catalog against direct root enumeration to h=200.
brute={}
for h in range(1,201):
 for u in range(h):
  d=h*h-u*u
  if 0<d<=20000 and square(h*h+d):brute.setdefault(d,set()).add(h*h)
for d in range(1,20001):
 assert {c for c,r,s in cat.get(d,[]) if c<=40000}==brute.get(d,set())
# Pure linear magic and collision guards, without claiming square entries.
for d in range(1,15):
 for e in range(1,15):
  C=d+e+1;M=magic(C,d,e)
  sums=[sum(r) for r in M]+[sum(M[i][j] for i in range(3)) for j in range(3)]+[sum(M[i][i] for i in range(3)),sum(M[i][2-i] for i in range(3))]
  assert sums==[3*C]*8
  assert guarded(C,d,e)==(e not in (d,2*d,F(d,2)))
# Common root-denominator scaling preserves a common area across triangles while multiplying it by L²; individual primitive reduction does not.
L=lcm(*[z.denominator for row in rows[:3] for z in map(F,row['triangle'])])
for row in rows[:3]:
 a,b,h=map(F,row['triangle']);assert (L*a).denominator==(L*b).denominator==(L*h).denominator==1
 assert L*a*L*b/2==L*L*N
out={'status':'PASS','scope':'24 positive multiples of (12,36) on E_6; integer d from 1 through 20000 only; no global exclusion', 'multiples':rows,'doubled_x_center_APs':[[str(z) for z in r] for r in aps(centers)],'integer_center_catalog':{str(d):r for d,r in cat.items()},'integer_APs_with_nine_distinct_guard':integer_ap,'hostiles':{'group_AP_not_center_AP':str(centers[0]-2*centers[1]+centers[2]),'torsion_x_AP_rejected':['(-6,0)','(0,0)','(6,0)'],'two_torsion_translates_duplicate_doubles':True,'individual_scale_not_common_area':True},'common_triangle_denominator_first_three':L}
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('multiples','integer_center_catalog')},indent=2))
print('integer gaps with centers:',len(cat),'centers:',sum(map(len,cat.values())),'max per gap:',max(map(len,cat.values())))
