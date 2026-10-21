#!/usr/bin/env python3
"""Independent exact standard-library audit; does not call the producer's code."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
from itertools import combinations, product
from math import factorial, gcd
from pathlib import Path
import json

SOURCE=Path(__file__).resolve().parent.parent/'session23b-flows-extra'
def need(ok,label):
    if not ok: raise RuntimeError(label)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def transpose(a):return [list(c)for c in zip(*a)]
def inverse_det(a):
    n=len(a);z=[[Q(x)for x in row]+[Q(i==j)for j in range(n)]for i,row in enumerate(a)];det=Q(1)
    for j in range(n):
        k=next((i for i in range(j,n)if z[i][j]),None)
        if k is None:return None,Q(0)
        if k!=j:z[j],z[k]=z[k],z[j];det=-det
        p=z[j][j];det*=p;z[j]=[x/p for x in z[j]]
        for i in range(n):
            if i!=j:
                p=z[i][j]
                if p:z[i]=[x-p*y for x,y in zip(z[i],z[j])]
    return [row[n:]for row in z],det

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def add(a,b):return trim([(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return trim(z)
def poly(exponents):
    z=[0]*(max(exponents)+1)
    for i,c in exponents.items():z[i]=c
    return trim(z)
def factor(d):return [1]+[0]*(d-1)+[-1]
def denominator(ds):return reduce(mul,(factor(d)for d in ds),[1])
def divide(a,b):
    a=list(map(Q,trim(a)));b=list(map(Q,trim(b)));c=[Q(0)]*max(1,(len(a)-len(b)+1))
    while a!=[0] and len(a)>=len(b):
        k=len(a)-len(b);v=a[-1]/b[-1];c[k]=v
        for j,y in enumerate(b):a[k+j]-=v*y
        a=trim(a)
    return trim(c),a

def pgcd(a,b):
    while b!=[0]:a,b=b,divide(a,b)[1]
    return scale(a,1/a[-1])
def reduced(n,d):
    g=pgcd(n,d);n,r=divide(n,g);need(r==[0],'numerator exact reduction');d,r=divide(d,g);need(r==[0],'denominator exact reduction')
    c=1/d[0];return scale(n,c),scale(d,c)
def gf(num,weights,N):
    z=[0]*(N+1)
    for i,v in enumerate(num):
        if i<=N:z[i]=v
    for w in weights:
        for i in range(w,N+1):z[i]+=z[i-w]
    return z

def plain(x):
    if isinstance(x,Q):return int(x)if x.denominator==1 else str(x)
    if isinstance(x,dict):return {str(k):plain(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return list(map(plain,x))
    return x

manifest=json.loads((SOURCE/'manifest.json').read_text())
for name,record in manifest['files'].items():
    b=(SOURCE/name).read_bytes();need(len(b)==record['bytes']and sha256(b).hexdigest()==record['sha256'],'manifest '+name)
need((SOURCE/'normal.json').read_bytes()==(SOURCE/'optimized.json').read_bytes(),'receipt equality')
J=json.loads((SOURCE/'normal.json').read_text())
need(J['status']=='PASS','producer status')
rays=[(1,0,0,0,0,0)]+[(1,)+tuple(int(j==i)for j in range(5))for i in range(5)]+[(1,)+tuple(int(j in(i,(i+2)%5))for j in range(5))for i in range(5)]+[(2,1,1,1,1,1)]
ineq=[(0,)+tuple(int(i==j)for j in range(5))for i in range(5)]+[(1,)+tuple(-int(j in(i,(i+1)%5))for j in range(5))for i in range(5)]
need([list(r)for r in rays]==J['triangulation']['rays'],'primitive ray list')
need(all(reduce(gcd,r)==1 for r in rays),'primitive rays')
need(all(dot(h,r)>=0 for h in ineq for r in rays),'ray feasibility')
# Independent Fraction-only exhaustive active-constraint test.
vertices=set();nonsingular=0
for active in combinations(range(10),5):
    inv,det=inverse_det([(1,0,0,0,0,0)]+[ineq[i]for i in active])
    if not det:continue
    nonsingular+=1;v=tuple(row[0]for row in inv)
    if all(dot(h,v)>=0 for h in ineq):vertices.add(v)
need(vertices=={tuple(Q(x,r[0])for x in r)for r in rays},'all 252 active sets')
need(nonsingular==J['vertex_certificate']['nonsingular_active_sets'],'nonsingular active sets')
# All primitive rays map to integral balanced flows; the proof below establishes saturation.
for t,*y in rays:
    x=[t-y[i]-y[(i-1)%5]for i in range(5)]
    need(all(z>=0 for z in x),'short arcs nonnegative')
    need(all(x[i]+y[i]==x[(i-1)%5]+y[(i-2)%5]for i in range(5)),'literal arc balance')
C=[tuple(v['rays'])for v in J['triangulation']['cells']]
need(len(C)==25 and len(set(C))==25 and all(len(c)==6 for c in C),'25 distinct six-ray cones')
Inv=[];dets=[];incidence=defaultdict(list)
for ci,c in enumerate(C):
    inv,det=inverse_det(transpose([rays[r]for r in c]));need(abs(det)==1,'unimodular cone '+str(ci));Inv.append(inv);dets.append(det)
    for j,r in enumerate(c):incidence[tuple(v for v in c if v!=r)].append((ci,j))
# Pairwise exact strict separators prove intersections equal their common faces.
separators=[]
for a,b in combinations(range(25),2):
    common=set(C[a])&set(C[b]);free=[i for i,r in enumerate(C[a])if r not in common]
    rows=[[dot(Inv[a][i],rays[r])for i in free]for r in C[b]if r not in common]
    weights=None
    for limit in (2,3,4,6,8):
        for v in product(range(1,limit+1),repeat=len(free)):
            if all(dot(row,v)<0 for row in rows):weights=v;break
        if weights is not None:break
    need(weights is not None,'exact strict separation '+str((a,b)))
    h=[sum(weights[k]*Inv[a][i][j]for k,i in enumerate(free))for j in range(6)]
    need(all(dot(h,rays[r])==0 for r in common),'separation equality')
    need(all(dot(h,rays[r])>0 for r in C[a]if r not in common),'first strict separation')
    need(all(dot(h,rays[r])<0 for r in C[b]if r not in common),'second strict separation')
    separators.append({'pair':[a,b],'common_rays':sorted(common),'normal':h})
# No internal boundary: every facet appears twice on opposite sides or is on dC.
facet_records=[];boundary=internal=0
for face,adj in sorted(incidence.items()):
    need(len(adj)in(1,2),'facet incidence at most two')
    a,j=adj[0];normal=Inv[a][j]
    need(all(dot(normal,rays[r])==0 for r in face),'facet vanishing')
    support=[i for i,h in enumerate(ineq)if all(dot(h,rays[r])==0 for r in face)]
    if len(adj)==1:
        need(bool(support),'no uncovered internal boundary');boundary+=1
    else:
        need(not support,'no duplicated exterior boundary')
        b,k=adj[1];need(dot(normal,rays[C[b][k]])<0,'opposite sides of internal facet');internal+=1
    facet_records.append({'rays':face,'incident_cells':[a for a,j in adj],'boundary_inequalities':support})
need(boundary+2*internal==150,'facet accounting')
# w is strictly interior, avoids every facet hyperplane, and yields integral shifts.
w=[sum((1<<i)*r[j]for i,r in enumerate(rays))for j in range(6)]
need(all(dot(h,w)>0 for h in ineq),'interior vector')
weights=[5*r[0]-sum(r[1:])for r in rays];groups=defaultdict(Counter);opens=[]
for ci,c in enumerate(C):
    alpha=[dot(row,w)for row in Inv[ci]];need(all(alpha),'generic vector')
    op=tuple(r for r,a in zip(c,alpha)if a<0);shift=sum(weights[r]for r in op);opens.append(op)
    groups[tuple(sorted(weights[r]for r in c))][shift]+=1
    saved=J['triangulation']['cells'][ci]
    need(saved['determinant']==dets[ci]and saved['offset_coordinates']==alpha and saved['open_ray_indices']==list(op)and saved['weight_shift']==shift,'record matches independently computed cell')
A=poly(dict(groups[(3,3,4,4,5,5)]));B=poly(dict(groups[(3,3,3,4,5,5)]))
need(sum(groups[(3,3,4,4,5,5)].values())==15 and sum(groups[(3,3,3,4,5,5)].values())==10,'groups')
P=poly({0:1,3:2,4:3,5:1,6:3,7:1,8:3,9:2,12:1})
N=poly({0:1,4:1,5:1,6:1,7:1,8:1,12:1})
need(add(mul(A,factor(3)),mul(B,factor(4)))==mul(P,factor(5)),'all-order numerator identity')
DR=denominator((3,3,3,4,4,5));DK=denominator((3,5,6,8))
need(pgcd(P,DR)==[1]and pgcd(N,DK)==[1],'reduced raw and fixed fractions')
# Independently verify reduced orbit/merged denominators, recurrence orders, periods.
FN,FD=reduced(add(mul(P,factor(5)),scale(denominator((3,3,3,4,4)),4)),scale(mul(DR,factor(5)),5))
HN,HD=reduced(add(mul(FN,DK),mul(N,FD)),scale(mul(FD,DK),2))
phis={1:[-1,1],2:[1,1],3:[1,1,1],4:[1,0,1],5:[1,1,1,1,1],6:[1,-1,1],8:[1,0,0,0,1]}
profiles={}
for name,d in [('R',DR),('f',FD),('k',DK),('merged',HD)]:
    left=d;profile={}
    for i,p in phis.items():
        multiplicity=0
        while True:
            quotient,rem=divide(left,p)
            if rem!=[0]:break
            multiplicity+=1;left=quotient
        if multiplicity:profile[i]=multiplicity
    need(len(left)==1,'complete cyclotomic factorization')
    profiles[name]={'pole_orders':profile,'minimal_recurrence_order':len(d)-1,'numerator_degree':len({'R':P,'f':FN,'k':N,'merged':HN}[name])-1,'minimal_quasiperiod':reduce(lambda a,b:a*b//gcd(a,b),profile,1)}
need([profiles[n]['minimal_recurrence_order']for n in('R','f','k','merged')]==[22,26,22,32],'orders')
# Complete direct integer enumeration, plus exhaustive half-open ownership at low weight.
M=45;bound=25;raw=[0]*(M+1);ownership_points=0;boundary_points=0
for t in range(2*M//5+1):
 for a in range(t+1):
  for b in range(t-a+1):
   for c in range(t-b+1):
    for d in range(t-c+1):
     for e in range(min(t-d,t-a)+1):
      m=5*t-a-b-c-d-e
      if m>M:continue
      raw[m]+=1
      if m>bound:continue
      z=(t,a,b,c,d,e);owners=[];on_boundary=False
      for ci,cell in enumerate(C):
       lam=[dot(row,z)for row in Inv[ci]]
       if any(x<0 for x in lam):continue
       on_boundary|=any(x==0 for x in lam)
       if all(v>0 if r in opens[ci]else v>=0 for r,v in zip(cell,lam)):owners.append(ci)
      need(len(owners)==1,'unique half-open lattice owner '+str(z));ownership_points+=1;boundary_points+=int(on_boundary)
need(raw==gf(P,(3,3,3,4,4,5),M),'raw enumeration through 45')
need(raw[:41]==J['independent_counts']['raw_coefficients'],'producer raw controls')
KMAX=100;kc=[0]*(KMAX+1)
for t in range(2*KMAX//5+1):
 for b in range(t//2+1):
  for a in range(t-b+1):
   for c in range(t-a+1):
    m=5*t-2*a-2*b-c
    if m<=KMAX:kc[m]+=1
need(kc==gf(N,(3,5,6,8),KMAX),'independent supplied fixed formula through 100')
fc=[];hc=[]
for m,r in enumerate(raw):
    total=r+4*(m//5+1 if m%5==0 else 0);need(total%5==0,'ordinary Burnside integral');f=total//5
    need(f>=kc[m]and (f+kc[m])%2==0,'anti Burnside integral');fc.append(f);hc.append((f+kc[m])//2)
Rlead=Q(sum(P),factorial(5)*3**3*4**2*5);Klead=Q(sum(N),factorial(3)*3*5*6*8)
need(Rlead==Q(17,259200)and Klead==Q(7,4320),'growth constants')
need(Klead/(Rlead/5)==Q(2100,17),'fixed fraction asymptotic')
# Exact expanded denominators are recurrence certificates, not fitted recurrences.
for coefficients,num,den in [(raw,P,DR),(fc,FN,FD),(kc,N,DK),(hc,HN,HD)]:
    convolution=mul(coefficients,den)
    need(convolution[:len(coefficients)]==(num+[0]*len(coefficients))[:len(coefficients)],'exact recurrence control')
out={'status':'PASS','source_manifest_verified':manifest,'independent_method':'Fraction Gaussian elimination, exact pairwise strict separators, facet incidence/side certificate, dense rational polynomial arithmetic, bounded direct enumeration; no producer execution or triangulation recursion imported','vertices':{'active_sets':252,'nonsingular_sets':nonsingular,'vertices':[list(v)for v in sorted(vertices)]},'cones':{'count':25,'determinants':dets,'pairwise_separators':separators,'facets':facet_records,'boundary_facets':boundary,'internal_facets':internal,'generic_vector':w,'grouped_shifts':[{'weights':k,'shifts':dict(v)}for k,v in sorted(groups.items())]},'controls':{'raw_bound':M,'raw_coefficients':raw,'ownership_bound':bound,'ownership_points':ownership_points,'boundary_points':boundary_points,'fixed_bound':KMAX,'fixed_coefficients':kc,'ordinary_coefficients':fc,'merged_coefficients':hc},'rational_consequences':{'profiles':profiles,'recurrence_denominators':{'R':DR,'f':FD,'k':DK,'merged':HD},'raw_leading':Rlead,'ordinary_leading':Rlead/5,'fixed_leading':Klead,'fixed_fraction':Klead/(Rlead/5)}}
print(json.dumps(plain(out),sort_keys=True,indent=2))
