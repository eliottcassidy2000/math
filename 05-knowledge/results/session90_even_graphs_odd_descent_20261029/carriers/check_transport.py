#!/usr/bin/env python3
"""Exact anti-involution transport, mass and nonflat nine-vertex control.
Self-contained, standard library; active checks survive python -O.
Default exhausts paired classes n<=3 and direct permutation Aut searches.
"""
from itertools import permutations, product
from collections import Counter
from fractions import Fraction
from math import factorial
import json
CHECKS=0

def need(v,msg='check failed'):
 global CHECKS
 CHECKS+=1
 if not v:raise AssertionError(msg)

def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):
 b=[0]*len(a)
 for i,j in enumerate(a):b[j]=i
 return tuple(b)
def odd_sqrt(a):
 """O(N) array operations, accepts precisely permutations with odd cycles."""
 seen=[False]*len(a)
 for v in a:
  if not isinstance(v,int) or not 0<=v<len(a) or seen[v]:raise ValueError('not a permutation')
  seen[v]=True
 out=[-1]*len(a)
 for i in range(len(a)):
  if out[i]>=0:continue
  cyc=[i];j=a[i]
  while j!=i:cyc.append(j);j=a[j]
  L=len(cyc)
  if L%2==0:raise ValueError('even cycle')
  for k,v in enumerate(cyc):out[v]=cyc[(k+(L+1)//2)%L]
 return tuple(out)
def conj(g,a):return compose(compose(g,a),inverse(g))
def transport(t,s):return odd_sqrt(compose(t,s))
def holonomy(u,t,s):
 return compose(inverse(transport(u,s)),compose(transport(u,t),transport(t,s)))
def acts(rows,p,anti=False):
 n=len(rows)
 return all(((rows[i]>>j)&1)==((rows[p[j] if anti else p[i]]>>(p[i] if anti else p[j]))&1) for i in range(n) for j in range(i+1,n))
def relabel(rows,p):
 out=[0]*len(rows)
 for i,row in enumerate(rows):
  for j in range(len(rows)):
   if row>>j&1:out[p[i]]|=1<<p[j]
 return tuple(out)
def lift(m,n):
 rows=[0]*(2*n+1)
 def arc(i,j,b):rows[i if b else j]|=1<<(j if b else i)
 for i in range(n):
  arc(1+i,1+n+i,1);d=(m>>(i*n+i))&1
  arc(0,1+i,d);arc(0,1+n+i,1-d)
 for i in range(n):
  for j in range(i+1,n):
   x=m>>(i*n+j)&1;y=m>>(j*n+i)&1;r=x^y
   arc(1+i,1+j,x);arc(1+n+i,1+n+j,1-x)
   arc(1+i,1+n+j,r);arc(1+j,1+n+i,r)
 return tuple(rows),(0,)+tuple(range(n+1,2*n+1))+tuple(range(1,n+1))
def paired_action(m,p):
 n=len(p);out=0
 for i in range(n):out|=((m>>(i*n+i))&1)<<(p[i]*n+p[i])
 for i in range(n):
  for j in range(i+1,n):
   x=m>>(i*n+j)&1;y=m>>(j*n+i)&1
   a,b=p[i],p[j]
   if a>b:a,b=b,a;x^=1;y^=1
   out|=x<<(a*n+b);out|=y<<(b*n+a)
 return out

def automorphisms(rows):
 # Exact full S_N search, degree filtering only prunes impossible candidates.
 degrees=tuple(r.bit_count() for r in rows)
 return tuple(p for p in permutations(range(len(rows))) if all(degrees[i]==degrees[p[i]] for i in range(len(rows))) and acts(rows,p))

def check_object(rows,tau,A,full_triangles=True):
 N=len(rows);e=tuple(range(N));As=set(A)
 need(e in As and len(A)%2==1,'odd group')
 need(acts(rows,tau,True) and compose(tau,tau)==e,'anti involution')
 C=tuple(a for a in A if compose(a,tau)==compose(tau,a))
 I=tuple(sorted({conj(a,tau) for a in A}))
 # Independently enumerate the anti coset and filter squares.
 need(set(I)=={compose(a,tau) for a in A if compose(compose(a,tau),compose(a,tau))==e},'anti list')
 need(len(I)*len(C)==len(A),'centralizer mass')
 K={compose(s,tau) for s in I}
 need(K=={a for a in A if conj(tau,a)==inverse(a)},'twisted coordinates')
 need(len(K)==len(I),'K/I bijection')
 need(all(odd_sqrt(k) in K for k in K),'square-root closure')
 need({compose(k,c) for k in K for c in C}==As,'polar factorization')
 pairs=0;nontrivial=0;witness=None
 for s,t in product(I,repeat=2):
  g=transport(t,s);need(g in As,'transport is ordinary auto')
  need(conj(g,s)==t,'transport conjugates')
  need(transport(s,t)==inverse(g),'inverse transport')
  need(odd_sqrt(compose(g,g))==g,'unique odd root')
  pairs+=1
 for s,t,u in product(I,repeat=3):
  if not full_triangles and s!=tau:continue
  h=holonomy(u,t,s)
  need(h in As and compose(h,s)==compose(s,h),'centralizer holonomy')
  if h!=e:
   nontrivial+=1
   if witness is None:witness={'sigma':s,'tau':t,'upsilon':u,'holonomy':h}
 commute=all(compose(x,y)==compose(y,x) for x,y in product(K,repeat=2))
 need((nontrivial==0)==commute,'flat iff K commutes')
 if commute:
  need(all(compose(x,y) in K for x,y in product(K,repeat=2)),'flat K subgroup')
 return {'vertices':N,'aut':len(A),'centralizer':len(C),'anti_involutions':len(I),'transport_pairs':pairs,'nontrivial_triangles':nontrivial,'witness':witness},I,K

def run():
 global CHECKS
 reports=[];all_classes=0;transport_pairs=0;mass=[]
 for n in range(4):
  ps=tuple(permutations(range(n)));seen=set();hist=Counter();weighted=Fraction();ordinary=Fraction();tri=0
  for m in range(1<<(n*n)):
   key=min(paired_action(m,p) for p in ps)
   if key in seen:continue
   seen.add(key);W,tau=lift(key,n);A=automorphisms(W)
   r,I,K=check_object(W,tau,A)
   need(r['nontrivial_triangles']==0,'all odd SC <=7 are flat')
   all_classes+=1;transport_pairs+=r['transport_pairs'];hist[(r['aut'],r['centralizer'],r['anti_involutions'])]+=1
   weighted+=Fraction(len(I),len(A));ordinary+=Fraction(1,len(A))
   # Full isomorphism-to-marked correction after every vertex relabel through N=5.
   if n<=2:
    for p in permutations(range(len(W))):
     Wp=relabel(W,p);sigma=conj(p,tau)
     for t in I:
      target=conj(p,t);g=transport(target,sigma);fixed=compose(g,p)
      need(g==conj(p,transport(t,tau)),'transport equivariant under arbitrary relabeling')
      need(relabel(W,fixed)==Wp,'corrected isomorphism preserves target')
      need(conj(fixed,tau)==target,'corrected isomorphism intertwines marks')
   tri+=len(I)**3
  need(len(seen)==(1,2,8,88)[n],'pair counts')
  need(weighted==Fraction(1<<(n*n),factorial(n)),'marked mass formula')
  mass.append({'n':n,'classes':len(seen),'marked_mass':str(weighted),'ordinary_mass':str(ordinary),'triangles':tri,'symmetry_profiles':[{'aut':a,'centralizer':c,'anti':i,'classes':v} for (a,c,i),v in sorted(hist.items())]})
 # Canonical C3[C3,C3,C3], independently full S9 search, then compare wreath group.
 W=tuple(sum(1<<(3*k+l) for k in range(3) for l in range(3) if ((k-i)%3==1 if i!=k else (l-j)%3==1)) for i in range(3) for j in range(3))
 tau=tuple(3*((-i)%3)+(-j)%3 for i in range(3) for j in range(3))
 explicit={tuple(3*((i+k)%3)+(j+b[i])%3 for i in range(3) for j in range(3)) for k in range(3) for b in product(range(3),repeat=3)}
 A=automorphisms(W);need(set(A)==explicit,'full S9 group equals wreath81')
 r,I,K=check_object(W,tau,A);need((r['aut'],r['centralizer'],r['anti_involutions'])==(81,3,27),'wreath profile')
 need(r['nontrivial_triangles']>0,'wreath nonflat')
 b={s:transport(s,tau) for s in I}
 based={(t,s):compose(b[t],inverse(b[s])) for s,t in product(I,repeat=2)}
 for s,t,u in product(I,repeat=3):
  need(compose(based[u,t],based[t,s])==based[u,s],'based transport is flat')
 for s,t in product(I,repeat=2):
  defect=compose(inverse(b[t]),compose(transport(t,s),b[s]))
  need(compose(defect,tau)==compose(tau,defect),'based defect in centralizer')
  need(compose(b[t],compose(defect,inverse(b[s])))==transport(t,s),'exact based reconstruction')
 K=tuple(sorted(K));ks=set(K);e=tuple(range(9))
 need(all(compose(compose(x,y),x) in ks for x,y in product(K,repeat=2)),'twisted subgroup')
 op=lambda x,y:odd_sqrt(compose(compose(x,compose(y,y)),x))
 table={ (x,y):op(x,y) for x,y in product(K,repeat=2)}
 need(all(z in ks for z in table.values()),'loop closed')
 need(all(len({table[x,y] for y in K})==len(K) and len({table[y,x] for y in K})==len(K) for x in K),'latin loop')
 need(all(table[e,x]==x and table[x,e]==x for x in K),'identity')
 need(all(inverse(table[x,y])==table[inverse(x),inverse(y)] for x,y in product(K,repeat=2)),'AIP')
 assoc_witness=None
 for x,y,z in product(K,repeat=3):
  need(table[x,table[y,table[x,z]]]==table[table[x,table[y,x]],z],'left Bol')
  if assoc_witness is None and table[table[x,y],z]!=table[x,table[y,z]]:assoc_witness={'x':x,'y':y,'z':z}
 need(assoc_witness is not None,'nonassociative Bruck loop27')
 r['bruck_order']=len(K);r['bruck_nonassociative_witness']=assoc_witness
 # Hostile: even cycle square root cannot be supplied by odd-cycle routine.
 try:odd_sqrt((1,0))
 except ValueError:need(True)
 else:need(False,'even root guard absent')
 return {'status':'PASS','active_checks':CHECKS,'scope':'all natural paired classes n<=3; full ordinary automorphism permutation search through order7 and wreath9; all triangle holonomies for each','small_mass_profiles':mass,'small_classes':all_classes,'small_transport_pairs':transport_pairs,'wreath9':r}
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
