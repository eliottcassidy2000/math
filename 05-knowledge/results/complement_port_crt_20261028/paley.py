"""Actual Paley atom actions, marked multiplier obstruction, and moment recursions."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from math import factorial,prod
from itertools import product
import random,json
import ports
HERE=Path(__file__).resolve().parent
CHECKS=0

def need(v,m):
 global CHECKS
 CHECKS+=1
 if not v:raise RuntimeError(m)

def squares(p):return sorted({a*a%p for a in range(1,p)})
def affine(p):return [tuple((a*x+b)%p for x in range(p)) for a in squares(p) for b in range(p)]
def atom(p):
 H=squares(p)
 return tuple(sum(1<<y for y in range(p) if (y-x)%p in H) for x in range(p))
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inverse(a):return tuple(a.index(i) for i in range(len(a)))
def root(a):
 v=tuple(range(len(a)));order=0
 while True:
  v=compose(v,a);order+=1
  if v==tuple(range(len(a))):break
 need(order%2==1,'unique odd square root')
 v=tuple(range(len(a)))
 for _ in range((order+1)//2):v=compose(v,a)
 need(compose(v,v)==a,'literal square root')
 return v
def neg(p,k):
 if not k:return (0,)
 m=p**(k-1);t=neg(p,k-1)
 return tuple((-a%p)*m+t[j] for a in range(p) for j in range(m))
def lift(p,a,b,children):
 m=len(children[0]);return tuple(((a*x+b)%p)*m+children[x][j] for x in range(p) for j in range(m))
def random_auto(p,k,rng):
 if not k:return (0,)
 return lift(p,rng.choice(squares(p)),rng.randrange(p),[random_auto(p,k-1,rng) for _ in range(p)])
def random_central(p,k,rng):
 if not k:return (0,)
 t=neg(p,k-1);child=[None]*p;child[0]=random_central(p,k-1,rng)
 for x in squares(p):
  child[x]=random_auto(p,k-1,rng);child[-x%p]=compose(t,compose(child[x],t))
 return lift(p,rng.choice(squares(p)),0,child)
def characters(g,p,k):
 if not k:return []
 m=p**(k-1);b=g[0]//m;a=(g[m]//m-b)%p
 child=[characters(tuple(g[x*m+j]%m for j in range(m)),p,k-1) for x in range(p)]
 return [a]+[prod(c[d] for c in child)%p for d in range(k-1)]
def local_scalar(p,k,depth,a):
 if depth==0:return lift(p,a,0,[tuple(range(p**(k-1)))]*p)
 children=[tuple(range(p**(k-1))) for _ in range(p)]
 children[0]=local_scalar(p,k-1,depth-1,a)
 return lift(p,1,0,children)
def mul(a,b):
 return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(len(a))]
def power(a,n):
 out=[Q(1)]+[Q(0)]*(len(a)-1)
 for _ in range(n):out=mul(out,a)
 return out
def moments(p,depth=8,R=6):
 h=(p-1)//2;f=[Q(1,factorial(r)) for r in range(R+1)];g=f[:];rows=[]
 for k in range(depth+1):
  aa=[f[r]*factorial(r) for r in range(R+1)];cc=[g[r]*factorial(r) for r in range(R+1)]
  need(all(v.denominator==1 for v in aa+cc),'integer Burnside moments')
  need(aa[1:4]==[1,2*k+1,6*k*k+2*(p-2)*k+1],'Paley first three moments')
  need(cc[1:3]==[2*k+1,8*k*k+4*k*(h-1)+1],'Paley centralizer moments')
  rows.append({'depth':k,'automorphism_moments':list(map(int,aa)),'centralizer_moments':list(map(int,cc))})
  linked=power([f[r]*2**r for r in range(R+1)],h)
  linked=[x/h for x in linked];linked[0]+=Q(h-1,h)
  g=mul(g,linked)
  fp=power(f,p)
  f=[Q(p-3,p-1)*x+Q(2,p*(p-1))*y for x,y in zip(f,fp)];f[0]+=Q(2,p)
 return rows

def main(write=True):
 global CHECKS
 CHECKS=0;out=[];rng=random.Random(71112)
 for p in (3,7,11):
  T=atom(p);A=affine(p);t=neg(p,1);h=(p-1)//2
  full=ports.oracle.isomorphisms(T,T,True)
  need(set(A)==set(full),'full automorphism group at verified Paley base')
  need(ports.core.is_strong(T),'strong Paley atom')
  # Test all nontrivial candidate modules, independently of decomposition code.
  for mask in range(1,(1<<p)-1):
   if mask.bit_count()==1:continue
   members=[x for x in range(p) if mask>>x&1]
   need(any(len({(T[y]>>x)&1 for x in members})>1 for y in range(p) if not(mask>>y&1)),'indecomposable Paley base')
  C=[g for g in A if compose(g,t)==compose(t,g)]
  fixed=Counter(sum(g[i]==i for i in range(p)) for g in A)
  expected=Counter({0:p-1,1:p*(h-1),p:1});expected+=Counter()
  need(fixed==expected,'Paley offspring law from actual automorphisms')
  marks=[compose(g,compose(t,inverse(g))) for g in A];marks=sorted(set(marks))
  loops=set()
  for s in marks:
   for u in marks:
    bs=root(compose(s,t));bu=root(compose(u,t));edge=root(compose(u,s))
    loops.add(compose(inverse(bu),compose(edge,bs)))
  need(loops=={tuple(range(p))},'base Paley marking connection is flat')
  for k in (1,2,3):
   tt=neg(p,k)
   for trial in range(5):
    g=random_auto(p,k,rng);v=random_auto(p,k,rng);c=random_central(p,k,rng)
    need(characters(compose(g,v),p,k)==[a*b%p for a,b in zip(characters(g,p,k),characters(v,p,k))],'depth multipliers are homomorphisms')
    need(characters(compose(tt,compose(g,tt)),p,k)==characters(g,p,k),'negation preserves multiplier characters')
    need(compose(c,tt)==compose(tt,c),'constructive centralizer element')
   for d in range(k):
    for a in squares(p):
     c=local_scalar(p,k,d,a);target=[1]*k;target[d]=a
     need(characters(c,p,k)==target and compose(c,tt)==compose(tt,c),'centralizer surjects onto independent depth multipliers')
  mm=moments(p)
  for r in range(7):
   need(Q(sum(sum(g[i]==i for i in range(p))**r for g in A),len(A))==mm[1]['automorphism_moments'][r],'base all-order moments')
   need(Q(sum(sum(g[i]==i for i in range(p))**r for g in C),len(C))==mm[1]['centralizer_moments'][r],'base centralizer moments')
  out.append({'prime':p,'atom':T,'full_base_aut_order':len(A),'base_centralizer_order':len(C),
              'base_triangle_errors':len(loops),'fixed_count_distribution':sorted(fixed.items()),
              'proved_holonomy_index_depth_1_to_8':[h**k for k in range(1,9)],
              'proved_holonomy_orders_depth_1_to_3':[(p*h)**(((p**k-1)//(p-1)-k)//2) for k in range(1,4)],
              'moments':mm})
 result={'status':'PROVED recurrences; FINITE-EXACT literal base and character controls','rows':out,'checks':CHECKS}
 if write:(HERE/'paley.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'base_controls':[(r['prime'],r['full_base_aut_order'],r['base_centralizer_order']) for r in out]}))
 return result
if __name__=='__main__':main()
