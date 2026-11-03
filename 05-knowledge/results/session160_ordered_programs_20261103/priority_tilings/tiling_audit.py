from fractions import Fraction
from math import lcm
from itertools import product,combinations
P=list(map(Fraction,['17/143','133/85','17/19','23/17','1015/69','23/29','31/23','38073/62','31/37','41/31','86/287','41/43','13/41']))
D=[p.denominator for p in P]
def req(condition,msg):
 if not condition: raise RuntimeError(msg)
def statuses(n):return tuple(int(n%d==0) for d in D)
def feasible(bits):
 L=lcm(*(d for d,b in zip(D,bits) if b));return all(b or L%d for d,b in zip(D,bits))
def scan(bits):
 v=(Fraction(1),Fraction(0))
 for e,f in zip(bits,P):
  v=(v[0]*(1-e), v[0]*e*f+v[1])
 return v
def nextn(n):
 for f in P:
  out=n*f
  if out.denominator==1:return out.numerator
 return None
ct=both=0
for e in product((0,1),repeat=13):
 test=all(e[i]<=e[j] for i,j in [(1,3),(4,6),(7,9),(10,12)])
 req(feasible(e)==test,('criterion',e))
 a=next((i for i,b in enumerate(e) if b),None)
 req(scan(e)==((1,0) if a is None else (0,P[a])),('scan',e))
 ct+=test
 both+=test and feasible(tuple(1-b for b in e))
req((ct,both)==(2592,512),('counts',ct,both))
req(statuses(23)==statuses(46) and nextn(23)==31 and nextn(46)==62 and statuses(31)!=statuses(62),'collision')
req(not feasible(tuple(1-b for b in statuses(41))),'complement')
ch=[(0,2),(0,3),(1,3)]
def adj(bits):
 b=dict(zip(ch,bits))
 return {(i,j) if j==i+1 or b[(i,j)] else (j,i) for i,j in combinations(range(4),2)}
a,b=adj((0,0,0)),adj((0,0,1));p=(3,0,1,2)
req({(p[i],p[j]) for i,j in a}==b,'isomorphism')
print('PASS 8192 patterns, exact matrix output, all finite image counts and collisions')
