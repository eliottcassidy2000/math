#!/usr/bin/env python3
"""Independent final-extension checks; no producer imports. Requires SymPy."""
import itertools,json,math,hashlib
from pathlib import Path
import sympy as s
P=Path(__file__).resolve().parent; W=P.parent

def need(t,m):
 if not t:raise RuntimeError(m)
def arrays(n):
 Q={i*i%n for i in range(n)}; out=set()
 for c,u,v in itertools.product(range(n),repeat=3):
  a=tuple(x%n for x in (c+u,c-u-v,c+v,c-u+v,c,c+u-v,c-v,c+u+v,c-u))
  if all(x in Q for x in a):out.add(a)
 return out
pairs=list(itertools.combinations(range(9),2))
def mask(a):return sum((a[i]==a[j])<<k for k,(i,j) in enumerate(pairs))
T={n:arrays(n) for n in (2,4,8,16,32,64)}
need(T[4]=={(0,)*9,(1,)*9},'M4')
need(T[8]=={tuple(4*x for x in a) for a in T[2]}|{(1,)*9},'M8')
need(T[16]=={tuple(1+8*x for x in a) for a in T[2]}|{(0,)*9,(4,)*9},'M16')
need({mask(a) for a in T[2]}=={mask(a) for a in T[8]}=={mask(a) for a in T[16]},'full mask equality')
need(all(len(set(a))<9 for n in (32,64) for a in T[n]),'32,64 Parker')
w=(41,6,31,16,26,36,21,46,11); roots=(29,34,41,46,26,44,39,36,31)
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
need(len(set(w))==9 and all(r*r%50==x for r,x in zip(roots,w)),'50 roots')
need(all(sum(w[i] for i in L)%50==28 for L in lines),'50 magic')
need(len({x%2 for x in w})==2 and all(len({x%2 for x in a})==1 for a in T[4]),'no 100 lift')
# Literal polynomial rational functions, independently reduced by SymPy.
z=s.symbols('z')
L1=lambda t:t/(1-t)**2
L2=lambda t:t*(1+t)/(1-t)**3
L3=lambda t:t*(1+4*t+t*t)/(1-t)**4
F=L3(z)-8*L2(z)+7*L1(z)-16*L2(z*z)+26*L1(z*z)+24*L1(z**3)+48*L1(z**4)+40*L1(z**5)+48*L1(z**6)
num,den=s.fraction(s.cancel(F)); orders={1:4,2:3,3:2,4:2,5:2,6:2}
D=s.prod(s.cyclotomic_poly(d,z)**k for d,k in orders.items())
need(s.expand(den-D)==0,'exact reduced denominator')
Q=[27,37,119,201,346,438,639,649,729,635,576,400,300,160,94,34,14,2]
need(s.expand(num-8*z**9*sum(v*z**i for i,v in enumerate(Q)))==0,'numerator')
need(s.gcd(num,den)==1,'coprime numerator denominator')
need(math.lcm(*orders)==60 and s.rem(den,s.cyclotomic_poly(60,z),z)!=0,'period and absent 60 pole')
def C(n):
 g=lambda d:math.gcd(n,d)
 return 4*g(4)+4*g(6)+2*g(2)**2+2*g(5)+g(2)*g(4)-2-4*g(2)
need(all(C(n)==7+13*(n%2==0)+8*(n%3==0)+12*(n%4==0)+8*(n%5==0)+8*(n%6==0) for n in range(1,121)),'indicator')
# Recompute Smith aggregation rather than trusting source count coefficients.
rows=((2,0),(0,2),(2,2),(2,-2),(1,2),(1,-2),(2,1),(2,-1)); agg={}
for k in range(2,9):
 for ss in itertools.combinations(rows,k):
  d1=math.gcd(*(x for r in ss for x in r)); det=math.gcd(*(a*d-b*c for (a,b),(c,d) in itertools.combinations(ss,2)))
  key=(d1,det//d1);agg[key]=agg.get(key,0)+(-1)**k
agg={k:v for k,v in agg.items() if v}
need(agg=={(1,1):-2,(1,2):-4,(1,4):4,(1,5):2,(1,6):4,(2,2):2,(2,4):1},'Smith inclusion exclusion')
# Diagram affine coordinates are tested as equivalence relations, not colours.
diagram_dir=W/'session55-diagrams'
if not diagram_dir.exists():diagram_dir=W/'diagrams'
c=json.loads((diagram_dir/'countdata.json').read_text())['crt']
for vals,coords in ((c['mod5'],[(i-j)%3 for i in range(3) for j in range(3)]),(c['mod7'],[(i+j+1)%3 for i in range(3) for j in range(3)])):
 need(all((vals[i]==vals[j])==(coords[i]==coords[j]) for i,j in pairs),'affine partition attribution')
need(c['ordinary_row_sums']==[10,45,80] and c['modular_line_sum']==10,'counts vs integer contrast')
out={'status':'PASS','local_arrays':{n:len(v) for n,v in T.items()},'local_masks':{n:len({mask(a) for a in v}) for n,v in T.items()},'rational_numerator_degree':s.degree(num,z).item() if hasattr(s.degree(num,z),'item') else int(s.degree(num,z)),'denominator_degree':int(s.degree(den,z)),'minimal_period':60,'smith_coefficients':{str(k):v for k,v in agg.items()},'sympy_version':s.__version__}
print(json.dumps(out,indent=2,sort_keys=True))
