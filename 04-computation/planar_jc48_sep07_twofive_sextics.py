#!/usr/bin/env python3
"""Exact normal-form boundaries and complete infinity-(2,13) A6 control.

Universal deformation/topology claims require the separate analytic proof.
"""
import sympy as S
from itertools import permutations
from collections import defaultdict, Counter
import hashlib,json
ncheck=0
def check(x,label):
 global ncheck
 ncheck+=1
 if not bool(x):raise RuntimeError(label)
def eq(a,b,label):check(S.cancel(a-b)==0,label)
z,a,b,c,d,e,t=S.symbols('z a b c d e t')
U=t**4+a*t**3+b*t**2
V=c*t**6+e*t**5+d*a*t**3+d*b*t**2
eq(S.expand(V-d*U+d*U**2/b**2).coeff(t,5),e+2*a*d/b,'finite fifth coefficient')
X=z*z*(1+a*z+b*z*z)/(c+e*z+d*a*z**3+d*b*z**4)
Z=z**6/(c+e*z+d*a*z**3+d*b*z**4)
eq(S.limit((Z-c*c*X**3)/z**7,z,0),(2*e-3*c*a)/c**2,'infinity seventh coefficient')
xn=S.cancel(X.subs({a:1,c:2,e:3}));zn=S.cancel(Z.subs({a:1,c:2,e:3}))
rho=S.cancel(zn/xn**3-4);k4=-6-24*b
tau=S.cancel(rho/xn-k4)
eq(S.limit(tau/z,z,0),7+12*b+8*d,'ninth coefficient supplier')
d0=-(7+12*b)/8
xx=S.cancel(xn.subs(d,d0));zz=S.cancel(zn.subs(d,d0));tt=S.cancel(tau.subs(d,d0))
k5=2*(3*b+1)*(20*b+3)
eq(S.limit(tt/xx,z,0),k5,'next even coefficient')
eta=S.cancel(tt/xx-k5)
eq(S.limit(eta/z,z,0),(3-20*b)/4,'eleventh coefficient supplier')
b0=S.Rational(3,20);x0=xx.subs(b,b0);z0=zz.subs(b,b0);et0=eta.subs(b,b0)
k6=-S.Rational(7069,250)
eq(S.limit(et0/x0,z,0),k6,'last even coefficient')
phi=S.cancel(et0/x0-k6)
eq(S.limit(phi/z,z,0),S.Rational(3,100),'thirteenth cannot vanish')
eq(S.limit((z0-4*x0**3-k4.subs(b,b0)*x0**4-k5.subs(b,b0)*x0**5-k6*x0**6)/z**13,z,0),S.Rational(3,6400),'actual odd thirteenth coefficient')
eq((e+2*a*d/b).subs({a:1,c:2,e:3,d:d0}),-S.Rational(7,4)/b,'finite cusp persists on eleventh family')
# U=t^4 case: V=c t^6+e t^5+f t^3+g t^2; e=0 gives odd infinity9.
f,g=S.symbols('f g')
vx=c*t**6+e*t**5+f*t**3+g*t*t
eq(S.expand(vx**2-g*g*t**4).coeff(t,5),2*f*g,'exceptional finite cusp criterion')
xc=z*z/(c+f*z**3+g*z**4);zc=z**6/(c+f*z**3+g*z**4)
eq(S.limit((zc-c*c*xc**3)/z**9,z,0),2*f/c**2,'exceptional infinity ninth coefficient')
for odd,nodes,margin in [(7,5,2),(9,4,0),(11,3,-2),(13,2,-4)]:
 check(2+(odd-1)//2+nodes==10,'complete genus ledger')
 check(16-(6+odd+1)==margin,'Nori linear margin')

weights={'L':1};edges=set()
centres=[('E1',('L',)),('E2',('L','E1')),('E3',('L','E2'))]
centres += [('E'+str(j),('E'+str(j-1),)) for j in range(4,8)]
centres += [('E8',('E6','E7'))]
for name,parents in centres:
 for p in parents:weights[p]-=1
 if len(parents)==2:edges.remove(tuple(sorted(parents)))
 weights[name]=-1
 for p in parents:edges.add(tuple(sorted((name,p))))
names=['L']+['E'+str(i) for i in range(1,9)]
check([weights[x] for x in names]==[-2,-2,-2,-2,-2,-2,-3,-2,-1],'thirteenth graph weights')
mat=S.Matrix([[-weights[x] if x==y else -int(tuple(sorted((x,y))) in edges) for y in names] for x in names])
check(mat.det()==-1,'total infinity determinant')
check(mat*S.Matrix([-6,-4,-8,-12,-10,-8,-6,-5,-10])==S.Matrix([0]*8+[1]),'marked H1 sign')
one=tuple(range(6))
def mul(p,q):return tuple(p[q[i]] for i in range(6))
def inv(p):return tuple(p.index(i) for i in range(6))
def pw(p,n):
 q=one
 for _ in range(n%60):q=mul(q,p)
 return q
def typ(p):
 seen=set();ls=[]
 for i in range(6):
  if i in seen:continue
  j=i;k=0
  while j not in seen:seen.add(j);k+=1;j=p[j]
  ls.append(k)
 return tuple(sorted(ls,reverse=True))
A6=[p for p in permutations(range(6)) if sum(p[i]>p[j] for i in range(6) for j in range(i+1,6))%2==0]
squares=defaultdict(list);cubes=defaultdict(list)
for p in A6:squares[pw(p,2)].append(p);cubes[pw(p,3)].append(p)
rows=[];fixed_counts=Counter();control=None
for c in sorted(set(squares)&set(cubes)):
 for l in squares[c]:
  for a in cubes[c]:
   bb=pw(a,2);d=mul(inv(mul(l,bb)),pw(c,2))
   if mul(c,d)!=mul(d,c):continue
   g=mul(inv(c),pw(d,2));h=mul(pw(c,-2),pw(d,3));f=mul(pw(c,-5),pw(d,7))
   for e in squares[f]:
    mu=mul(inv(mul(h,e)),f)
    if typ(mu)!=(3,1,1,1):continue
    fib=dict(zip(names,[l,a,bb,c,d,g,h,e,f]))
    inc={'L':[c],'E1':[bb],'E2':[a,c],'E3':[l,bb,d],'E4':[c,g],'E5':[d,h],'E6':[g,f],'E7':[f],'E8':[h,e,mu]}
    equations=[]
    for v in names:
     product=one
     for q in inc[v]:product=mul(product,q)
     equations.append(pw(fib[v],-weights[v])==product)
    equations += [mul(fib[x],fib[y])==mul(fib[y],fib[x]) for x,y in edges]
    equations += [mul(f,mu)==mul(mu,f)]
    if not all(equations):continue
    fixed=[i for i in range(6) if all(v[i]==i for v in fib.values()) and mu[i]==i]
    check(bool(fixed),'every marked image fixes a label')
    check(c==one or typ(c)==(5,1),'central-fibre analytic alternatives')
    check(f==one or (c==one and typ(f)==(3,1,1,1)),'terminal fibre control')
    fixed_counts[len(fixed)]+=1
    rows.append([fib,mu])
    if len(fixed)==1 and control is None:control=[fib,mu]
check(len(rows)==1120,'complete marked assignment count')
check(dict(fixed_counts)=={3:40,2:360,1:720},'complete support distribution')
check(control is not None,'nonabelian positive boundary control exists')
report={'status':'FINITE-EXACT normal forms and actual infinity13 marked group; global family transfer requires separate proof','nodes_by_infinity':{'7':5,'9':4,'11':3,'13':2},'normalized_even_coefficients':[str(k4),str(k5),str(k6)],'last_odd_coefficient':'3/6400','weights':weights,'edges':sorted(edges),'marked_assignments':len(rows),'fixed_label_counts':dict(fixed_counts),'positive_control':control,'gates':ncheck}
print(json.dumps(report,indent=2,sort_keys=True))
print('assignment_sha256='+hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest())
