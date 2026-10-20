#!/usr/bin/env python3
import sympy as sp,json,itertools
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
checks=0
def req(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def a(n):
 if n%3==0:return -8*n**4+16*n**5+32*n**6-64*n**7
 if n%3==1:return 8*n**4+16*n**5-96*n**6-64*n**7
 return -n+12*n*n-60*n**3+160*n**4-240*n**5+192*n**6-64*n**7
x,z=sp.symbols('x z');polys=[-8*x**4+16*x**5+32*x**6-64*x**7,8*x**4+16*x**5-96*x**6-64*x**7,-x+12*x*x-60*x**3+160*x**4-240*x**5+192*x**6-64*x**7]
chi=sp.Poly((x-1)**8*(x*x+x+1)**7,x);coeff=list(reversed([int(v) for v in chi.all_coeffs()]))
for r in range(3):req(sp.expand(sum(coeff[j]*polys[(r+j)%3].subs(x,3*z+r+j) for j in range(23)))==0,'global symbolic recurrence')
H=sp.Matrix(22,22,lambda i,j:a(i+j));smith=smith_normal_form(H,domain=ZZ);diag=[abs(int(smith[i,i])) for i in range(22)]
req(all(diag[i+1]%diag[i]==0 for i in range(21)),'Smith divisibility')
req(H.det()==-2**206*3**191*5**22*7**8*61**7,'determinant')
req(sp.prod(diag)==abs(H.det()),'Smith determinant')
C=sp.zeros(22)
for j in range(21):C[j+1,j]=1
for i in range(22):C[i,21]=-coeff[i]
req(H*C==C.T*H,'companion intertwining')
def vp(n,p,K=1000):
 if n==0:return K
 v=0
 while n%p==0:n//=p;v+=1
 return v
def local_smith(p,K,ncols=22):
 mod=p**K;M=[[int(H[i,j])%mod for j in range(ncols)] for i in range(22)];vals=[]
 for k in range(ncols):
  v,i,j=min((vp(M[i][j],p,K),i,j) for i in range(k,22) for j in range(k,ncols))
  req(v<K,'local pivot nonzero')
  M[k],M[i]=M[i],M[k]
  for row in M:row[k],row[j]=row[j],row[k]
  power=p**v;unit=M[k][k]//power;inv=pow(unit,-1,p**(K-v))
  for i in range(k+1,22):
   factor=(M[i][k]//power)*inv%(p**(K-v))
   M[i]=[(M[i][j]-factor*M[k][j])%mod for j in range(ncols)]
   req(M[i][k]==0,'row elimination')
  for j in range(k+1,ncols):
   factor=(M[k][j]//power)*inv%(p**(K-v))
   for i in range(22):M[i][j]=(M[i][j]-factor*M[i][k])%mod
   req(M[k][j]==0,'column elimination')
  vals.append(v)
 return vals
vals={}
for p in [2,3,5,7,61]:
 expected=[vp(d,p) for d in diag];actual=local_smith(p,30);req(actual==expected,'independent local Smith');vals[p]=actual
prefixes={}
for r in [6,12,14,18,20,22]:
 actual=local_smith(2,30,r);req(actual==vals[2][:r],'prefix image certificate');prefixes[r]=actual
H6=H[:6,:6];b=sp.Matrix([a(i+6) for i in range(6)]);rec=[]
for k in range(1,8):
 mod=2**k;v=H6.inv_mod(mod)*b;c=[int(q%mod) for q in v]
 residual=[(a(n+6)-sum(c[j]*a(n+j) for j in range(6)))%mod for n in range(22)]
 req((not any(residual))==(k<=6),'precision threshold')
 rec.append({'k':k,'coefficients':c,'nonzero_initial_residuals':[[i,v] for i,v in enumerate(residual) if v]})
for k in range(1,13):
 mod=2**k;P=3*mod;seq=[a(n)%mod for n in range(P)]
 req(all(a(n+P)%mod==seq[n] for n in range(P)),'period upper control')
 for d in sp.divisors(P):
  if d<P:req(any(seq[n]!=seq[(n+d)%P] for n in range(P)),'proper divisor hostile')
seen=set();hist={}
for v in itertools.product(range(2),repeat=6):
 if v in seen:continue
 orbit=[];w=v
 while w not in seen:seen.add(w);orbit.append(w);w=w[1:]+w[:1]
 req(w==v,'rotation cycle');hist[len(orbit)]=hist.get(len(orbit),0)+1
req(hist=={1:2,2:1,3:2,6:9},'binary cycle inventory')
print(json.dumps({'status':'PASS','checks':checks,'smith_diagonal':diag,'prime_valuations':vals,'prefix_two_adic_valuations':prefixes,'six_registers':rec,'binary_cycles':hist,'scope':'Exact symbolic global recurrence, exact integer Smith plus independent prime-adic pivot, finite period/cycle controls; general clock proven in note'},indent=2))
