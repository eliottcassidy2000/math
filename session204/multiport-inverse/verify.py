import sympy as S,json
D,u,v=S.symbols('D u v');X=D+u;Y=D+v
moves=[(2*D+u,u,(D+7*u)/6-36),(2*D+v,2*u-v-36,(D+13*u)/6-v-78)]
def need(x,msg):
 if not x:raise RuntimeError(msg)
for b,(d,a,c) in enumerate(moves):
 need(S.expand(d+a-(2*X-36*b))==0,'X quotient')
 need(S.expand(d+c-(S.Rational(13,6)*X-36-42*b))==0,'Y quotient')
K,b=S.symbols('K b');need(S.expand((72*K-36*b)-(78*K-36-42*b)-36-(-6*K+6*b))==0,'u increment')
count=0;todo=[(563,-23,-17,15,0,[],[15])]
while todo:
 d,a,c,k,t,w,ks=todo.pop();n=len(w)
 need(d==36*k+23+6*t and a==-23-6*t,'prefix identity')
 need(k==15*2**n-int(''.join(map(str,w))or'0',2),'binary identity')
 if n:
  need(0<=t<=15*2**(n-1)-8,'T bound')
  need(504*2**n+59<=d<=585*2**n-25,'band')
 count+=1
 if n==12:continue
 for b in (0,1):
  if b==0:dd,aa,cc=2*d+a,a,S.Rational(d+7*a,6)-36
  else:dd,aa,cc=2*d+c,2*a-c-36,S.Rational(d+13*a,6)-c-78
  inc=(7 if n==0 else ks[-2]-w[-1])*b
  todo.append((dd,aa,cc,2*k-b,t+inc,w+[b],ks+[2*k-b]))
need(count==8191,'enumeration coverage')
print(json.dumps(dict(status='PASS',exact_symbolic_identities=5,finite_regression_words=count,maximum_length=12,scope='All-word algebraic prefix formula and dyadic gaps for fixed family only.')))
