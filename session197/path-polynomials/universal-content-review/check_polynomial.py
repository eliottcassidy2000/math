import json,math,random

def need(b,s):
 if not b:raise RuntimeError(s)
def mul(p,q):
 o=[0]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):o[i+j]+=x*y
 return o
def basis(k,parity):
 p=[0]*(2-parity)+[1]
 for j in range(1,k+1):p=mul(p,[-j*j,0,1])
 return p
def value(p,x):return sum(c*x**j for j,c in enumerate(p))
def newton(p):
 return [sum((-1)**(r-j)*math.comb(r,j)*value(p,j) for j in range(r+1)) for r in range(len(p))]
def gcd(xs):
 z=0
 for x in xs:z=math.gcd(z,x)
 return z
def coordinates(p,parity):
 p=list(p);d=(len(p)-1+1)//2;u=[0]*d
 for k in range(d-1,-1,-1):
  b=basis(k,parity);u[k]=p[len(b)-1]
  for j,x in enumerate(b):p[j]-=u[k]*x
 need(not any(p),'integer central basis exact');return u
def main():
 columns=entries=controls=0;R=random.Random(1971439)
 for parity in (0,1):
  for k in range(17):
   p=basis(k,parity);f=newton(p);gamma=math.factorial(2*k+1)*(1 if parity else k+1)
   for r,v in enumerate(f):
    ch=lambda a,b:math.comb(a,b) if 0<=b<=a else 0
    N=ch(k,2*k+1-r) if parity else ch(k+1,2*k+2-r)+ch(k,2*k+2-r)
    need(v==gamma*N,'central Newton binomial column')
    if r<=k:need(v==0,'triangular zeros')
    if r==k+1:need(v==gamma,'unit diagonal')
    if N:need(math.factorial(r-1)%math.factorial(k)==0,'weighted coefficient integer')
    entries+=1
   need(gamma%(math.factorial(k)**2)==0,'gamma double factorial divisibility');columns+=1
  for degree_index in range(1,12):
   deg=2*degree_index-parity
   for rep in range(100):
    p=[0]*(deg+1)
    for j in range(2-parity,deg+1,2):p[j]=R.randrange(-100,101)
    if p[-1]==0:p[-1]=1
    f=newton(p)[1:];u=coordinates(p,parity)
    gu=[math.factorial(2*k+1)*(1 if parity else k+1)*v for k,v in enumerate(u)]
    w=[math.factorial(r)*v for r,v in enumerate(f)]
    g,G=gcd(f),gcd(w)
    need(g==gcd(gu),'ordinary central content')
    need(G==gcd([math.factorial(k)*v for k,v in enumerate(gu)]),'universal central content')
    need(G==gcd(w[:degree_index]),'short Newton bank')
    need(G%g==0 and g**3%(G*G)==0,'sharp valuation upper bound')
    for prime in (2,3,5,7,11,13,17,19):
     need((g%prime==0)==(G%prime==0),'prime support')
     need((g%(prime**2)==0)==(G%(prime**2)==0),'square divisibility')
    controls+=1
 p=[0,0,210,0,214,0,34,0,1];u=coordinates(p,0);f=newton(p)[1:];w=[math.factorial(r)*v for r,v in enumerate(f)]
 need(gcd(f)==9 and gcd(w)==27 and gcd(f)**3==gcd(w)**2,'strict physical exponent equality')
 print(json.dumps({'status':'PASS','central_columns':columns,'exact_Newton_entries':entries,'signed_monomial_polynomials':controls,'strict_central_coordinates':u,'strict_g':9,'strict_G':27},sort_keys=True,indent=2))
if __name__=='__main__':main()
