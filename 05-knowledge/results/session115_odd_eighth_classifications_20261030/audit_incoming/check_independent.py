import cmath, math, json
from collections import defaultdict

def ones(r,q,A):
    n=r; count=0
    for _ in range(A):
        if n%2: n=(q*n+1)//2; count+=1
        else: n//=2
    return count

def phase(x): return cmath.exp(2j*math.pi*x)
maxerr=0; universes=0; entries=0
for h in range(1,6):
 for q in (3,5,7):
  T=h+5; H=2**h
  for c in range(1,H,2):
   # Literal cost-truncated characters, grouped by frequency modulo 1/2^h.
   groups=defaultdict(lambda:[0j]*(T+1))
   for A in range(1,T+1):
    for r in range(1,2**A,2):
     m=ones(r,q,A)
     freq=r*2**(T-A)
     key=freq % 2**(T-h)
     groups[key][m]+=2**(-A)*phase(c*r/2**A)
   literal=[[sum(v[m]*v[n].conjugate() for v in groups.values()) for n in range(T+1)] for m in range(T+1)]
   F=[0j]*(T+1)
   for A in range(1,h+1):
    for r in range(1,2**A,2): F[ones(r,q,A)]+=2**(-A)*phase(c*r/2**A)
   predicted=[[F[m]*F[n].conjugate() for n in range(T+1)] for m in range(T+1)]
   for j in range(1,T-h+1):
    u=c*pow(pow(q,j,H),-1,H)%H
    coeff=[0j]*(h+1)
    for r in range(H): coeff[ones(r,q,h)]+=phase(u*r/H)
    weight=4**(-h)*sum(math.comb(L-1,j-1)*4**(-L) for L in range(j,T-h+1))
    for a in range(h+1):
     for b in range(h+1): predicted[j+a][j+b]+=weight*coeff[a]*coeff[b].conjugate()
   err=max(abs(literal[m][n]-predicted[m][n]) for m in range(1,T+1) for n in range(1,T+1))
   assert err<1e-11,(q,h,c,err)
   maxerr=max(maxerr,err); universes+=1; entries+=T*T
# Compute E3 coefficients from literal forward shortcut parity, not printed polynomial.
remainders=0
for c in (1,3,5,7):
 for M in range(3,15):
  K=0
  for b in range(3):
   u=c*pow(pow(3,M-b,8),-1,8)%8
   B=sum(phase(u*r/8) for r in range(8) if ones(r,3,3)<=b)
   K+=3**b*abs(B)**2/64
  chi=1 if c in (1,7) else -1
  target=.25+(-1)**M*3*math.sqrt(2)*chi/64
  assert abs(K-target)<1e-13
  decoded=88/105+1024*math.sqrt(2)/315*(-1)**M*(K-.25)
  assert abs(decoded-(1/c+1/(8-c)))<1e-13
  remainders+=1
print(json.dumps({'status':'PASS','finite_gram_universes':universes,'finite_gram_entries':entries,'h_range':[1,5],'q':[3,5,7],'cost_cutoff':'h+5','max_absolute_error':maxerr,'remainder_and_decoder_controls':remainders,'scope':'Floating point independent finite controls; symbolic all-depth proof audited separately. Not a convergence or novelty certificate.'},indent=2))
