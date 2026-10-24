import json
from math import comb

def req(ok,msg):
    if not ok:raise RuntimeError(msg)
N=10000
c=[1,2]
for n in range(1,N):
    q,r=divmod((26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1],(n+1)**2)
    req(r==0,'nonintegral recurrence');c.append(q)
for n in range(241):
    F=sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2,n+1))
    req(F==sum(c[k]*c[n-k] for k in range(n+1)),'independent binomial square')
def prediction(n):
    return 1 if n==0 else 2 if n&(n-1)==0 else 0
req(all(v%4==prediction(n) for n,v in enumerate(c)),'dyadic mod4 support')
req(all(v%2==0 for v in c[1:]),'positive parity')
req(all((n*c[n])%2==0 for n in range(N+1)),'Euler derivative even')
# Expand finite product directly, allowing all cross terms before reduction.
p=[1]+[0]*N;power=1
while power<=N:
    for n in range(N,power-1,-1):p[n]=(p[n]+2*p[n-power])%4
    power*=2
req(all(v==c[n]%4 for n,v in enumerate(p)),'finite product')
# Actual naive-lift relation coefficientwise, not general composition equivalence.
for n in range(N+1):
    rhs=(c[n//2] if n%2==0 else 2*c[(n-1)//2])%4
    req(c[n]%4==rhs,'naive lift after parity')
# A=1+t has nonzero discrepancy2t³; its Euler derivative is not even.
a=[1,0,1,2]; naive=[1,0,1,0]
req(any((x-y)%4 for x,y in zip(a,naive)),'naive lift countermodel')
req(c[0]==1 and c[0]%4!=0,'n0 exception')
req(c[1]%4==2,'not all divisible by4')
req(c[2]%8==6 and c[2]%8!=2,'unchanged mod8 product rejected')
req((c[4]-c[1])%8==4,'mod8 two-step hostile')
# Mutation of a power-of-two coefficient must violate the result under -O too.
req((c[8]+2)%4!=prediction(8),'mutation detected')
print(json.dumps({'scope':'FINITE-EXACT support of separately proved corollary','coefficient_bound':N,'independent_binomial_square_bound':240,'positive_nonzero_mod4_indices':[n for n in range(1,N+1) if c[n]%4],'countermodels':{'naive_lift_A_1_plus_t_difference_at_t3':2,'constant_at_n0':1,'c2_mod8':c[2]%8,'c4_minus_c1_mod8':(c[4]-c[1])%8},'status':'PASS'},sort_keys=True,indent=2))
