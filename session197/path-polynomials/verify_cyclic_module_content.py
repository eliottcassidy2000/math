#!/usr/bin/env python3
from verify_universal_module_content import *
from verify_cyclic_gram_forests import jacobi

def kernel(c):
 f=[factorial(j)*v for j,v in enumerate(c)];f[0]=0
 J=jacobi(f)
 return [[factorial(i+1)*factorial(j+1)*x for j,x in enumerate(row)]for i,row in enumerate(J)]
def kgcd(K):return reduce(gcd,[x for row in K for x in row])
def formal_kernel(P):
 f=newton(P);J=jacobi(f)
 return [[factorial(i+1)*factorial(j+1)*x for j,x in enumerate(row)]for i,row in enumerate(J)]
def dot(K,x,y):return sum(K[i][j]*v*w for i,v in enumerate(x)for j,w in enumerate(y))
def run():
 physical=bank_counts=formal=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   V=tournament(n,code);c=cover_polynomial(V);P=bin_polynomial(V);u=coordinates(P,n%2);K=kernel(c);g,G,_,_=contents(c)
   a=[factorial(k)**2*factorial(2*k+1)*(1 if n%2 else k+1)*v for k,v in enumerate(u)];C=reduce(gcd,a)
   require(C==kgcd(K)and C%G==0 and c[1]%C==0 and G*G%(g*C)==0 and g*g%C==0 and G**4%(C**3)==0,'physical cyclic content and valuation bounds');physical+=1
   L=max(2,(n+1)//2-1);vectors=[cover_polynomial(transitive(j))[1:]for j in range(1,L+1)]
   vals=[dot(K,x,y)for x in vectors for y in vectors];require(reduce(gcd,vals)==C,'short physical square bank gcd')
   if n<=4:
    for i in range(1,L+1):
     for j in range(1,L+1):require(h_subsets(cyclic(transitive(i),transitive(j),V))[-1]==vals[(i-1)*L+j-1],'actual cyclic bank count');bank_counts+=1
 rng=Random(1971432)
 for parity in (0,1):
  for d in range(1,9):
   for turn in range(40):
    u=[rng.randrange(-9,10)for k in range(d)];u[0]=2*rng.randrange(-8,9)+1
    P=add(*(scale(basis(k,parity),v)for k,v in enumerate(u)));K=formal_kernel(P)
    a=[factorial(k)**2*factorial(2*k+1)*(1 if parity else k+1)*v for k,v in enumerate(u)]
    require(kgcd(K)==gcd(*a),'formal odd-evaluation cyclic content')
    x,y=u[0],(factorial(3)*(1 if parity else 2)*u[1]if d>1 else 0)
    require((K[0][0],K[0][1])==((3*x+y,2*x+2*y)if parity else(5*x+y,6*x+2*y)),'exact determinant-four initial boundary');formal+=1
 hostile=[0,1,0,1];require(kgcd(formal_kernel(hostile))==4,'even-evaluation hostile kernel')
 V=(382,420,88,706,810,972,18,597,204,327);validate(V);c=cover_polynomial(V);P=bin_polynomial(V);g,G,f,w=contents(c);K=kernel(c);C=kgcd(K)
 require(c[1:]==[7857,70551,185274,211989,124725,40851,7672,814,45,1]and(g,G,C)==(9,27,81),'literal physical three-level hierarchy')
 require(P==[0,0,2772,0,4020,0,1000,0,64,0,1]and coordinates(P,0)==[7857,10449,2043,94,1],'literal full bin profile')
 graphs=[cyclic(transitive(1),transitive(1),V),cyclic(transitive(1),transitive(2),V),substitute(bank(4),0,V)]
 counts=[h_subsets(a)[-1]for a in graphs]
 require(counts==[164673,462591,94405041],'three literal large context counts')
 require(gcd(*counts[:2])==81 and gcd(c[1],counts[2])==27 and counts[2]%81==27,'strict physical restriction')
 require(exterior_runs(bank(4),0)==[91,101,44,6],'four-run matching context coefficients')
 print(json.dumps({'status':'PASS','physical_profiles':physical,'literal_small_bank_counts':bank_counts,'signed_formal_odd_boundary_controls':formal,'even_boundary_hostile_content':4,'witness':{'input_rows':V,'input_covers':c[1:],'bin_profile':P,'central_coordinates':coordinates(P,0),'ordinary_content':g,'universal_content':G,'cyclic_content':C,'kernel_diagonal':[K[i][i]for i in range(len(K))],'kernel_off_diagonal':[K[i][i+1]for i in range(len(K)-1)],'graphs':[{'rows':a,'H':h}for a,h in zip(graphs,counts)]}},indent=2,sort_keys=True))
if __name__=='__main__':run()
