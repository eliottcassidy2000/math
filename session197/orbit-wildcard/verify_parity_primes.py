import itertools as it, math, json
from verify import gate,beta,clone_count,runs,cyclic

def valuation(n,p):
 e=0
 while n%p==0:n//=p;e+=1
 return e
checks=0
for D in range(1,6):
 delta=[None]+[math.prod(b**e for b in range(1,D+1))*math.prod((j-i)*(j+i) for i in range(1,D+1) for j in range(i+1,D+1)) for e in (1,2)]
 for coeff in it.product(range(3),repeat=D):
  if not coeff[0] or not coeff[-1]:continue
  for scale in (1,7,289):
   a=[0]+[scale*x for x in coeff];b=beta(a);g=math.gcd(*a[1:])
   gate(math.gcd(*b[1:])==g,'unimodular coefficient gcd')
   gall=math.gcd(*(math.factorial(k)*a[k] for k in range(1,D+1)))
   gate(gall==math.gcd(*(clone_count(a,j) for j in range(1,D+1))),'unconditioned Stirling gcd')
   gs=[]
   for e in (1,2):
    ge=math.gcd(*(clone_count(a,e+2*t) for t in range(D)));gs.append(ge)
    gate(ge%g==0 and g*delta[e]%ge==0,'integral determinant bound')
    gate(ge%gall==0 and delta[e]%(ge//gall)==0,'relative parity divisor bound')
    for p in (2,3,5,7,11,13,17,19):
     if p>2*D-1:gate(valuation(ge,p)==valuation(g,p),'large prime valuation')
   gate(math.gcd(*gs)==gall,'union parity gcd');checks+=1
physical=[]
for q in (3,5,7):
 rows=cyclic(q,list(range(1,(q+1)//2)));a=runs(rows,0);D=len(a)-1
 O=[j for j in range(1,q) if rows[0]>>j&1];I=[j for j in range(1,q) if rows[j]&1]
 best=0
 def match(t,used,count):
  global best
  if t==len(O):best=max(best,count);return
  match(t+1,used,count)
  for j in I:
   if not(used>>j&1) and rows[O[t]]>>j&1:match(t+1,used|(1<<j),count+1)
 match(0,0,0);gate(D==best+1,'physical matching degree')
 physical.append({'q':q,'D':D,'matching':best,'prime_cutoff':2*D-1,'largest_certificate_order':q-1+2*D})
print(json.dumps({'status':'PASS','formal_scaled_vectors':checks,'physical_templates':physical},sort_keys=True,indent=2))
