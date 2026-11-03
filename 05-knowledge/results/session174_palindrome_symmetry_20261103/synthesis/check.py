import json
checks=0
def check(c):
 global checks
 checks+=1
 if not c:raise RuntimeError('synthesis check failed')
def valuation(n,p):
 k=0
 while n%p==0:n//=p;k+=1
 return k
inert=[5,7,13,23,29,31,37,47]
for q in inert:check(pow(q-2,(q-1)//2,q)==q-1)
for a in range(1,101):
 for b in range(-100,101):
  H=3*a*a-2*a*b+b*b;x=b-a
  check(H==x*x+2*a*a)
  for q in inert:check(valuation(H,q)%2==0)
  if a%2 and b%2==0:check(H%8==3)
  check(H not in (7,21))
solutions=[(x,a) for x in range(49) for a in range(49) if (x*x+2*a*a)%49==0]
check(set(solutions)=={(x,a) for x in range(0,49,7) for a in range(0,49,7)})
check(3*9**2-2*9*84+84**2==5787)
print(json.dumps({'status':'PASS','explicit_checks':checks,'a_range':[1,100],'b_range':[-100,100],'inert_primes':inert,'norm_zero_pairs_mod49':len(solutions)},indent=2))
