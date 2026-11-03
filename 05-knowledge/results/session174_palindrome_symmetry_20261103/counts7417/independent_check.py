#!/usr/bin/env python3
import json, math
from pathlib import Path
p=7417
def ck(x):
 if not x: raise RuntimeError('audit failure')
cert=json.loads(Path(__file__).with_name('certificate.json').read_text())
rows=[]
for m in range(1,7):
 counts=[0]*p
 for half in range(10**m,10**(m+1)):
  s=str(half);counts[int(s+s[-2::-1])%p]+=1
 row={'digits':2*m+1,'total':sum(counts),'occupied':sum(x>0 for x in counts),'minimum':min(counts),'maximum':max(counts),'divisible_by_7417':counts[0]}
 ck(row==cert['palindrome_censuses'][m-1]); rows.append(row)
units={}
for n,expected in [(240,4),(1001,60),(2222,100)]:
 orders=[]
 for a in range(n):
  if math.gcd(a,n)!=1:continue
  x=a%n;k=1
  while x!=1:x=x*a%n;k+=1
  orders.append(k)
 ck(math.lcm(*orders)==expected)
 units[n]={'order':len(orders),'exponent':expected}
x=1; orbit=[]
while x not in orbit:orbit.append(x);x=x*10%p
ck(len(orbit)==2472 and x==1 and orbit[1236]==p-1)
ck(next(k for k in range(1,2473) if (1+pow(10,k,p))%p==0)==1236)
ck(next(r for r in range(1,2473) if (pow(10,r,p)+pow(10,-r,p))%p==0)==618)
print(json.dumps({'status':'PASS','independent_method':'direct half-word enumeration through 13 digits, not DP','rows':rows,'unit_groups':units,'first_zero_pair_distance':618,'first_zero_exponent_gap':1236},indent=2))
