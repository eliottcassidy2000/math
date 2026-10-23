"""Small exact controls for native odd-depth target coverage bounds."""
from collections import Counter
import json

def need(c,msg):
 if not c:raise ArithmeticError(msg)
def step(n,s):
 z=3*n+s;a=(z&-z).bit_length()-1
 return z//2**a,a
rows=[];checks=Counter()
for X in (1,2,16,64,256,1024):
 for r in range(5):
  bounds=[(3*(2**i)*X+1).bit_length()-1 for i in range(r)]
  prod=bound=1
  for a in bounds:prod*=a;bound+=prod
  coarse=sum((X.bit_length()-1+r+3)**j for j in range(r+1))
  need(bound<=coarse,'coarse bound')
  for s in (-1,1):
   for y in (1,5,17):
    hits=[];word_owners={}
    for n in range(1,X+1,2):
     z=n;w=[];reaches=z==y
     if reaches:word_owners[()]=n
     for i in range(r):
      z,a=step(z,s);w.append(a)
      need(a<=bounds[i],'valuation bound')
      need(z<=2**(i+1)*X,'height bound')
      checks['step_bounds']+=1
      if z==y:
       need(tuple(w) not in word_owners or word_owners[tuple(w)]==n,'word target injectivity')
       word_owners[tuple(w)]=n;reaches=True
     if reaches:hits.append(n)
    need(len(hits)<=bound,'coverage bound')
    checks['coverage_cases']+=1
    if X==1024 and r==4:rows.append({'X':X,'depth':r,'sign':s,'target':y,'actual_starts':len(hits),'product_bound':bound,'coarse_bound':coarse})
print(json.dumps({'counts':dict(checks),'largest_controls':rows},indent=2))
