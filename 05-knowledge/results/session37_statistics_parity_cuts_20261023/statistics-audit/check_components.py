#!/usr/bin/env python3
import json
from pathlib import Path
P=Path(__file__).resolve().parent

def check(b,m):
 if not b:raise AssertionError(m)
def mul(a,b):
 out={}
 for(n,k),x in a.items():
  for(m,l),y in b.items():
   if n+m<=7:out[n+m,k+l]=out.get((n+m,k+l),0)+x*y
 return out

def main():
 rows=json.loads((P/'graph_signs.json').read_text());chars=json.loads((P/'characters.json').read_text())['rows'];prod={(0,0):1}
 for r in rows:
  n=r['n']
  check(r['exterior_degrees']==chars[n]['exterior_invariant_degrees'],('exterior',n))
  check(r['royle_degrees']==chars[n]['walsh_invariant_degrees'],('Walsh',n))
  if n:
   b,f=r['bosonic_connected'],r['fermionic_connected'];e=r['royle_connected']
   check(b-f==(int(n==1)-int(n==3)),('balance',n))
   check(e==b+f-(rows[n//2]['fermionic_connected']if n%2==0 else 0),('dyadic',n))
  for k,count in enumerate(r['connected_exterior_degrees']):
   if not count:continue
   factor={(0,0):1}
   if k%2:factor[n,k]=1
   else:
    for j in range(1,7//n+1):factor[j*n,j*k]=1
   for _ in range(count):prod=mul(prod,factor)
 for r in rows:
  n=r['n'];check([prod.get((n,k),0)for k in range(n*(n-1)//2+1)]==r['exterior_degrees'],('bigraded product',n))
 # n3 grading cannot realize tournament converse on invariants.
 check(sum((-1)**k*a for k,a in enumerate(rows[3]['exterior_degrees']))==0,'exterior parity')
 check(sum((-1)**k*a for k,a in enumerate(rows[3]['royle_degrees']))==2,'tournament converse')
 receipt={'status':'PASS','graph_isomorphism_types_checked':sum(r['graphs']for r in rows),'max_n':7,'bigraded_product':True,'connected_balance_and_dyadic_recursion':True,'signed_stabilizer_matches_character':True,'n3_converse_intertwiner_obstruction':{'exterior_parity_trace':0,'tournament_converse_trace':2}}
 (P/'components.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
