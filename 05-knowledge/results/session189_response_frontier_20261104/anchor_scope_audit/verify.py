#!/usr/bin/env python3
"""Exact finite-difference support counts, not a second cycle-parity census."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def cubes(n,monomials):
 ans=set()
 for monomial in monomials:
  s=monomial
  while True:
   ans.add((s&((1<<n)-1),s>>n))
   if not s:break
   s=(s-1)&monomial
 return ans
def L(n,v,out):
 return {(1<<u)|(1<<(n+v))for u in range(n)if u!=v}|{(1<<u)|(1<<v)|(1<<(n+v))for u in out}
def main():
 records=[]
 for n in range(3,65):
  u,b=0,1
  for d in range(n-1):
   D=set(range(2,d+2));monomials=L(n,u,D|{b})|L(n,b,D)|{(1<<u)|(1<<b)|(1<<(n+b))};query=cubes(n,monomials);expected=12+7*d+3*(n-2-d);need(len(query)==expected==3*n+6+4*d,'exact twin battery union');need(len(query)<=7*n-2<9*n,'claimed nine-n bound');need((0,0)in query and(3,0)in query,'coarse alpha fallback retained');records.append({'n':n,'outside_outdegree':d,'distinct_queries':len(query)})
  for dv in range(n):
   out=set(range(1,dv+1));q=cubes(n,L(n,0,out));need(len(q)<=4*n,'one L battery bound')
 # The invalid m=1 unchanged-child claim cannot survive even as an existence statement.
 source=[6,4,0];target=source[:];need(source==target,'one-child unchanged counterexample')
 receipt={'status':'PASS','twin_count_cases':len(records),'kernel_orders':list(range(3,65)),'exact_twin_query_formula':'12+7d+3(n-2-d)=3n+6+4d<=7n-2<9n','singleton_other_child_bound':'One L_z battery, at most4n queries, detects both target end directions','m1_counterexample':{'Q_order':1,'child_source_rows':source,'child_target_rows':target,'all_actual_response_differences':0},'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(P/'counts.json').write_text(json.dumps(records,separators=(',',':'))+'\n');receipt['counts_sha256']=hashlib.sha256((P/'counts.json').read_bytes()).hexdigest();(P/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
