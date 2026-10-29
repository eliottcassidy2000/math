"""Exhaustive exact set cover, restricted to recorded primes <=71."""
from pathlib import Path
from itertools import combinations
import json
ROOT=Path(__file__).resolve().parent

def main():
 c=json.loads((ROOT/'certificate.json').read_text());pairs=list(combinations(range(len(c['ids'])),2));full=(1<<len(pairs))-1
 masks={int(p):sum(1<<k for k,(i,j) in enumerate(pairs) if row[i] is not None and row[j] is not None and row[i]!=row[j]) for p,row in c['traces'].items()}
 primes=sorted(p for p,m in masks.items() if m);tested=0;found=None
 for size in range(len(primes)+1):
  for subset in combinations(primes,size):
   tested+=1;mask=0
   for p in subset:mask|=masks[p]
   if mask==full:found=subset;break
  if found is not None:break
 if found is None:raise ValueError('bank fails cover')
 forced={}
 for k,(i,j) in enumerate(pairs):
  only=[p for p,m in masks.items() if (m>>k)&1]
  if len(only)==1:forced.setdefault(only[0],[c['ids'][i],c['ids'][j]])
 if set(forced)!=set(found):raise ValueError('forced-prime lower bound does not match optimum')
 result={'forced_prime_witness_pairs':forced,'scope':'Minimum cardinality only among recorded primes <=71; only common good unequal traces count.','minimum_cardinality':len(found),'one_minimum_bank':list(found),'subsets_tested_including_solution':tested,'zero_coverage_primes':[p for p in sorted(masks) if not masks[p]],'pairs':len(pairs)}
 print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
