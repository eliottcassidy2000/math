#!/usr/bin/env python3
"""Independent rooted permutation DFS for every certificate, plus inverse tails.
All checks are explicit exceptions and remain active under python -O.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import csv,json,hashlib
HERE=Path(__file__).resolve().parent

def check(t,msg):
 if not t: raise RuntimeError(msg)

def read(n,kind):
 return {int(x['value']):(list(map(int,x['rows'].split(','))),int(x['root'])) for x in csv.DictReader(open(HERE/f'n{n}_{kind}.tsv'),delimiter='\t')}

def tournament(rows):
 n=len(rows)
 check(all(0<=r<(1<<n) and not(r>>i&1) for i,r in enumerate(rows)),'diagonal/range')
 check(all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i,j in combinations(range(n),2)),'pair orientation')

def rooted_paths(rows,root,closing=False):
 n=len(rows);full=(1<<n)-1
 if closing and n<3:return 0
 def dfs(S,u):
  if S==full:return int(bool(rows[u]>>root&1)) if closing else 1
  out=rows[u]&~S;total=0
  while out:
   bit=out&-out;out-=bit;total+=dfs(S|bit,bit.bit_length()-1)
  return total
 return dfs(1<<root,root)

def converse(rows):
 return [sum(1<<j for j in range(len(rows)) if rows[j]>>i&1) for i in range(len(rows))]

def close_endpoint(rows,r):
 n=len(rows)
 return [a|((i==r)<<n) for i,a in enumerate(rows)]+[((1<<n)-1)^(1<<r)]

def pad_last(rows,extra):
 # Last vertex must have indegree one. Replace it by a transitive chain.
 n=len(rows);v=n-1;ins=[i for i in range(v) if rows[i]>>v&1]
 check(len(ins)==1,'padding vertex is not indegree one')
 out=[a&((1<<v)-1) for a in rows[:v]]
 for u in range(v):
  if rows[u]>>v&1:out[u]|=((1<<(extra+1))-1)<<v
 for i in range(extra+1):
  out.append((rows[v]&((1<<v)-1))|sum(1<<(v+j) for j in range(i+1,extra+1)))
 return out

def inherited_comparison(cycle_supports):
 # Separate inherited raw class table, not an input to the generator/census.
 import gzip
 expected_counts={6:56,7:456,8:6880,9:191536}
 for n in range(6,10):
  path=HERE/(f'classes{n}.tsv'+('.gz' if n==9 else ''))
  if not path.exists():continue
  f=gzip.open(path,'rt') if n==9 else path.open()
  rows=list(csv.DictReader(f,delimiter='\t'));vals=set()
  check(len(rows)==expected_counts[n],'inherited count mismatch')
  for x in rows:
   numerator=int(x['closed'])*int(x['Aut'])
   check(numerator%n==0,'nonintegral cycle quotient')
   vals.add(numerator//n)
  check(vals==cycle_supports[n],'independent inherited spectrum mismatch')


def main():
 Cs={};Es={};checks=0
 for n in range(1,10):
  Cs[n]=read(n,'cycles');Es[n]=read(n,'endpoints')
  for c,(rows,r) in Cs[n].items():
   tournament(rows);check(rooted_paths(rows,0,True)==c,f'cycle certificate {n} {c}');checks+=1
  for e,(rows,r) in Es[n].items():
   tournament(rows);check(rooted_paths(converse(rows),r)==e,f'endpoint certificate {n} {e}');checks+=1
 check(set(Cs[7])==set(range(16))|{17,24},'order7 spectrum')
 check(set(Cs[8])==set(range(55)),'order8 spectrum')
 holes9={215,219,223,226,229,*range(232,243),*range(244,249)}
 check(set(Cs[9])==set(range(250))-holes9,'order9 spectrum')
 check(set(range(594))<=set(Es[9]) and 594 not in Es[9],'endpoint contiguous prefix')
 for n in range(3,10):check(set(Cs[n])<=set(Es[n]),f'finite C_n subset E_n fails at {n}')
 inherited_comparison({n:set(v) for n,v in Cs.items()})
 inverses=[];tail_checks=0
 for k in range(1,594):
  candidates=[n for n in range(1,10) if k in Cs[n]]
  threshold=min(candidates) if candidates else 10
  # Choose the earliest known endpoint certificate. Its closed witness and all
  # transitive-chain paddings furnish every size at least tail_start.
  m=min(n for n in range(1,10) if k in Es[n] and n>=2)
  tail_start=m+1
  check(all(k in Cs[n] for n in range(threshold,min(tail_start,10))),f'gap before certified tail {k}')
  rows,r=Es[m][k];closed=close_endpoint(rows,r)
  for extra in (0,1,2):
   graph=closed if not extra else pad_last(closed,extra)
   tournament(graph);check(rooted_paths(graph,0,True)==k,f'physical closure/padding {k} +{extra}')
   tail_checks+=1
  inverses.append({'k':k,'exact_min_order':threshold,'all_larger_orders':True,'endpoint_witness_order':m,'endpoint_root':r,'endpoint_rows':rows})
 (HERE/'inverse_sizes_1_593.json').write_text(json.dumps(inverses,indent=2)+'\n')
 summary={'verified_saved_witnesses':checks,'physical_endpoint_closure_and_padding_checks':tail_checks,'exact_inverse_counts':593,'endpoint_prefix_through_order9':593,'cycle_supports':{n:sorted(c) for n,c in Cs.items()},'order9_cycle_holes_below_max':sorted(holes9),'boundaries':'positive tournament orders; directed Hamiltonian cycles modulo cyclic rotation, not reversal; HC=0 for orders1,2'}
 (HERE/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
