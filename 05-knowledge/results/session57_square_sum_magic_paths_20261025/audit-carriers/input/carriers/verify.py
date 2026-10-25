#!/usr/bin/env python3
"""Self-contained finite certificate. Python standard library only; no network."""
from itertools import combinations
from math import isqrt
import json
from pathlib import Path
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def require(condition, message='mathematical certificate check failed'):
 if not condition:raise RuntimeError(message)
def square(x):return x>=0 and isqrt(x)**2==x
def qr(p):return {x*x%p for x in range(p)}
def sums(m):return [sum(m[i] for i in l) for l in LINES]
def array(c,a,b):return [a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a]
def path_certificate(m,p=None):
 out=[]
 for l in LINES:
  s=0;path=[(0,0)]
  for k,i in enumerate(l,1):
   s+=m[i]
   if p:s%=p
   path.append((k,s))
  out.append({'cell_labels':l,'increments':[m[i] for i in l],'vertices':path})
 return out

def ham(n):
 adj={x:tuple(y for y in range(1,n+1) if y!=x and square(x+y)) for x in range(1,n+1)}
 count=0;first=None
 def dfs(path,seen):
  nonlocal count,first
  if len(path)==n:
   count+=1
   if first is None:first=path[:]
   return
  for y in adj[path[-1]]:
   bit=1<<(y-1)
   if not seen&bit:dfs(path+[y],seen|bit)
 for x in range(1,n+1):dfs([x],1<<(x-1))
 require(count % 2 == 0)
 if first:
  require(sorted(first) == list(range(1, n + 1)))
  require(all((square(x + y) for x, y in zip(first, first[1:]))))
 return {'n':n,'edges':sum(map(len,adj.values()))//2,'oriented_hamiltonian_paths':count,'unoriented_hamiltonian_paths':count//2,'first_path':first}

def main():
 out={}
 for name,p,roots in [('source_F29',29,[9,11,1,6,0,14,12,16,8]),('source_F59',59,[20,12,7,2,1,23,22,25,29])]:
  m=[x*x%p for x in roots];T=3*m[4]%p
  require(len(set(m)) == 9 and set(m) <= qr(p))
  require(all((s % p == T for s in sums(m))))
  opposite=[(m[i]+m[8-i])%p for i in range(4)]
  require(opposite == [2 * m[4] % p] * 4)
  out[name]={'p':p,'roots':roots,'residues':m,'total_mod_p':T,'integer_residue_line_sums':sums(m),'opposite_sum':opposite[0],'opposite_sum_is_quadratic_residue':opposite[0] in qr(p),'paths':path_certificate(m,p)}
 m=[4,29,15,27,16,5,17,3,28]
 require(len(set(m)) == 9 and set(m) <= qr(59) and (sums(m) == [48] * 8))
 require(m == array(16, 4, 15))
 out['weak_no_wrap_F59']={'residues':m,'integer_line_sums':sums(m),'actual_integer_square_cells':[i for i,x in enumerate(m) if square(x)],'roots_mod59':[next(r for r in range(59) if r*r%59==x) for x in m]}
 require(out['weak_no_wrap_F59']['actual_integer_square_cells'] == [0, 4])
 m=[8,1,6,3,5,7,4,9,2]
 require(sums(m) == [15] * 8)
 out['LoShu']={'entries':m,'integer_line_sums':sums(m),'row_adjacencies_in_square_sum_graph':[square(m[i]+m[j]) for l in LINES[:3] for i,j in zip(l,l[1:])]}
 out['square_sum_graphs']=[ham(n) for n in range(15,26)]
 q=qr(29);h=[t for t in combinations(sorted(q),3) if sum(t)%29==0]
 es={e for t in h for e in combinations(t,2)};false=(0,1,4)
 require(all((e in es for e in combinations(false, 2))) and sum(false) % 29 != 0)
 out['shadow_false_triangle']={'p':29,'total':0,'triangle':false,'edge_witnesses':[(e,next(t for t in h if set(e)<=set(t))) for e in combinations(false,2)]}
 pairs=[(i,j) for i,j in combinations(range(8),2) if set(LINES[i])&set(LINES[j])]
 require(len(pairs) == 22)
 require([sum((i in l for l in LINES)) for i in range(9)] == [3, 2, 3, 2, 4, 2, 3, 2, 3])
 out['incidence']={'lines':LINES,'line_intersection_edges':pairs,'line_intersection_density':'11/14','cell_line_incidence_edges':24,'incidence_cycle_rank':8}
 # Equal-total triples with total itself square need not have square entries.
 # Scale Lo Shu by 15: total 225=15^2; this is not a magic square OF squares.
 out['square_total_variant']={'entries':[15*x for x in m],'total':225}
 require(sums(out['square_total_variant']['entries']) == [225] * 8 and square(225))
 require(not all(map(square, out['square_total_variant']['entries'])))
 return out
if __name__=='__main__':
 result=main()
 target=Path(__file__).with_name('certificate.json')
 data=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if '--check' in __import__('sys').argv:
  require(target.read_text() == data, 'Stored certificate differs from recomputation')
  print('PASS: all mathematical assertions and exact stored certificate verified')
 else:
  target.write_text(data)
  print('PASS: wrote',target)
 print('Hamiltonian counts:',[(x['n'],x['unoriented_hamiltonian_paths']) for x in result['square_sum_graphs']])
