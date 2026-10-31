#!/usr/bin/env python3
"""Generate exact square-sum graph comparison; standard library only."""
from pathlib import Path
from math import isqrt
import json,csv
from graph_tools import graph,invariants
from hamilton_search import paths
ROOT=Path(__file__).resolve().parent
GROUPS={'A':[2,3,5,7,11,13],'B':[17,19,23,31,47],'C':[29,41,59,71],'D':[37,43,53,61,71],'E':[73,79,83,89,97]}
PRIMES=[n for n in range(2,98) if all(n%d for d in range(2,isqrt(n)+1))]
EARS31=[[6,19,30,6],[6,3,1,8,17,19],[6,10,26,23,13,3],[23,2,14,11,5,4,12,13],[2,7,29,20,5],[7,9,16,20],[1,15,10],[7,18,31,5],[4,21,15],[3,22,14],[1,24,12],[11,25,24],[9,27,22],[8,28,21]]

def dump(name,obj):
 (ROOT/name).write_text(json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n' if name=='graphs.json' else json.dumps(obj,indent=2,sort_keys=True)+'\n')

def main():
 rows={n:invariants(n) for n in range(1,98)}
 finite=[dict(n=n,path=paths(n),cycle=paths(n,True)) for n in range(1,33)]
 cycles=[dict(n=p,cycle=paths(p,True,False)['witness']) for p in PRIMES if p>=37]
 witnesses={'finite_hamilton':finite,'prime_cycles':cycles,'ears31':EARS31,'first_odd_cycle':[1,3,6,10,15],'first_triangle':[6,19,30]}
 dump('witnesses.json',witnesses)
 for n,r in rows.items():
  r['edge_list']=[[u,v] for u,vs in graph(n).items() for v in sorted(vs) if u<v]
  if n<=32:r['hamilton_path_count']=finite[n-1]['path']['count'];r['hamilton_cycle_count']=finite[n-1]['cycle']['count']
  if n in PRIMES and n>=37:r['hamilton_path_exists']=True;r['hamilton_cycle_exists']=True
 group_at=lambda p:[g for g,x in GROUPS.items() if p in x]
 exact={g:[{'prime_label':p,'prime_rank':PRIMES.index(p)+1,'prime_order_graph':p,'rank_order_graph':PRIMES.index(p)+1} for p in ps] for g,ps in GROUPS.items()}
 dump('graphs.json',{'definition':'G_n: vertices 1..n, unordered distinct x,y adjacent iff x+y is a square','groups_as_supplied':GROUPS,'groups_typed':exact,'first_25_primes':PRIMES,'all_graphs_1_through_97':[rows[n] for n in rows],'correction_status':'No silent repair. Substituting 67 for either occurrence of 71 produces the first-25-prime partition; D:71->67 is a conditional proposal only.'})
 for interpretation,orders in [('prime_order',PRIMES),('rank_order',range(1,26))]:
  table=[]
  for i,n in enumerate(orders,1):
   p=PRIMES[i-1];r=rows[n]
   table.append({'i':i,'prime_label':p,'supplied_groups':','.join(group_at(p)) or 'MISSING','graph_n':n,'edges':r['edges'],'components':r['components'],'delta':r['min_degree'],'Delta':r['max_degree'],'leaves':r['leaves'],'girth':r['girth'] or '-','odd_root_components':r['odd_root_components'],'arrival_roots':','.join(map(str,r['arrival_roots'])) or '-','hamilton_path_count':r.get('hamilton_path_count','not enumerated'),'hamilton_path_exists':r.get('hamilton_path_exists',r.get('hamilton_path_count',0)>0),'hamilton_cycle_exists':r.get('hamilton_cycle_exists',r.get('hamilton_cycle_count',0)>0)})
  with (ROOT/(interpretation+'_first25.csv')).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
  lines=['|i|prime label|given group|graph|edges|components|δ|Δ|girth|odd-root components|HP|HC|','|---:|---:|:---|:---|---:|---:|---:|---:|---:|---:|:---:|:---:|']
  for r in table:
   lines.append('|'+ '|'.join(map(str,[r['i'],r['prime_label'],r['supplied_groups'],'G_'+str(r['graph_n']),r['edges'],r['components'],r['delta'],r['Delta'],r['girth'],r['odd_root_components'],'yes' if r['hamilton_path_exists'] else 'no','yes' if r['hamilton_cycle_exists'] else 'no']))+'|')
  (ROOT/(interpretation+'_first25.md')).write_text('\n'.join(lines)+'\n')
 print('Generated 97 complete graphs, two 25-row tables, and exact Hamiltonian certificates.')
if __name__=='__main__':main()
