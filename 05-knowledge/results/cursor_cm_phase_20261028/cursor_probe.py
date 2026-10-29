"""Search cursor faithfulness using exact endpoint automorphisms and orbit refinements."""
from pathlib import Path
import importlib.util,json,argparse
from itertools import permutations,product
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('ports',ROOT/'05-knowledge/results/joint_port_germs_20261028/port_core.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
o=c.oracle
def endpoint_group(t):
 return [g for g in o.isomorphisms(t,t,True) if g[0]==0 and g[-1]==len(t)-1]
def orbits(G,n):
 remaining=set(range(n));ans=[]
 while remaining:
  v=min(remaining);b=tuple(sorted({g[v] for g in G}));ans.append(b);remaining.difference_update(b)
 return ans
def partition_maps(blocks,n):
 for parts in product(*(tuple(permutations(b)) for b in blocks)):
  p=list(range(n))
  for a,b in zip(blocks,parts):
   for x,y in zip(a,b):p[x]=y
  yield tuple(p)
def scan_pair(t,u):
 if len(t)!=len(u):return False
 n=len(t)
 for a,b in ((t,u),(c.tau(t),c.tau(u))):
  maps=[g for g in o.isomorphisms(a,b,True) if g[0]==0 and g[-1]==n-1]
  if not maps or not all(any(g[i]==i for g in maps) for i in range(n)):return False
 return True
def analyze(t):
 n=len(t);G=endpoint_group(t)
 if len(G)==1:return None
 u=c.tau(t);H=endpoint_group(u)
 if len(H)==1:return None
 A=orbits(G,n);B=orbits(H,n)
 intersections=[sorted(set(a)&set(b)) for a in A for b in B if len(set(a)&set(b))>1]
 candidates=0;hostile=None
 for perm in partition_maps(A,n):
  v=c.relabel(t,perm)
  if v==t or not all(c.edge(v,i,i+1) for i in range(n-1)):continue
  candidates+=1
  if scan_pair(t,v):
   hostile={'left':t,'right':v,'permutation':perm};break
 return {'rows':t,'groups':[len(G),len(H)],'orbits':[A,B],'nontrivial_intersections':intersections,
         'candidate_paths':candidates,'collision':hostile}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--max-order',type=int,default=7)
 args=parser.parse_args()
 rows=[];examples=[];collisions=[]
 for n in range(1,args.max_order+1):
  count=1<<((n-1)*(n-2)//2);dual=0;overlap=0;trials=0
  for mask in range(count):
   result=analyze(c.presentation(n,mask))
   if result is None:continue
   dual+=1;overlap+=bool(result['nontrivial_intersections']);trials+=result['candidate_paths']
   if len(examples)<12:examples.append(result)
   if result['collision']:collisions.append(result)
  row={'n':n,'states':count,'both_endpoint_groups_nontrivial':dual,'overlapping_orbits':overlap,'candidate_paths':trials}
  rows.append(row);print(json.dumps(row),flush=True)
 out={'status':'FINITE-EXACT','rows':rows,'examples':examples,'collisions':collisions}
 (HERE/'cursor_probe.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'collision_count':len(collisions)}),flush=True)
if __name__=='__main__':main()
