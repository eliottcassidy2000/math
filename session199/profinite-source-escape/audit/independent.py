from collections import deque
from pathlib import Path
import hashlib, json, importlib.util

# This loads the independent Walsh checker, not the target packet implementation.
spec=importlib.util.spec_from_file_location('audit',Path(__file__).parent/'verify_walsh.py')
a=importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def modular_bfs(M):
 start=(1%M,0,2%M); parents={start:None}; queue=deque([start])
 while queue:
  D,E,U=queue.popleft()
  for label,t in [('O',(D,0,2*U%M)),('Z',((D-E+U)%M,1,(2*U-1)%M))]:
   if t not in parents:parents[t]=((D,E,U),label);queue.append(t)
 return parents

def recipe(parents,target):
 word=[]
 while parents[target] is not None:
  target,letter=parents[target];word.append(letter)
 return ''.join(reversed(word))

pairs=0;ghosts=[]
for M in range(2,102,2):
 parents=modular_bfs(M)
 for D in range(1,M,2):
  for U in range(1,M,2):
   need((D,1,U) in parents, 'independent modular BFS missing unit pair')
   pairs+=1
 m=M//2;word=recipe(parents,(1,1,3%M))+'O';y=a.state(word);p=a.old(y)
 need(all((x-z)%m==0 for x,z in zip(p,(1,0,6,1))),'phantom witness')
 if m%24==0:
  need(p[0]>=m+1 and p[2]>=2*m+2 and 5*p[0]+11*p[2]-4>=27*m+23,'height bound')
 ghosts.append({'m':m,'word_length':len(word),'D':p[0],'U':p[2]})
# Independent positive-coordinate source bank and independent DP for all sources.
bank={a.b:''}; todo=[a.b]
for y in todo:
 for letter,A in a.GEN.items():
  z=a.mul(A,y)
  if a.old(z)[2]<=11 and z not in bank:bank[z]=bank[y]+letter;todo.append(z)
need(len(bank)==14,'cap11 cardinality')
for y,word in bank.items():
 p,s=a.direct_profile(*a.build(word));need(p==a.old(y) and s==1,'direct physical bank')
support={15}
for x in bank:
 for y in bank:
  for P in [a.PN,a.PM]:
   h=a.dot(x,a.mul(P,y))
   if h<=67:support.add(h)
expected={13,15,23,25,29,41,43,45,49,53,55}
need(support==expected and 67 not in support,'exact complete C2 support')
b=a.state('Z');p=(0,4,0,1)
need(a.dot(b,a.mul(a.PN,p))==67,'phantom value')
report={'status':'PASS','proof_sha256':hashlib.sha256((Path(__file__).parents[1]/'PROOF.md').read_bytes()).hexdigest(),'independent_modular_BFS_moduli':'every even M from2 through100','unit_interface_pairs':pairs,'phantom_profiles_reconstructed':len(ghosts),'all_14_cap11_sources_direct_graph_checked':True,'C2_ordered_pair_comparisons':14*14*2,'C2_support_through67':sorted(support),'manual_all_order_audit':'Compiler, all-moduli sequence, real-height escape, sparse-family construction, and bounded-height equivalence all valid. No correction required.','scope':'No positive-image coverage inference.'}
(Path(__file__).parent/'INDEPENDENT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
