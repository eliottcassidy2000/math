"""Bounded SAT scout for necessary endpoint and shared-cubic congruences, not full equality."""
import argparse,json,time
from pathlib import Path
from hashlib import sha256
from pysat.solvers import Solver

def main():
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=12);p.add_argument('--budget',type=int,default=20000);p.add_argument('--output',type=Path,required=True)
 args=p.parse_args();n=args.n
 if n<4 or n%2:raise RuntimeError('positive even n>=4')
 source=Path(__file__).read_bytes()
 solver=Solver(name='glucose4')
 nxt=1;clauses=0
 def fresh():
  nonlocal nxt
  nxt+=1;return nxt
 def add(c):
  nonlocal clauses
  solver.add_clause(c);clauses+=1
 add([1]) # true literal
 def land(a,b):
  if a==1:return b
  if b==1:return a
  if a==-1 or b==-1 or a==-b:return -1
  if a==b:return a
  z=fresh();add([-z,a]);add([-z,b]);add([z,-a,-b]);return z
 def xor2(a,b):
  if a==-1:return b
  if b==-1:return a
  if a==1:return -b
  if b==1:return -a
  if a==b:return -1
  if a==-b:return 1
  z=fresh()
  for c in ([-a,-b,-z],[a,b,-z],[a,-b,z],[-a,b,z]):add(c)
  return z
 def xor(xs):
  level=list(xs)
  if not level:return -1
  while len(level)>1:
   level=[xor2(level[i],level[i+1]) if i+1<len(level) else level[i] for i in range(0,len(level),2)]
  return level[0]
 a=[[0]*n for _ in range(n)]
 for u in range(n):
  for v in range(u+1,n):
   z=fresh();a[u][v]=z;a[v][u]=-z
 # Every tournament has a Hamilton path, so a labeling with these arcs exists.
 for u in range(n-1):add([a[u][u+1]])
 full=(1<<n)-1;s=[{} for _ in range(full)]
 for mask in range(1,full):
  vs=[u for u in range(n) if mask>>u&1]
  for u in vs:
   rest=mask^(1<<u)
   s[mask][u]=1 if not rest else xor(land(a[u][v],s[rest][v]) for v in vs if v!=u)
  # Inherited Redei parity is redundant mathematically but aids propagation.
  add([xor(s[mask].values())])
 for u in range(n):
  add([-xor(s[m][u] for m in range(1,full) if m>>u&1)])
 for v in range(n):
  for u in range(n):
   if u==v:continue
   add([-xor(s[m][u] for m in range(1,full) if m>>u&1 and not m>>v&1
             and m.bit_count()%2==0 and m.bit_count()<=n-4)])
 solver.conf_budget(args.budget)
 start=time.perf_counter();answer=solver.solve_limited()
 elapsed=time.perf_counter()-start
 rows=None
 if answer:
  model=set(solver.get_model());rows=[sum(1<<v for v in range(n) if v!=u and a[u][v] in model) for u in range(n)]
 result=dict(status='SAT necessary-filter candidate' if answer else 'UNSAT solver report, not independently certified' if answer is False else 'UNKNOWN: conflict budget exhausted',
             n=n,conflict_budget=args.budget,variables=nxt,clauses=clauses,solver_stats=solver.accum_stats(),
             rows=rows,elapsed_seconds=round(elapsed,3),
             scope='Hamilton-path labeling; all endpoint parities balanced and all shared cubic residues zero. No C-skew condition or full-response equality imposed.',
             source_sha256=sha256(source).hexdigest())
 solver.delete();args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps(result))
if __name__=='__main__':main()

