"""Independent bit-parity, closure, literal-query, and integral-root-frame controls."""
import argparse,json,sys
from pathlib import Path
from hashlib import sha256
from itertools import permutations
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'homogeneous_pair_response_20261103'))
import check
need=check.need

def parities(a):
 n=len(a);full=(1<<n)-1;s=[0]*(full+1);e=[0]*(full+1)
 incoming=[sum(1<<u for u in range(n) if a[u]>>v&1) for v in range(n)]
 for m in range(1,full+1):
  todo=m
  while todo:
   bit=todo&-todo;todo-=bit;u=bit.bit_length()-1;rest=m^bit
   if not rest or (a[u]&s[rest]).bit_count()%2:s[m]|=bit
   if not rest or (incoming[u]&e[rest]).bit_count()%2:e[m]|=bit
 return s[-1],e[-1]

def prime_closure(a):
 n=len(a);full=(1<<n)-1
 for u in range(n):
  for v in range(u+1,n):
   mask=(1<<u)|(1<<v)
   while True:
    add=0
    for w in range(n):
     if not mask>>w&1 and (a[w]&mask) not in (0,mask):add|=1<<w
    if not add:break
    mask|=add
   if mask!=full:return False
 return True

def skew(a):
 return [[0 if u==v else 1 if a[u]>>v&1 else -1 for v in range(len(a))] for u in range(len(a))]

def curl(a):
 s=skew(a);r=len(a)-1
 return [[s[i][j]-s[i][r]+s[j][r] for j in range(r)] for i in range(r)]

def realize(k):
 r=len(k)
 if any(len(row)!=r for row in k):return None
 if any(k[i][i] for i in range(r)):return None
 if any(k[i][j]!=-k[j][i] or abs(k[i][j]) not in (1,3) for i in range(r) for j in range(i+1,r)):return None
 adj=[[] for _ in range(2*r)];rev=[[] for _ in adj]
 for i in range(r):
  for j in range(i+1,r):
   for b in (0,1):
    for c in (0,1):
     if abs(k[i][j]+(2*b-1)-(2*c-1))!=1:
      for u,v in [(2*i+b,2*j+(1-c)),(2*j+c,2*i+(1-b))]:
       adj[u].append(v);rev[v].append(u)
 seen=set();order=[]
 def visit(u):
  seen.add(u)
  for v in adj[u]:
   if v not in seen:visit(v)
  order.append(u)
 for u in range(2*r):
  if u not in seen:visit(u)
 comp=[-1]*(2*r)
 def assign(u,c):
  comp[u]=c
  for v in rev[u]:
   if comp[v]<0:assign(v,c)
 number=0
 for u in reversed(order):
  if comp[u]<0:assign(u,number);number+=1
 if any(comp[2*i]==comp[2*i+1] for i in range(r)):return None
 signs=[1 if comp[2*i+1]>comp[2*i] else -1 for i in range(r)]
 a=[0]*(r+1)
 for i,t in enumerate(signs):
  if t==1:a[i]|=1<<r
  else:a[r]|=1<<i
 for i in range(r):
  for j in range(i+1,r):
   t=k[i][j]+signs[i]-signs[j]
   need(t in (-1,1),'2SAT assignment')
   if t==1:a[i]|=1<<j
   else:a[j]|=1<<i
 need(curl(a)==k,'curl inverse')
 return a

def one(a,mask):
 n=len(a);b=a[:]+[0]
 for u in range(n):
  if mask>>u&1:b[u]|=1<<n
  else:b[n]|=1<<u
 return b

def triple_count(s):
 total=0;n=len(s)
 for a,b,c,d,e in permutations(range(n),5):
  wac=s[a][d]-s[a][e]-s[b][d]+s[b][e]
  if wac!=-4:continue
  total+=(s[a][b]-s[a][c]+s[b][c],s[b][d]-s[b][e]-s[c][d]+s[c][e])==(-3,2)
  total+=(s[a][c]-s[b][c]+s[b][a],s[c][d]-s[c][e]-s[a][d]+s[a][e])==(-3,2)
 return total

def root_triple(s):
 roots=list(permutations(range(len(s)),2))
 g=[[(i==k)-(i==l)-(j==k)+(j==l) for k,l in roots] for i,j in roots]
 w=[[s[i][k]-s[i][l]-s[j][k]+s[j][l] for k,l in roots] for i,j in roots]
 total=0
 for a in range(len(roots)):
  for b in range(len(roots)):
   if g[a][b]!=-1 or w[a][b]!=-3:continue
   total+=sum(g[a][c]==g[b][c]==0 and w[a][c]==-4 and w[b][c]==2 for c in range(len(roots)))
 return total

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'audit.json');args=p.parse_args()
 source=Path(__file__).read_bytes();data=json.loads((HERE/'switch.json').read_text())
 prime_count=0;controls=[]
 for rec in data['records']:
  m=rec['mask'];a=[x^((m^(4095 if m>>u&1 else 0))&~(1<<u)) for u,x in enumerate(data['base_rows'])]
  need(parities(a)==(rec['S_parity'],rec['E_parity']),'independent parity orbit')
  isprime=prime_closure(a);need(isprime==rec['prime'],'closure vs subset primality')
  prime_count+=isprime
  if m in (0,1,3,11):
   need(check.h(a)==rec['H'],'independent integer H')
   s=skew(a);need(triple_count(s)==root_triple(s)==660,'two independent root-triple paths')
   state=[1]*12;walks=[]
   for length in range(17):
    count=sum(state);walks.append(count)
    need(count==(12 if not length else 66*5**(length-1)),'walk count')
    state=[sum(state[v] for v in range(12) if a[u]>>v&1) for u in range(12)]
   rebuilt=realize(curl(a));need(rebuilt is not None and check.h(rebuilt)==rec['H'],'lattice H')
   if isprime:need(rebuilt==a,'strong labelled curl inverse')
   controls.append(dict(mask=m,rows=a,H=rec['H'],selected_genus3_coefficient=660,walk_counts=walks,
                        simple_paths_by_vertices=rec['simple_paths_by_vertices']))
 # All small labelled inputs: inverse lattice feasibility and H preservation.
 lattice_count=0;literal_current_checks=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   a=check.tournament(n,code);b=realize(curl(a))
   need(b is not None and check.literal_h(a)==check.literal_h(b),'all small lattice fibers')
   lattice_count+=1
   h,s,e,c=check.sidecar(a)
   for u in (range(n) if n<=4 else [0]):
    expected=2*(e[u]-s[u])+sum(c[u])-sum(row[u] for row in c)
    actual=check.literal_h(one(a,1<<u))-check.literal_h(one(check.opposite(a),1<<u))
    need(actual==expected,'literal current formula');literal_current_checks+=1
 impossible=[[0,3,3],[-3,0,3],[-3,-3,0]]
 need(realize(impossible) is None,'infeasible root form hostile')
 # Score-parity regularity alone does not explain the conference result.
 hostile=[110,148,248,194,137,26,50,97]
 sp,ep=parities(hostile)
 need(all(x.bit_count()%2==1 for x in hostile) and sp^ep==197,'score parity hostile')
 need(Path(__file__).read_bytes()==source,'source changed')
 result=dict(status='PASS independent algorithms; no universal prime12 claim',source_sha256=sha256(source).hexdigest(),
             switch_receipt_sha256=sha256((HERE/'switch.json').read_bytes()).hexdigest(),
             full_orbit_parity_checks=2048,full_orbit_closure_checks=2048,prime_count=prime_count,
             controls=controls,small_labelled_lattice_inputs=lattice_count,
             literal_one_port_current_checks=literal_current_checks,infeasible_form=impossible,
             score_parity_hostile=dict(rows=hostile,boundary_mask=197))
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('controls','source_sha256')},indent=2))
if __name__=='__main__':main()

