"""Exact integer marked-path reciprocity and complete physical controls."""
from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,tempfile
from paths import matrix,plus,minus,product,full
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent

def need(x,s):
 if not x:raise RuntimeError(s)
def rows(n,code):
 T=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else(j,i);T[u]|=1<<v
 return T
def tt(n):return [sum(1<<v for v in range(u+1,n)) for u in range(n)]
def op(T):return [((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T)]
def cyclic(A):
 n=len(A);return [r|1<<n for r in A]+[1<<(n+1),(1<<n)-1]
def augment(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for v in range(n):
  for w,bits in((n,a),(n+1,c)):
   if bits>>v&1:U[v]|=1<<w
   else:U[w]|=1<<v
 return tuple(U)
def identities(A):
 n=len(A);N=1<<n;G=matrix(A);C=[plus(*(G[u][v] for u in range(n))) for v in range(n)];R=[plus(*G[u]) for u in range(n)];H=plus(*C);U=list(H);U[0]=1;rev=lambda p:[z*(-1)**s.bit_count() for s,z in enumerate(p)];one=[1]+[0]*(N-1);need(product(U,rev(U))==one,'coefficientwise integer reciprocity')
 for p in range(n):need(plus(C[p],product(U,rev(R[p])))==[0]*N,'integer column-row duality')
 HH=product(H,H);out=[]
 for p in range(n):
  row=[]
  for q in range(n):
   E=minus(G[p][q],G[q][p]);EO=[v if s.bit_count()%2 else 0 for s,v in enumerate(E)];D=minus(product(C[p],R[q]),product(C[q],R[p]));need(D==[-2*v for v in product(U,EO)],'every integer marked coefficient');delta=D[-1]-full(H,D)+2*full(H,E)+full(HH,E);alt=2*(full(HH,EO)-EO[-1])+2*full(H,E)+full(HH,E);need(delta==alt,'exact quartic specialization');row.append(delta)
  out.append(row)
 return out

def main():
 controls=[];small=[]
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   A=rows(n,code);D=identities(A);controls.append(dict(rows=A,edge_code=code,quartic_matrix=D))
   if n<=4:small.append((A,D))
 P7=[sum(1<<((i+d)%7) for d in(1,2,4)) for i in range(7)];larger=[]
 for name,A in [('TT6',tt(6)),('TT7',tt(7)),('cyclic7',[sum(1<<((i+d)%7) for d in(1,2,3)) for i in range(7)]),('Paley7',P7)]:larger.append(dict(name=name,rows=A,quartic_matrix=identities(A)))
 pmat=larger[-1]['quartic_matrix'];need({pmat[p][q] for p in range(7) for q in range(7) if p!=q}=={-8,8},'Paley exact mod8 boundary')
 need(all(all(v==0 for v in r['quartic_matrix'][p]) for r in controls if len(r['rows'])==2 for p in range(2)),'pair quartic zero boundary')
 graphs=[];lookup={}
 def val(T,a,c):
  G=augment(T,a,c)
  if G not in lookup:lookup[G]=len(graphs);graphs.append(G)
  return lookup[G]
 literal=[]
 for A,D in small:
  n=len(A);T=cyclic(A);U=cyclic(op(A))
  for p in range(n):
   for q in range(n):
    if p==q:continue
    qs=[]
    for k in range(16):a=(1<<p if k&1 else 0)|(1<<n if k&2 else 0);c=(1<<q if k&4 else 0)|(1<<(n+1) if k&8 else 0);qs.append(dict(a=a,c=c,sign=(-1)**k.bit_count(),ids=[val(V,a,c) for V in(T,U)]))
    literal.append(dict(child_rows=A,p=p,q=q,expected=D[p][q],queries=qs))
 T=cyclic(P7);U=cyclic(op(P7));n=len(T);masks=[s for s in range(1<<n) if s.bit_count()<=2];bank=[dict(a=a,c=c,ids=[val(V,a,c) for V in(T,U)]) for a in masks for c in masks];extra=[dict(a=(41*k+13)&511,c=(67*k+29)&511,ids=[val(V,(41*k+13)&511,(67*k+29)&511) for V in(T,U)]) for k in range(32)]
 with tempfile.TemporaryDirectory(prefix='marked-path-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'all exact actual graphs')
 for r in literal:
  for q in r['queries']:q['values']=[values[i] for i in q.pop('ids')]
  need(sum(q['sign']*(q['values'][0]-q['values'][1]) for q in r['queries'])==r['expected'],'every literal small-child quartic')
 table={}
 for r in bank:r['values']=[values[i] for i in r.pop('ids')];table[r['a'],r['c']]=r['values']
 def subs(s):
  a=s
  while True:
   yield a
   if not a:break
   a=(a-1)&s
 for r in bank:r['coefficients']=[sum((-1)**(r['a'].bit_count()+r['c'].bit_count()-a.bit_count()-c.bit_count())*table[a,c][j] for a in subs(r['a']) for c in subs(r['c'])) for j in range(2)]
 for r in extra:r['values']=[values[i] for i in r.pop('ids')];need(r['values']==[sum(z['coefficients'][j] for z in bank if not z['a']&~r['a'] and not z['c']&~r['c']) for j in range(2)],'complete named polynomial reconstruction')
 for p in range(7):
  for q in range(7):
   if p==q:continue
   z=next(r for r in bank if r['a']==(1<<p|128) and r['c']==(1<<q|256));need(z['coefficients'][0]-z['coefficients'][1]==pmat[p][q],'all Paley marked quartics')
 witness=next((r for r in bank if (r['values'][0]-r['values'][1])%4==2),None);need(witness is not None,'kernel blindness is not full-response collision');counts={d:sum(z['coefficients'][0]!=z['coefficients'][1] for z in bank if z['a'].bit_count()+z['c'].bit_count()==d) for d in range(5)}
 paths=[ROOT/f for f in('PROOF.md','verify.py','paths.py','count.cpp')]+[BASE/'cyclic_child_disjoint'/'PROOF.md',BASE/'response'/'PROOF.md']
 return dict(status='PASS',complete_small_algebra_controls=controls,larger_controls=larger,literal_quartic_controls=literal,Paley_complete_named_bank=dict(source=T,target=U,coordinates=bank,full_mask_controls=extra,nonzero_difference_counts=counts,actual_residue_two_witness=witness),distinct_actual_graphs=len(graphs),source_sha256={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
