"""Independent permutation path enumeration and literal named-response reconstruction."""
from pathlib import Path
from itertools import permutations
import json,hashlib,tempfile,subprocess
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'marked_path_duality'
def need(c,m):
 if not c:raise RuntimeError(m)
def add(*Ps):
 D={}
 for P in Ps:
  for k,v in P.items():D[k]=D.get(k,0)+v
 return {k:v for k,v in D.items() if v}
def scale(k,P):return {m:k*v for m,v in P.items() if k*v}
def mul(A,B):
 D={}
 for m,v in A.items():
  for s,w in B.items():
   if not m&s:D[m|s]=D.get(m|s,0)+v*w
 return {k:v for k,v in D.items() if v}
def rev(P):return {m:(-1)**m.bit_count()*v for m,v in P.items()}
def paths(T):
 n=len(T);G=[[{} for _ in T] for _ in T]
 for k in range(1,n+1):
  for w in permutations(range(n),k):
   if all(T[u]>>v&1 for u,v in zip(w,w[1:])):
    P=G[w[0]][w[-1]];m=sum(1<<v for v in w);P[m]=P.get(m,0)+1
 return G

def algebra(T):
 n=len(T);full=(1<<n)-1;G=paths(T);C=[add(*(G[u][v] for u in range(n))) for v in range(n)];R=[add(*r) for r in G];H=add(*C);U=add({0:1},H);HH=mul(H,H);need(mul(U,rev(U))=={0:1},'every integer reciprocity coefficient');out=[];checked=1
 for p in range(n):need(not add(C[p],mul(U,rev(R[p]))),'every column-row coefficient');checked+=1
 for p in range(n):
  row=[]
  for q in range(n):
   E=add(G[p][q],scale(-1,G[q][p]));O={m:v for m,v in E.items() if m.bit_count()%2};D=add(mul(C[p],R[q]),scale(-1,mul(C[q],R[p])));need(not add(D,scale(2,mul(U,O))),'every marked integer identity coefficient');checked+=1
   K=add(D,scale(-1,mul(H,D)),scale(2,mul(H,E)),mul(HH,E));L=add(scale(2,mul(add(HH,{0:-1}),O)),mul(add(scale(2,H),HH),E));need(K==L,'complete kernel polynomials, not just top support');row.append(K.get(full,0))
  out.append(row)
 return out,checked

def aug(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  for bits,w in [(a,n),(c,n+1)]:
   if bits>>v&1:R[v]|=1<<w
   else:R[w]|=1<<v
 return tuple(R)
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def cyc(T):return [r|(1<<len(T)) for r in T]+[1<<(len(T)+1),(1<<len(T))-1]
def subs(s):
 v=s
 while True:
  yield v
  if not v:break
  v=(v-1)&s

def main():
 d=json.loads((SRC/'normal.json').read_text());algebra_records=[];identities=0
 for r in d['complete_small_algebra_controls']+d['larger_controls']:
  D,k=algebra(r['rows']);need(D==r['quartic_matrix'],'every independent permutation-derived marked matrix');identities+=k;algebra_records.append([r['rows'],D])
 expected={}
 def take(T,U,qs):
  for q in qs:
   for side,W in enumerate([T,U]):
    G=aug(W,q['a'],q['c']);v=q['values'][side];need(G not in expected or expected[G]==v,'consistent repeated graph');expected[G]=v
 for r in d['literal_quartic_controls']:
  T=cyc(r['child_rows']);U=cyc(op(r['child_rows']));take(T,U,r['queries']);v=sum(q['sign']*(q['values'][0]-q['values'][1]) for q in r['queries']);need(v==r['expected'],'all physical fourth differences')
 b=d['Paley_complete_named_bank'];T,U=b['source'],b['target'];take(T,U,b['coordinates']);take(T,U,b['full_mask_controls']);lookup={(q['a'],q['c']):q['values'] for q in b['coordinates']};hist={k:0 for k in range(5)}
 for q in b['coordinates']:
  a,c=q['a'],q['c'];v=[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*lookup[x,y][side] for x in subs(a) for y in subs(c)) for side in [0,1]];need(v==q['coefficients'],'every complete named coefficient');hist[a.bit_count()+c.bit_count()]+=v[0]!=v[1]
 for q in b['full_mask_controls']:
  need(q['values']==[sum(r['coefficients'][side] for r in b['coordinates'] if r['a']&q['a']==r['a'] and r['c']&q['c']==r['c']) for side in [0,1]],'all full-mask polynomial reconstructions')
 need({str(k):v for k,v in hist.items()}==b['nonzero_difference_counts'],'full-degree histogram')
 graphs=list(expected);need(len(graphs)==d['distinct_actual_graphs'],'complete physical graph set');wire='\n'.join(str(len(G))+'\n'+' '.join(map(str,G)) for G in graphs)+'\n'
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'p';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'physical.cpp'),'-o',str(exe)],check=True);values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()))
 need(values==[expected[g] for g in graphs],'all actual counts by independent predecessor DP');P=algebra_records[-1][1];need({P[p][q] for p in range(7) for q in range(7) if p!=q}=={-8,8},'Paley nonzero integer mod8-blind boundary');w=b['actual_residue_two_witness'];need((w['values'][0]-w['values'][1])%4==2,'full response remains visible')
 pins={f:hashlib.sha256((SRC/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','paths.py','count.cpp','normal.json']}
 print(json.dumps({'status':'PASS','source_pins':pins,'audit_pins':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['audit.py','physical.cpp']},'complete_algebra_controls':len(algebra_records),'coefficientwise_identity_count':identities,'algebra_sha256':hashlib.sha256(repr(algebra_records).encode()).hexdigest(),'literal_quartics':len(d['literal_quartic_controls']),'distinct_actual_graphs':len(graphs),'physical_values_sha256':hashlib.sha256(repr(values).encode()).hexdigest(),'named_coordinates':len(b['coordinates']),'nonzero_difference_histogram':hist,'Paley_quartic_matrix':P,'actual_residue_two_witness':w},sort_keys=True,indent=2))
if __name__=='__main__':main()
