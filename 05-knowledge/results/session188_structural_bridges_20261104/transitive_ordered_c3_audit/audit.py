"""Independent complete exterior-word kernel and literal count audit."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import hashlib,json,math,subprocess,tempfile,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'transitive_ordered_c3'
def need(ok,label):
 if not ok:raise RuntimeError(label)
def add(*ps):
 out={}
 for P in ps:
  for e,c in P.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items() if c}
def scale(P,c):return {e:v*c for e,v in P.items() if v*c}
CAP=(3,3,8,8)
def mul(*ps):
 out={(0,0,0,0):1}
 for P in ps:
  nxt={}
  for a,x in out.items():
   for b,y in P.items():
    e=tuple(u+v for u,v in zip(a,b))
    if all(z<=c for z,c in zip(e,CAP)):nxt[e]=nxt.get(e,0)+x*y
  out={e:c for e,c in nxt.items() if c}
 return out

def exterior_kernel():
 # Direct paths in X,Y,F,L, rather than the producer eliminated-interface formula.
 R=(6,4,8,3);G=[[{} for _ in R] for _ in R]
 @lru_cache(None)
 def f(counts,last):
  if not counts[last]:return (0,)*4
  cs=list(counts);cs[last]-=1;cs=tuple(cs)
  if not any(cs):return tuple(int(i==last) for i in range(4))
  out=[0]*4
  for p in range(4):
   if R[p]>>last&1:
    for first,z in enumerate(f(cs,p)):out[first]+=z
  return tuple(out)
 for counts in product(*(range(c+1) for c in CAP)):
  if not any(counts):continue
  for last in range(4):
   for first,z in enumerate(f(counts,last)):
    if z:G[first][last][counts]=z
 E=[add(*(G[u][v] for u in range(4))) for v in range(4)];S=[add(*G[u]) for u in range(4)];one={(0,0,0,0):1};I=(0,2,3);alpha=add(one,*(E[v] for v in I));beta=add(one,S[1]);tau=add(*(G[1][v] for v in I));rho=[G[1][v] for v in(0,1)];ell=[add(*(G[u][v] for v in I)) for u in(0,1)]
 kernel=add(scale(mul(alpha,beta,add(G[0][1],scale(G[1][0],-1))),2),scale(mul(add(one,scale(tau,2)),add(mul(E[0],S[1]),scale(mul(E[1],S[0]),-1))),-1),mul(alpha,add(mul(rho[0],S[1]),scale(mul(rho[1],S[0]),-1))),mul(beta,add(mul(E[0],ell[1]),scale(mul(E[1],ell[0]),-1))))
 weights=[0,1,2,2];out={}
 for(s,t,a,b),c in kernel.items():out[a,b]=out.get((a,b),0)+c*weights[s]*weights[t]
 return {e:c for e,c in out.items() if c},len(kernel)
def q(n,k):return sum((-1)**(k-j)*math.comb(k,j)*j**n for j in range(k+1))
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def extend(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for u in range(n):
  U[u if a>>u&1 else n]|=1<<(n if a>>u&1 else u)
  U[u if c>>u&1 else n+1]|=1<<(n+1 if c>>u&1 else u)
 return tuple(U)
def collapsed(T,a,b,am,cm):
 U=extend(T,am,cm);n=len(T);sets=[((1<<a)-1)<<9,((1<<b)-1)<<(9+a)];need(all(module(U,M) for M in sets),'actual queried chains remain homogeneous');names=list(range(9))+[n,n+1,9,9+a]
 return tuple(sum(1<<j for j,v in enumerate(names) if U[u]>>v&1) for u in names)
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,D in enumerate(blocks):
  for u,row in enumerate(D):T[off[i]+u]=(row<<off[i])|sum(((1<<len(B))-1)<<off[j] for j,B in enumerate(blocks) if Q[i]>>j&1)
 return tuple(T)
def tt(n):return tuple(((1<<n)-1)^((1<<(u+1))-1) for u in range(n))
def main():
 paths=[SRC/x for x in('PROOF.md','verify.py','algebra.py','run_count.cpp','literal.cpp')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());poly,terms=exterior_kernel();expected={(s,t):c for s,t,c in data['coefficient_polynomial']['terms']};need(poly==expected and len(poly)==23,'all23 coefficients from independent full exterior words')
 need(all(s<=4 and t<=4 for(s,t),c in poly.items() if c>0),'all positive exponent positions');need([(s,t,c) for(s,t),c in sorted(poly.items()) if c<0 and(s<4 or t<4)]==[(5,3,-96)],'unique discarded negative');positive=sum(c*q(6,s)*q(6,t) for(s,t),c in poly.items() if c>0);negative=sum(c*q(6,s)*q(6,t) for(s,t),c in poly.items() if c<0 and s>=4 and t>=4);need((positive,negative)==(2080753736,-2560896000),'exact all-order bound constants')
 Q=(4,13,8,17,7);C3=(2,4,1);W=(10,28,25,16,1);mapping=(3,2,4,0,1);need(all(bool(Q[i]>>j&1)==bool(W[mapping[i]]>>mapping[j]&1) for i in range(5) for j in range(5)),'physical W5 transport');need(not any(module(Q,M) for M in range(1,31) if M.bit_count()>1),'prime quotient')
 requests={};literal={};controls=[]
 def save(T,a,b,am,cm,value):
  key=(a,b,collapsed(T,a,b,am,cm));need(key not in requests or requests[key]==value,'consistent exact weighted value');requests[key]=value
  if len(T)<=14:
   lk=(tuple(T),am,cm);need(lk not in literal or literal[lk]==value,'consistent full-name value');literal[lk]=value
 for rec in data['actual_physical_controls']:
  a,b=rec['a_size'],rec['b_size'];T=blow(Q,[C3,C3,C3,tt(a),tt(b)]);U=blow(op(Q),[C3,C3,C3,tt(a),tt(b)]);L=flip(T,7);need((T,U,L)==tuple(tuple(rec[s]) for s in('source_rows','quotient_target_rows','local_target_rows')),'all physical named graphs');v=sum(c*q(a,s)*q(b,t) for(s,t),c in poly.items());d1=d2=0;need(len(rec['queries'])==16,'all physical finite-difference corners')
  for z in rec['queries']:
   need(z['sign']==(-1)**(z['a'].bit_count()+z['c'].bit_count()),'quartic interpolation sign');save(T,a,b,z['a'],z['c'],z['source']);save(U,a,b,z['a'],z['c'],z['quotient_target']);save(L,a,b,z['a'],z['c'],z['local_target']);d1+=z['sign']*(z['source']-z['local_target']);d2+=z['sign']*(z['source']-z['quotient_target'])
  need(v==d1==d2==rec['actual_coefficient'],'full local and quotient marked coefficients');controls.append([a,b,v])
 bank=data['root_family_whole_child_bank'];T=tuple(bank['source_rows']);blocks=[7,56,448,63<<9,63<<15];need(T==blow(Q,[C3,C3,C3,tt(6),tt(6)]),'family source');battery=list(map(tuple,bank['battery']));quartic=[tuple(x) for x in bank['quartic']];compressed=[2,4]+blocks[1:];anchor=[(B|(1 if T[0]&B else 0),1) for B in compressed]+[(1,0)];need(anchor==list(map(tuple,bank['anchor'])) and len(compressed)==6 and sum(bool(T[0]&B) for B in compressed)==2,'correct exposed-anchor parities');need(set(battery)=={(a,c) for a,c,s in quartic}|set(anchor) and len(battery)==22,'exact fixed22-query union');base=bank['baseline_values'];at={x:i for i,x in enumerate(battery)}
 for(a,c),v in zip(battery,base):save(T,6,6,a,c,v)
 need(len(bank['whole_child_controls'])==32,'all32 whole-child controls');seen=set()
 for rec in bank['whole_child_controls']:
  bits=rec['selected_child_bits'];need(bits not in seen,'unique family');seen.add(bits);U=op(T)
  for j,B in enumerate(blocks):
   if bits>>j&1:U=flip(U,B)
  need(U==tuple(rec['target_rows']),'exact original-module XOR labels');vals=rec['values'];need(vals!=base and vals[at[(0,0)]]==base[at[(0,0)]],'scalar equal and response different')
  for(a,c),v in zip(battery,vals):save(U,6,6,a,c,v)
  if bits&1:need(sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in quartic)==-573454264,'restored A coefficient')
  else:need(2 in[(base[at[x]]-vals[at[x]])%4 for x in anchor],'unrestored A anchor')
 with tempfile.TemporaryDirectory() as td:
  td=Path(td)
  for name in('weighted','literal'):subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/(name+'.cpp')),'-o',str(td/name)],check=True,capture_output=True)
  wire=str(len(requests))+'\n'+'\n'.join('11 '+str(a)+' '+str(b)+' '+' '.join(map(str,T)) for a,b,T in requests)+'\n';out=list(map(int,subprocess.check_output([str(td/'weighted')],input=wire.encode()).split()));need(out==list(requests.values()),'every retained weighted count independently reproduced')
  wire=str(len(literal))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in literal)+'\n';out=list(map(int,subprocess.check_output([str(td/'literal')],input=wire.encode()).split()));need(out==list(literal.values()),'all small physical counts independently reproduced')
 need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired receipts equal');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','complete_exterior_word_kernel_terms':terms,'independent_full_polynomial':[[s,t,c] for(s,t),c in sorted(poly.items())],'positive_bound_base':positive,'retained_negative_bound_base':negative,'actual_quartic_controls':controls,'independent_weighted_actual_counts':len(requests),'independent_full_name_counts':len(literal),'whole_child_root_family_controls':32,'root_battery_queries':22,'scope':'All a,b>=6 sign proof analytically reviewed; finite family controls only whole-child subfamilies, no arbitrary strong outer lift.','source_sha256':pins,'audit_sha256':{x:hashlib.sha256((HERE/x).read_bytes()).hexdigest() for x in('audit.py','weighted.cpp','literal.cpp')}},sort_keys=True,indent=2))
if __name__=='__main__':main()
