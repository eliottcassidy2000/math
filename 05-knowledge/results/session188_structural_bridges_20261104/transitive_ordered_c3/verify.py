from pathlib import Path
from itertools import combinations
import hashlib,json,subprocess,tempfile
from algebra import expand,value,q,need
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;Q=(4,13,8,17,7);C3=(2,4,1)
def tt(n):return tuple(((1<<n)-1)^((1<<(u+1))-1) for u in range(n))
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def extend(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:U[u]|=1<<n
  else:U[n]|=1<<u
  if c>>u&1:U[u]|=1<<(n+1)
  else:U[n+1]|=1<<u
 return tuple(U)
def collapse(T,lo,hi,am,cm):
 n=len(T);V=extend(T,am,cm);masks=[((1<<lo)-1)<<9,((1<<hi)-1)<<(9+lo)];need(all(module(V,M) for M in masks),'colored chains homogeneous in actual query');order=list(range(9))+[n,n+1,9,9+lo]
 return tuple(sum(1<<j for j,v in enumerate(order) if V[u]>>v&1) for u in order)
def main():
 formula=expand();need(all(q(n,s)*q(6,4)<=q(6,s)*q(n,4) for n in range(6,31) for s in range(5)),'finite controls of inherited low Stirling ratios');need(all(q(n,s)*q(6,4)>=q(6,s)*q(n,4) for n in range(6,31) for s in range(4,7)),'finite controls of inherited high Stirling ratios');signs=[[a,b,value(a,b)] for a in range(1,21) for b in range(1,21)];need(all(z<0 for a,b,z in signs if a>=6 and b>=6),'all-order theorem finite controls')
 W=(10,28,25,16,1);mapping=(3,2,4,0,1);need(all(bool(Q[i]>>j&1)==bool(W[mapping[i]]>>mapping[j]&1) for i in range(5) for j in range(5)),'physical W5 label map');need(not any(module(Q,M) for M in range(1,31) if M.bit_count()>1),'prime quotient directly')
 requests=[];ids={};literal={};records=[]
 def add(T,lo,hi,a,c):
  key=(lo,hi,collapse(T,lo,hi,a,c))
  if key not in ids:ids[key]=len(requests);requests.append(key)
  if len(T)<=14:literal[(lo,hi,tuple(T),a,c)]=ids[key]
  return ids[key]
 def quartet():return [(1*(bool(b&1))+8*(bool(b&2)),2*(bool(b&4))+64*(bool(b&8)),(-1)**b.bit_count()) for b in range(16)]
 for lo,hi in[(1,1),(2,2),(3,2),(4,4),(5,7),(6,6),(6,7),(7,6),(8,8)]:
  T,_=blow(Q,[C3,C3,C3,tt(lo),tt(hi)]);U,_=blow(op(Q),[C3,C3,C3,tt(lo),tt(hi)]);L=flip(T,7);rec={'a_size':lo,'b_size':hi,'source_rows':T,'local_target_rows':L,'quotient_target_rows':U,'predicted_coefficient':value(lo,hi),'queries':[]}
  for a,c,s in quartet():rec['queries'].append({'a':a,'c':c,'sign':s,'source_id':add(T,lo,hi,a,c),'local_id':add(L,lo,hi,a,c),'quotient_id':add(U,lo,hi,a,c)})
  records.append(rec)
 lo=hi=6;T,off=blow(Q,[C3,C3,C3,tt(lo),tt(hi)]);n=len(T);blocks=[7,56,448,63<<9,63<<15];compressed=[2,4]+blocks[1:];d=sum(bool(T[0]&B) for B in compressed);need(len(compressed)==6 and d==2,'exact anchor parities');anchor=[(B|(1 if T[0]&B else 0),1) for B in compressed]+[(1,0)];qu=quartet();battery=list(dict.fromkeys([(a,c) for a,c,s in qu]+anchor));need(len(battery)==22,'exact22-query battery');bank={'source_rows':T,'blocks':blocks,'battery':battery,'quartic':qu,'anchor':anchor,'baseline_ids':[add(T,lo,hi,a,c) for a,c in battery],'whole_child_controls':[]}
 for bits in range(32):
  U=op(T)
  for t,M in enumerate(blocks):
   if bits>>t&1:U=flip(U,M)
  bank['whole_child_controls'].append({'selected_child_bits':bits,'target_rows':U,'value_ids':[add(U,lo,hi,a,c) for a,c in battery]})
 with tempfile.TemporaryDirectory(prefix='transitive-ordered-c3-') as td:
  td=Path(td)
  for name in('run_count','literal'):subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/(name+'.cpp')),'-o',str(td/name)],check=True,capture_output=True)
  wire=str(len(requests))+'\n'+'\n'.join(f'11 {lo} {hi} '+' '.join(map(str,T)) for lo,hi,T in requests)+'\n';vals=list(map(int,subprocess.check_output([str(td/'run_count')],input=wire.encode()).split()));need(len(vals)==len(requests),'all run-word outputs')
  graphs=[extend(T,a,c) for lo,hi,T,a,c in literal];wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';actual=list(map(int,subprocess.check_output([str(td/'literal')],input=wire.encode()).split()));need(actual==[vals[t] for t in literal.values()],'independent expanded-name subset DP controls')
 for rec in records:
  loc=quo=0
  for z in rec['queries']:z['source']=vals[z.pop('source_id')];z['local_target']=vals[z.pop('local_id')];z['quotient_target']=vals[z.pop('quotient_id')];loc+=z['sign']*(z['source']-z['local_target']);quo+=z['sign']*(z['source']-z['quotient_target'])
  need(loc==quo==rec['predicted_coefficient'],'actual complete physical quartic against exact polynomial');rec['actual_coefficient']=loc
 base=bank['baseline_values']=[vals[t] for t in bank.pop('baseline_ids')];at={q:i for i,q in enumerate(battery)}
 for f in bank['whole_child_controls']:
  z=f['values']=[vals[t] for t in f.pop('value_ids')];need(z!=base and z[at[(0,0)]]==base[at[(0,0)]],'whole-child root-family scalar equality and detection')
  if f['selected_child_bits']&1:
   delta=sum(s*(base[at[(a,c)]]-z[at[(a,c)]]) for a,c,s in qu);need(delta==value(6,6)<0,'restored A negative quartic');f['quartic_gap']=delta
  else:f['anchor_residues']=[(base[at[q]]-z[at[q]])%4 for q in anchor];need(2 in f['anchor_residues'],'unrestored A anchor')
 paths=[ROOT/n for n in('PROOF.md','verify.py','algebra.py','run_count.cpp','literal.cpp')]+[BASE/d/'PROOF.md' for d in('ordered_cyclic_module','cyclic_child_disjoint','covered_quotient','module_decomposition','alternating_prime_quotient')]+[BASE/'returns_next'/'UNIFORM.md']
 print(json.dumps({'status':'PASS','coefficient_polynomial':formula,'quotient':Q,'W5_map':mapping,'sign_controls':signs,'actual_physical_controls':records,'root_family_whole_child_bank':bank,'distinct_exact_run_word_counts':len(requests),'expanded_name_subset_dp_controls':len(literal),'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
