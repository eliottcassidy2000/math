"""Independent predecessor physical DP and two-path allocation audit."""
from pathlib import Path
from itertools import product,permutations
from functools import lru_cache
import hashlib,json,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'coarse_two_cover_affine'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tt(n):return [sum(1<<v for v in range(u+1,n)) for u in range(n)]
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 T=[0]*n
 for i,b in enumerate(B):
  for v,r in enumerate(b):
   T[off[i]+v]=r<<off[i]
   for j,c in enumerate(B):
    if Q[i]>>j&1:T[off[i]+v]|=((1<<len(c))-1)<<off[j]
 return T,off
def augment(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  for w,mask in((n,a),(n+1,c)):
   if mask>>v&1:R[v]|=1<<w
   else:R[w]|=1<<v
 return tuple(R)
def physical(mask,B,off):return sum(((1<<len(b))-1)<<o for i,(b,o) in enumerate(zip(B,off)) if mask>>i&1)
@lru_cache(None)
def paths(T,k):
 if not any(k):return ()
 result=[]
 def rec(w,left):
  if not any(left):result.append(w);return
  for v,z in enumerate(left):
   if z and(not w or T[w[-1]]>>v&1):u=list(left);u[v]-=1;rec(w+(v,),tuple(u))
 rec((),k);return tuple(result)
def pairs(T,k):
 cap=tuple(k)+(1,1);out=set()
 for a in product(*(range(z+1) for z in cap)):
  if not any(a) or a==cap:continue
  b=tuple(x-y for x,y in zip(cap,a))
  for P in paths(T,a):
   for Q in paths(T,b):out.add((P,Q))
 return out

def main():
 raw=(SOURCE/'normal.json').read_bytes();need(raw==(SOURCE/'optimized.json').read_bytes(),'source paired exact');j=json.loads(raw)
 for f,h in j['source_pins'].items():need(sha(SOURCE/f)==h,'source pin')
 for f,h in j['dependency_pins'].items():need(sha(ROOT.parent/f)==h,'dependency pin')
 graphs=[];lookup={};expected=[];occurrences=0
 def add(T,a,c,H,q2):
  nonlocal occurrences
  R=augment(T,a,c);occurrences+=1
  if R in lookup:need(expected[lookup[R]]==[H,q2],'duplicate actual graph data')
  else:lookup[R]=len(graphs);graphs.append(R);expected.append([H,q2])
 model_count=0
 for r in j['small_model_controls']:
  Q=r['quotient_rows'];m=len(Q);records=r['extracted_coefficients']
  for model in r['models']:
   B=[tt(n) for n in model['sizes']];T,o=blow(Q,B);U,_=blow(op(Q),B)
   for (a,c),v in zip(r['block_queries'],model['values_H_q2_source_target']):aa=physical(a,B,o);cc=physical(c,B,o);add(T,aa,cc,*v[:2]);add(U,aa,cc,*v[2:]);model_count+=1
  for qi,rec in enumerate(records):
   for side,z in enumerate(rec['coefficients']):
    g=[md['values_H_q2_source_target'][qi][2*side+1] for md in r['models']];need(z['B0']==g[0] and not g[0]%2,'base coefficient')
    for i in range(m):need(z['B2'][i]*2==g[1+2*i]-g[0] and z['B3'][i]*6==g[2+2*i]-3*g[1+2*i]+2*g[0] and z['B2'][i]%2==z['B3'][i]%2==0,'exact model extraction')
 canon=0
 def predicted(rec,types,unit):
  return [unit*(z['B0']+2*sum((t!='S')*b+(t in('O','C'))*c for t,b,c in zip(types,z['B2'],z['B3'])))%8 for z in rec]
 coef={(tuple(r['quotient_rows']),*v['block_query']):v['coefficients'] for r in j['small_model_controls'] for v in r['extracted_coefficients']}
 for r in j['canonical_controls']:
  Q=[2,4,1];B=r['child_rows'];T,o=blow(Q,B);U,_=blow(op(Q),B);unit=pow(3,r['types'].count('C'),8);need(unit==r['odd_H_product_mod8'],'canonical actual odd unit')
  for (a,c),(aa,cc),v,p in zip(product(range(8),repeat=2),r['physical_queries'],r['values_H_q2_source_target'],r['predicted_q2_source_target_mod8']):
   need((aa,cc)==(physical(a,B,o),physical(c,B,o)),'physical named block mask');add(T,aa,cc,*v[:2]);add(U,aa,cc,*v[2:]);need(predicted(coef[tuple(Q),a,c],r['types'],unit)==p==[v[1]%8,v[3]%8],'all canonical affine values');canon+=2
 for r in j['larger_child_controls']:
  Q=r['quotient_rows'];B=r['child_rows'];T,o=blow(Q,B);U,_=blow(op(Q),B);unit=1
  for q in r['profiles']:unit=unit*q[1]%8
  need(unit==r['odd_H_product_mod8'],'actual larger child odd product')
  for (a,c),(aa,cc),v,p in zip(r['block_queries'],r['physical_queries'],r['values_H_q2_source_target'],r['predicted_q2_mod8']):
   need((aa,cc)==(physical(a,B,o),physical(c,B,o)),'larger physical masks');add(T,aa,cc,*v[:2]);add(U,aa,cc,*v[2:]);need(predicted(coef[tuple(Q),a,c],r['types'],unit)==p==[v[1]%8,v[3]%8],'all actual child affine values')
 with tempfile.TemporaryDirectory(prefix='q2-affine-audit-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T)) for T in graphs)+'\n';v=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need([v[i:i+2] for i in range(0,len(v),2)]==expected,'every complete actual H and q2 value')
 total=boundary=0
 for rec in j['path_swap_instances']:
  T=augment(rec['quotient_rows'],*rec['block_query']);ww=pairs(T,rec['visits']);need(len(ww)==rec['count'] and len(ww)%2==0,'independent path allocation count')
  for P,Q in ww:need(P!=Q and (Q,P) in ww,'exact fixed-point-free swap');boundary+=int(P[-1]==Q[0])
  total+=len(ww)
 need(total==j['total_swapped_word_pairs'] and boundary>0,'all swap data and cut boundary')
 for rec in j['positive_word_witnesses'].values():need(pairs(augment(rec['quotient_rows'],*rec['block_query']),rec['visits'])=={tuple(map(tuple,w)) for w in rec['words']},'complete retained word witness')
 # Independent assignment to three slots, including identical port-free singleton paths.
 T=augment([0],0,0);cap=(3,1,1);three=set()
 for aa in product(*(range(z+1) for z in cap)):
  if not any(aa):continue
  left=tuple(z-a for z,a in zip(cap,aa))
  for bb in product(*(range(z+1) for z in left)):
   cc=tuple(z-b for z,b in zip(left,bb))
   if not any(bb) or not any(cc):continue
   for P in paths(T,aa):
    for Q in paths(T,bb):
     for R in paths(T,cc):three.add((P,Q,R))
 need(len(three)==9 and three=={tuple(map(tuple,w)) for w in j['three_cover_boundary']['all_nine_abstract_covers']},'all nine three-cover boundary words')
 return dict(status='PASS',source_sha256={f:sha(SOURCE/f) for f in('PROOF.md','verify.py','count.cpp','normal.json','optimized.json')},dependency_pins=j['dependency_pins'],checker_sha256=sha(Path(__file__)),independent_counter_sha256=sha(ROOT/'count.cpp'),actual_physical_occurrences=occurrences,distinct_actual_graphs=len(graphs),canonical_individual_q2_values=canon,all_model_extractions_checked=True,explicit_path_swap_instances=len(j['path_swap_instances']),explicit_word_pairs=total,same_label_cut_words=boundary,three_cover_boundary_count=len(three))
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
