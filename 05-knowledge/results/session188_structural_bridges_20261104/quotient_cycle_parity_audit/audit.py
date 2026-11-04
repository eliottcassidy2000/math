"""Independent root-zero cycle DP and complete coefficient/physical review."""
from pathlib import Path
from itertools import product
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent;SOURCE=BASE/'quotient_cycle_parity'
def need(x,s):
 if not x:raise RuntimeError(s)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(T):return [((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T)]
def blow(Q,B):
 o=[];n=0
 for b in B:o.append(n);n+=len(b)
 T=[0]*n
 for i,b in enumerate(B):
  for v,r in enumerate(b):
   T[o[i]+v]=r<<o[i]
   for j,c in enumerate(B):
    if Q[i]>>j&1:T[o[i]+v]|=((1<<len(c))-1)<<o[j]
 return T,o
def aug(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  for w,bits in((n,a),(n+1,c)):
   if bits>>v&1:R[v]|=1<<w
   else:R[w]|=1<<v
 return tuple(R)
def mask(a,B,o):return sum(((1<<len(b))-1)<<s for i,(b,s) in enumerate(zip(B,o)) if a>>i&1)
def profile(T,full=True):
 n=len(T);N=1<<n;d=[[0]*n for _ in range(N)];H=[0]*N
 for M in range(1,N):
  for v in range(n):
   if M>>v&1:
    S=M^(1<<v);d[M][v]=sum(d[S][u] for u in range(n) if T[u]>>v&1) if S else 1
  H[M]=sum(d[M])
 if not full:return H[-1]
 q=[[0]*N for _ in range(n+1)];q[0][0]=1
 for h in range(1,n+1):
  for M in range(1,N):
   S=M
   while S:q[h][M]+=H[S]*q[h-1][M^S];S=(S-1)&M
 return [q[h][-1] for h in range(n+1)]
def main():
 raw=(SOURCE/'normal.json').read_bytes();need(raw==(SOURCE/'optimized.json').read_bytes(),'producer modes');j=json.loads(raw)
 for f,h in j['source_pins'].items():need(sha(SOURCE/f)==h,'source pins')
 for f,h in j['dependency_pins'].items():need(sha(BASE/f)==h,'dependency pins')
 data=j['cycle_graphs'];need(len({tuple(r['rows']) for r in data})==len(data)==j['distinct_literal_cycle_graphs'],'unique full cycle graph set')
 with tempfile.TemporaryDirectory(prefix='quotient-cycle-audit-') as td:
  exe=Path(td)/'cycles';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'cycles.cpp'),'-o',str(exe)],check=True);wire=str(len(data))+'\n'+'\n'.join(str(len(r['rows']))+' '+' '.join(map(str,r['rows'])) for r in data)+'\n';z=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(z==[r['Hamilton_cycles'] for r in data],'every exact cycle count independently rooted at zero')
 Z={tuple(r['rows']):v for r,v in zip(data,z)};ham=json.loads((BASE/'coarse_mod8_affine'/'normal.json').read_text());two=json.loads((BASE/'coarse_two_cover_affine'/'normal.json').read_text());expected={}
 for h,t in zip(ham['small_model_controls'],two['small_model_controls']):
  need(h['quotient_rows']==t['quotient_rows'],'matching complete model domains')
  for a,b in zip(h['extracted_coefficients'],t['extracted_coefficients']):expected[tuple(h['quotient_rows']),*a['block_query']]=(a['coefficients'],b['coefficients'])
 relations=0
 for r in j['all_coefficient_controls']:
  Q=r['quotient_rows'];a,c=r['block_query'];need(expected.pop((tuple(Q),a,c))==(r['Hamilton_coefficients'],r['two_cover_coefficients']),'complete exact coefficient payload');zz=[Z[aug(P,a,c)] for P in(Q,op(Q))];need(zz==r['bare_cycle_counts'],'bare physical masks')
  for cyc,H,B in zip(zz,r['Hamilton_coefficients'],r['two_cover_coefficients']):
   need((B['B0']-2*H['K0']+4+4*(len(Q)%2)*cyc)%8==0,'base cycle identity')
   for k,b,c0 in zip(H['K2'],B['B2'],B['B3']):need((b-2*k-2*cyc)%4==0 and(c0-2*cyc)%4==0,'both coefficient identities');relations+=1
 need(not expected,'no coefficient rows omitted')
 larger={(r['name'],tuple(q)):r for r in two['larger_child_controls'] for q in r['block_queries']};canonical={tuple(r['types']):r for r in two['canonical_controls']}
 for r in j['physical_cycle_collapse_controls']:
  a,c=r['block_query'];B=r['child_rows'];Q=[2,4,1] if r['kind']=='canonical' else larger[r['kind'],(a,c)]['quotient_rows'];T,o=blow(Q,B);U,_=blow(op(Q),B);aa=mask(a,B,o);cc=mask(c,B,o);bare=[Z[aug(P,a,c)] for P in(Q,op(Q))];actual=[Z[aug(P,aa,cc)] for P in(T,U)];need(bare==r['bare_cycle_counts'] and actual==r['physical_cycle_counts'] and all((x-y)%2==0 for x,y in zip(bare,actual)),'complete literal cycle collapse');epsilon=sum(map(len,B))%2;unit=r['odd_H_product_mod8'];vv=r['physical_H_q2']
  for side,cyc in enumerate(bare):H,q=vv[2*side:2*side+2];need((q*pow(unit,-1,8)-2*H*pow(unit,-1,8)+4+4*epsilon*cyc)%8==0,'literal normalized two-cover identity')
 neg=j['odd_total_negative_control'];Q=[2,4,1];B=neg['child_rows'];T,o=blow(Q,B);U,_=blow(op(Q),B);a,c=neg['physical_query'];pro=[profile(aug(P,a,c)) for P in(T,U)];need(pro==neg['full_core_profiles'],'complete negative core profiles');ext=1<<(len(T)+1);ambient=[]
 for P in(T,U):R,_=blow(Q,[P,[0],[0]]);ambient.append(profile(aug(R,a|ext,c|ext),False))
 need(ambient==neg['ambient_H_counts'] and ambient[0]-ambient[1]==neg['ambient_gap']==-600,'actual hostile ambient');need(neg['bare_cycle_counts']==[1,2] and [pro[0][k]-pro[1][k] for k in(1,2,3)]==[-68,-396,-1116],'odd-order correction boundary')
 return dict(status='PASS',source_sha256={f:sha(SOURCE/f) for f in('PROOF.md','verify.py','cycles.cpp','normal.json','optimized.json')},dependency_pins=j['dependency_pins'],checker_sha256=sha(Path(__file__)),counter_sha256=sha(ROOT/'cycles.cpp'),distinct_exact_cycle_graphs=len(data),coefficient_relations=relations,complete_actual_collapse_controls=len(j['physical_cycle_collapse_controls']),negative_complete_core_profiles=pro,negative_exact_ambient_counts=ambient)
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
