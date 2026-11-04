from pathlib import Path
from itertools import combinations,product
import hashlib,json,math,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
TYPES=[('P',[0,1,2,0],0),('E',[0,1,6,4],0),('O',[0,1,6,6],1),('C',[0,1,2,2],1)]
def need(c,m):
 if not c:raise RuntimeError(m)
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def blow(Q,blocks):
 off=[];T=[]
 for B in blocks:off.append(len(T));T.extend([0]*len(B))
 for i,B in enumerate(blocks):
  for u,v in enumerate(B):T[off[i]+u]=(v<<off[i])|sum(((1<<len(C))-1)<<off[j] for j,C in enumerate(blocks) if Q[i]>>j&1)
 return T,off
def physical(a,blocks,off):return sum(((1<<len(B))-1)<<off[i] for i,B in enumerate(blocks) if a>>i&1)
def main():
 reps=json.loads((ROOT.parent/'prime_five_even_profiles/quotients.json').read_text());rs=[reps[0],reps[2]];half=[a for a in range(32) if a.bit_count()<=2];queries=list(product(half,repeat=2));lines=['2']
 for r in rs:lines.extend([' '.join(map(str,r['rows'])),str(len(queries))]+[f'{a} {c}' for a,c in queries])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'covers';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'covers.cpp'),'-o',str(exe)],check=True);ks=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).decode().splitlines()]
 need(len(ks)==1024 and all(len(v)==153 for v in ks),'complete source/target first-three kernels');shapes=[(1,)*5]
 for i in range(5):
  for k in[2,3]:t=[1]*5;t[i]=k;shapes.append(tuple(t))
 for i,j in combinations(range(5),2):
  for k,l in product([2,3],repeat=2):t=[1]*5;t[i]=k;t[j]=l;shapes.append(tuple(t))
 need(len(shapes)==51 and len(set(shapes))==51,'all relevant run shapes');records=[];alllook=[]
 for qi,r in enumerate(rs):
  bank=[];lookup={}
  for j,(a,c) in enumerate(queries):
   T,U=ks[qi*512+2*j:qi*512+2*j+2];bank.append([a,c,T,U]);lookup[a,c]=[T,U]
  states=[]
  for ids in product(range(4),repeat=5):
   if sum(TYPES[i][2] for i in ids)%2==0:continue
   w=[math.prod(TYPES[ids[i]][1][s[i]] for i in range(5))%8 for s in shapes];first=None;full=[]
   for a,c,T,U in bank:
    vv=[sum(v*(T[3*j+h]-U[3*j+h]) for j,v in enumerate(w))%8 for h in range(3)];full.extend(vv)
    if first is None and vv==[4,0,0]:first=[a,c]
   need(first is not None,'every odd-total type assignment has complete first-three separator');states.append({'types':''.join(TYPES[i][0] for i in ids),'query':first,'all256_profile_contrasts_sha256':hashlib.sha256(bytes(full)).hexdigest()})
  need(len(states)==512,'complete odd-total type domain');records.append({'quotient_code':r['code'],'quotient_rows':r['rows'],'kernel_bank_columns':['a','c','source153','target153'],'complete_kernel_bank':bank,'all512_odd_type_certificates':states});alllook.append(lookup)
 # Check the source's singleton-position certificate and middle-class nontrivial mod4 separator.
 local=json.loads((ROOT.parent/'prime_five_all_profiles/normal.json').read_text());singleton=[]
 for r in local['records']:
  need(len(r['singleton_pattern_mod4'])==31 and all(row[1][2]==2 for row in r['singleton_pattern_mod4']),'all singleton patterns have mod4 baseline gap');singleton.append({'code':r['code'],'records':r['singleton_pattern_mod4']})
 middle=next(r for r in local['records'] if r['code']==41);v=next(r for r in middle['full_run_bank'] if r[:2]==[0,1]);T,U=v[2:];need((T[0]-U[0]+2*sum(T[1+i]-U[1+i] for i in range(5)))%4==2,'middle quotient universal nontrivial mod4 query')
 # Actual mixed canonical children in a strong cyclic exterior.
 can=[tt(2),tt(4),tt(3),[2,4,1]];patterns=[(0,0,2,0,0),(3,0,1,0,0)];physicalrecords=[];graphs=[]
 for qi,(r,ids) in enumerate(zip(rs,patterns)):
  name=''.join(TYPES[i][0] for i in ids);a,c=next(s['query'] for s in records[qi]['all512_odd_type_certificates'] if s['types']==name);blocks=[can[i] for i in ids];w=[math.prod(TYPES[ids[i]][1][s[i]] for i in range(5))%8 for s in shapes];unit=pow(3,ids.count(3),8)
  for side in [0,1]:
   Q=r['rows'] if not side else[31^(1<<i)^v for i,v in enumerate(r['rows'])];T,off=blow(Q,blocks);n=len(T);T=[v|(1<<n) for v in T]+[1<<(n+1),(1<<n)-1];aa=physical(a,blocks,off)|(1<<(n+1));cc=physical(c,blocks,off)|(1<<(n+1));K=alllook[qi][a,c][side];core=[unit*sum(v*K[3*j+h] for j,v in enumerate(w))%8 for h in range(3)];expected=(3*core[0]+core[1])%8;graphs.append((T,aa,cc));physicalrecords.append({'quotient_code':r['code'],'types':name,'side':side,'old_rows':T,'query':[aa,cc],'block_query':[a,c],'outer_rows':[2,4,1],'odd_normalization_unit':unit,'core_first_three_mod8':core,'predicted_ambient_mod8':expected})
 lines=[]
 for T,a,c in graphs:lines.extend([f'{len(T)} 1',' '.join(map(str,T)),f'{a} {c}'])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);vv=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(len(vv)==4,'all direct strong-ambient physical counts')
 for r,v in zip(physicalrecords,vv):need(v%8==r['predicted_ambient_mod8'],'physical root-incidence lift and odd scalar');r['actual_ambient_mod16']=v
 need(all((vv[i]-vv[i+1])%8==4 for i in[0,2]),'both actual strong-ambient separators')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','covers.cpp','physical.cpp']};deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in ['prime_five_all_profiles/PROOF.md','prime_five_even_ambient/PROOF.md','mod8_transfer/PROOF.md','module_decomposition/PROOF.md','monotone_dual_queries/PROOF.md']}
 print(json.dumps({'status':'PASS','scope':'Every prime5 quotient with arbitrary nonempty children in every fixed ambient; branches and source-relative families stated separately.','run_shape_order':shapes,'kernel_coordinate_order':['q1','q2','q3'],'records':records,'singleton_pattern_mod4_certificates':singleton,'middle_nontrivial_mod4_query':[0,1],'physical_ambient_controls':physicalrecords,'source_pins':pins,'dependency_pins':deps},sort_keys=True,indent=2))
if __name__=='__main__':main()
