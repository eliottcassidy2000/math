from pathlib import Path
from itertools import product,combinations
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
TYPES=[('S',0,0,0),('P',2,0,0),('E',6,4,0),('O',6,6,1),('C',2,2,1)];PAIRS=list(combinations(range(5),2))
def need(c,m):
 if not c:raise RuntimeError(m)
def subsets(a):
 s=a
 while True:
  yield s
  if not s:break
  s=(s-1)&a

def weights(ids):
 ts=[TYPES[i] for i in ids];u=[t[1] for t in ts];return [1]+u+[2*t[3] for t in ts]+[4*bool(u[i])*bool(u[j]) for i,j in PAIRS]
def value(v,w):return sum(x*y for x,y in zip(v,w))%8
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def blow(Q,blocks):
 off=[];T=[]
 for B in blocks:off.append(len(T));T.extend([0]*len(B))
 for i,B in enumerate(blocks):
  for u,v in enumerate(B):T[off[i]+u]=(v<<off[i])|sum(((1<<len(C))-1)<<off[j] for j,C in enumerate(blocks) if Q[i]>>j&1)
 return T,off
def phys(a,blocks,off):return sum(((1<<len(B))-1)<<off[i] for i,B in enumerate(blocks) if a>>i&1)
def main():
 reps=json.loads((ROOT.parent/'prime_five_even_profiles/quotients.json').read_text());need([(r['code'],r['orbit_size']) for r in reps]==[(10,120),(41,120),(76,24)],'audited complete quotient representatives');inp='3\n'+'\n'.join(' '.join(map(str,r['rows'])) for r in reps)+'\n'
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'run';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'run_bank.cpp'),'-o',str(exe)],check=True);ls=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=inp.encode()).decode().splitlines()]
 need(len(ls)==6144 and all(len(v)==21 for v in ls),'full run-channel bank');half=[a for a in range(32) if a.bit_count()<=2];records=[];lookups=[];alltypes=list(product(range(5),repeat=5))
 for qi,r in enumerate(reps):
  bank=[];values={}
  for a in range(32):
   for c in range(32):
    T,U=ls[qi*2048+2*(32*a+c):qi*2048+2*(32*a+c)+2];need(all(v%2==0 for v in T[6:11]+U[6:11]),'all single triple-run coefficients even');bank.append([a,c,T,U]);values[a,c]=[T,U]
  coef={};cb=[]
  for a,c in product(half,repeat=2):
   vv=[[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*values[x,y][s][j] for x in subsets(a) for y in subsets(c))%16 for j in range(21)] for s in [0,1]];coef[a,c]=vv;cb.append([a,c,*vv])
  for a,c in product(range(32),repeat=2):
   for side in [0,1]:
    vv=[sum(v[side][j] for (x,y),v in coef.items() if x&a==x and y&c==y)%16 for j in range(21)];need(vv==values[a,c][side],'complete coarse cubic/quartic reconstruction for all channels')
  gaps=[[ (x-y)%8 for x,y in zip(T,U)] for a,c,T,U in bank];states=[];blind=[]
  for ids in alltypes:
   w=weights(ids);vv=[value(v,w) for v in gaps];first=next((i for i,z in enumerate(vv) if z),None);names=''.join(TYPES[i][0] for i in ids)
   if first is None:blind.append(names)
   states.append({'types':names,'nonzero_queries':sum(bool(x) for x in vv),'first_separator':None if first is None else[bank[first][0],bank[first][1],vv[first]],'full_contrast_sha256':hashlib.sha256(bytes(vv)).hexdigest()})
  need(blind==([] if qi<2 else['EEEEE']),'all3125 exact normalized type assignments')
  singleton=[]
  for M in range(1,32):
   found=None
   for a,c,T,U in bank:
    gap=(T[0]-U[0]+2*sum(T[1+i]-U[1+i] for i in range(5) if not M>>i&1))%4
    if gap:found=[a,c,gap];break
   need(found is not None,'all singleton-position patterns have mod4 baseline separator');singleton.append([M,found])
  records.append({**r,'full_run_bank_columns':['a','c','source21','target21'],'full_run_bank':bank,'complete_coefficient_bank':cb,'type_classification':states,'blind_types':blind,'singleton_pattern_mod4':singleton});lookups.append(values)
 can=[tt(1),tt(2),tt(4),tt(3),[2,4,1]];cases=[];physical=[]
 for qi,r in enumerate(reps):
  for ids in [(3,)*5,(4,)*5,(2,1,3,4,0)]:
   names=''.join(TYPES[i][0] for i in ids);state=next(s for s in records[qi]['type_classification'] if s['types']==names);a,c,_=state['first_separator'];blocks=[can[i] for i in ids];unit=pow(3,ids.count(4),8);w=weights(ids)
   for side in [0,1]:
    Q=r['rows'] if side==0 else[31^(1<<i)^v for i,v in enumerate(r['rows'])];T,off=blow(Q,blocks);qq=[phys(a,blocks,off),phys(c,blocks,off)];cases.append((T,qq));physical.append({'quotient_code':r['code'],'types':names,'side':side,'old_rows':T,'block_query':[a,c],'physical_query':qq,'odd_normalization_unit':unit,'predicted_mod8':unit*value(lookups[qi][a,c][side],w)%8})
 lines=[]
 for T,(a,c) in cases:lines.extend([f'{len(T)} 1',' '.join(map(str,T)),f'{a} {c}'])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);vals=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(len(vals)==len(physical),'all literal counts');
 for r,v in zip(physical,vals):need(v%8==r['predicted_mod8'],'actual mixed-child physical count and scalar restoration');r['actual_mod16']=v
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','run_bank.cpp','physical.cpp']};deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in ['five_profile_types/PROOF.md','prime_five_even_profiles/PROOF.md','parity_profile_transfer/PROOF.md']}
 print(json.dumps({'status':'PASS','scope':'Local prime-five quotient-converse separation for all nonempty child tournaments. No ambient/family inference in this packet.','channel_order':['K0']+[f'K{i}' for i in range(5)]+[f'L{i}' for i in range(5)]+[f'K{i}{j}' for i,j in PAIRS],'type_order':[t[0] for t in TYPES],'complete_states_per_quotient':3125,'records':records,'literal_physical_controls':physical,'source_pins':pins,'dependency_pins':deps},sort_keys=True,indent=2))
if __name__=='__main__':main()
