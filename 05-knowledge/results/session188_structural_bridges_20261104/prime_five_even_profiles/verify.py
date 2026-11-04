from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,tempfile,sys
ROOT=Path(__file__).resolve().parent
def need(c,m):
 if not c:raise RuntimeError(m)
def sub(a):
 s=a
 while True:
  yield s
  if not s:break
  s=(s-1)&a

def classify():
 pairs=list(combinations(range(5),2));ps=list(permutations(range(5)))
 def rows(code):
  T=[0]*5
  for b,(i,j) in enumerate(pairs):
   u,v=(i,j) if code>>b&1 else(j,i);T[u]|=1<<v
  return T
 def prime(T):
  return not any(all(T[v]&M in(0,M) for v in range(5) if not M>>v&1) for M in range(32) if 2<=M.bit_count()<5)
 codes={c for c in range(1024) if prime(rows(c))};left=set(codes);reps=[]
 while left:
  c=min(left);T=rows(c);orb={sum(bool(T[p[i]]>>p[j]&1)<<b for b,(i,j) in enumerate(pairs)) for p in ps};need(orb<=left,'disjoint complete prime orbits');left-=orb;reps.append({'code':c,'rows':T,'orbit_size':len(orb)})
 need(len(codes)==264 and reps==json.loads((ROOT/'quotients.json').read_text()),'all labelled prime5 tournaments covered');return reps

def main():
 reps=classify();inp=str(len(reps))+'\n'+'\n'.join(' '.join(map(str,r['rows'])) for r in reps)+'\n'
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'run';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'run_bank.cpp'),'-o',str(exe)],check=True);lines=[list(map(int,s.split())) for s in subprocess.check_output([str(exe)],input=inp.encode()).decode().splitlines()]
 need(len(lines)==6144 and all(len(v)==16 for v in lines),'complete full bank dimensions');batteries=[[(0,2)],[(0,1)],[(0,1),(2,18)]];half=[i for i in range(32) if i.bit_count()<=2];records=[];at=0
 for rid,r in enumerate(reps):
  values={};bank=[]
  for a in range(32):
   for c in range(32):
    T,U=lines[at:at+2];at+=2;need(all(T[1+3*i+j]==U[1+3*i+j]==0 for i in range(5) for j in [1,2]),'q3/q5 corrections vanish separately');values[a,c]=[T,U];bank.append([a,c,T,U])
  coef={};cb=[]
  for a in half:
   for c in half:
    v=[[sum((-1)**((a^x).bit_count()+(c^y).bit_count())*values[x,y][side][j] for x in sub(a) for y in sub(c))%(16 if j==0 else 2) for j in range(16)] for side in [0,1]];coef[a,c]=v;cb.append([a,c,*v])
  for a in range(32):
   for c in range(32):
    for side in [0,1]:
     vv=[sum(v[side][j] for (x,y),v in coef.items() if x&a==x and y&c==y)%(16 if j==0 else 2) for j in range(16)];need(vv==values[a,c][side],'all1024 values reconstructed by complete coarse cubic/quartic bank')
  states=[]
  for bits in range(32):
   gap=lambda a,c:(values[a,c][0][0]-values[a,c][1][0]+8*sum(((bits>>i)&1)*(values[a,c][0][1+3*i]-values[a,c][1][1+3*i]) for i in range(5)))%16
   full=[gap(a,c) for a in range(32) for c in range(32)];selected=[gap(a,c) for a,c in batteries[rid]];need(any(selected),'fixed battery covers every formal state');states.append({'epsilon_bits':bits,'nonzero_full_queries':sum(bool(x) for x in full),'battery_contrasts':selected,'full_contrast_sha256':hashlib.sha256(bytes(full)).hexdigest()})
  expected=[[(12,[1,0,1,0,1])],[(6,[1,1,1,0,0])],[(0,[0,0,1,1,0]),(8,[0,0,1,1,0])]][rid]
  formulas=[]
  for (a,c),(b,ds) in zip(batteries[rid],expected):
   T,U=values[a,c];need((T[0]-U[0])%16==b and [(T[1+3*i]-U[1+3*i])%2 for i in range(5)]==ds,'four explicit affine formulas');formulas.append({'query':[a,c],'base':b,'epsilon_coefficients_mod2':ds})
  records.append({**r,'affine_channel_order':['TT4_base']+[f'child{i}_q{k}_perturbation' for i in range(5) for k in [2,3,5]],'full_bank_columns':['a','c','source_channels','target_channels'],'complete_full_bank':bank,'coefficient_bank_columns':['a_support','c_support','source_coefficients','target_coefficients'],'complete_coarse_coefficient_bank':cb,'battery_formulas':formulas,'all32_state_controls':states})
 physical=json.loads(subprocess.check_output([sys.executable,str(ROOT/'physical_controls.py')]))
 for r in physical['records']:
  ix=next(i for i,q in enumerate(reps) if q['rows']==r['quotient_rows'] or [31^(1<<j)^v for j,v in enumerate(q['rows'])]==r['quotient_rows']);side=int(r['reverse']);lookup={(a,c):[T,U] for a,c,T,U in records[ix]['complete_full_bank']}
  need(r['values_mod16']==[lookup[tuple(q)][side][0] for q in r['block_queries']],'eight independent actual22-name counts')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','INTERPOLATION.md','verify.py','run_bank.cpp','physical.cpp','physical_controls.py','quotients.json']};deps={f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in ['even_cyclic_mod16/PROOF.md','parity_profile_transfer/PROOF.md']}
 print(json.dumps({'status':'PASS','scope':'Local prime-five quotient-converse separation with arbitrary even-order>=4 children; complete formal domain, no automatic mod16 ambient lift.','prime_labeled_quotients':264,'formal_profiles_per_child':8,'formal_assignments_per_quotient':32768,'remaining_bits_per_quotient':5,'records':records,'physical_controls':physical,'source_pins':pins,'dependency_pins':deps},sort_keys=True,indent=2))
if __name__=='__main__':main()
