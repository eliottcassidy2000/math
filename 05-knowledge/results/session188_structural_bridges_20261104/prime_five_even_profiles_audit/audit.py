#!/usr/bin/env python3
"""Independent complete physical-interpolation audit of prime5 affine banks."""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_five_even_profiles'
PINS={'PROOF.md':'102dc7b512c56a44586725d73dbd88af35a2f70db035103e462b6c33dc303673','INTERPOLATION.md':'289a96298108c3a49fa8efa7f8143330d86f35cf9156d9a31b801401c16b86d6','verify.py':'f0f9d313b8274449cb1e8bfa797bba060a8f10b75f0bf635db38b95b7c1fa49d','run_bank.cpp':'641110e85af5b040eef83f55cbc24e982edb1f4a1417f44991e18b2e4ef6ee82','normal.json':'af99f39f941eb72cdd4e48fec2d04da5f31ebf541989d89313f2da3b628a9abe','optimized.json':'af99f39f941eb72cdd4e48fec2d04da5f31ebf541989d89313f2da3b628a9abe'}
def need(x,msg):
 if not x:raise RuntimeError(msg)
PAIRS=list(combinations(range(5),2))
def rows(c):
 T=[0]*5
 for k,(i,j) in enumerate(PAIRS):a,b=(i,j) if c>>k&1 else(j,i);T[a]|=1<<b
 return T
def op(T):return [31^(1<<i)^r for i,r in enumerate(T)]
def classify():
 groups={};perms=list(permutations(range(5)))
 for code in range(1024):
  T=rows(code)
  if any(all(T[v]&M in(0,M) for v in range(5) if not M>>v&1) for M in range(1,31) if 2<=M.bit_count()<5):continue
  canon=min(sum(bool(T[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(PAIRS)) for p in perms);groups.setdefault(canon,[]).append(code)
 need(sorted(groups)==[10,41,76] and list(map(len,groups.values()))==[120,120,24],'complete264 labelled prime quotient classification');return groups

def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen source '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['source_pins'].items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'source '+f)
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((SRC.parent/f).read_bytes()).hexdigest()==h,'dependency '+f)
 groups=classify();reps=[rows(c) for c in sorted(groups)];need([(r['code'],r['rows'],r['orbit_size']) for r in d['records']]==[(c,rows(c),len(groups[c])) for c in sorted(groups)],'complete physical representative convention')
 with tempfile.TemporaryDirectory(prefix='prime5-physical-interpolation-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'interpolate.cpp'),'-o',str(exe)],check=True);wire='3\n'+'\n'.join(' '.join(map(str,Q)) for Q in reps)+'\n';proc=subprocess.run([str(exe)],input=wire.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);physical=[list(map(int,line.split())) for line in proc.stdout.decode().splitlines()];need(proc.stderr.decode().strip()=='physical_counts 374784 max_names 11','exact physical enumeration accounting')
 need(len(physical)==6144 and all(len(x)==16 for x in physical),'complete independent affine bank');allcoeff=[];allstates=[];offset=0
 for rid,r in enumerate(d['records']):
  bank=r['complete_full_bank'];need(len(bank)==1024,'full1024 query bank');values=[]
  for idx,(a,c,T,U) in enumerate(bank):need(a==idx//32 and c==idx%32 and [T,U]==physical[offset+2*idx:offset+2*idx+2],'every full affine vector equals literal interpolation');values.append([T,U])
  offset+=2048;coeffs=[]
  # Full ten-variable Mobius transform, including all forbidden high degrees.
  for side in(0,1):
   co=[x[side][:] for x in values]
   for bit in range(10):
    for mask in range(1024):
     if mask>>bit&1:
      for j in range(16):co[mask][j]=(co[mask][j]-co[mask^(1<<bit)][j])%(16 if j==0 else 2)
   need(all(all(v==0 for v in co[mask]) for mask in range(1024) if (mask>>5).bit_count()>2 or (mask&31).bit_count()>2),'no hidden higher-degree data')
   re=[x[:] for x in co]
   for bit in range(10):
    for mask in range(1024):
     if mask>>bit&1:
      for j in range(16):re[mask][j]=(re[mask][j]+re[mask^(1<<bit)][j])%(16 if j==0 else 2)
   need(re==[x[side] for x in values],'complete full-vector zeta reconstruction');coeffs.append(co)
  half=[x for x in range(32) if x.bit_count()<=2];expected=[]
  for a in half:
   for c in half:expected.append([a,c,coeffs[0][32*a+c],coeffs[1][32*a+c]])
  need(expected==r['complete_coarse_coefficient_bank'],'all256 complete coefficient coordinates for every channel/orientation');allcoeff.append(expected)
  need(all(T[1+3*i+j]==U[1+3*i+j]==0 for a,c,T,U in bank for i in range(5) for j in(1,2)),'q3/q5 derivative parity separately')
  formulas=r['battery_formulas'];states=[]
  for bits,rec in enumerate(r['all32_state_controls']):
   gaps=[(T[0]-U[0]+8*sum(((bits>>i)&1)*(T[1+3*i]-U[1+3*i]) for i in range(5)))%16 for a,c,T,U in bank];battery=[gaps[32*z['query'][0]+z['query'][1]] for z in formulas];need(rec['epsilon_bits']==bits and sum(bool(v) for v in gaps)==rec['nonzero_full_queries'] and battery==rec['battery_contrasts'] and any(battery),'complete32-state exact-vector separation');need(hashlib.sha256(bytes(gaps)).hexdigest()==rec['full_contrast_sha256'],'additional retained state digest')
   for z,value in zip(formulas,battery):need(value==(z['base']+8*sum(((bits>>i)&1)*v for i,v in enumerate(z['epsilon_coefficients_mod2'])))%16,'explicit affine battery formulas')
   states.append(gaps)
  need(len(states)==32,'all q2-bit assignments');allstates.append(states)
 # Reconstruct every retained22-name graph/mask and verify its residue via small physical interpolation.
 count=0
 for r in d['physical_controls']['records']:
  rid=next(i for i,Q in enumerate(reps) if r['quotient_rows'] in(Q,op(Q)));side=int(r['reverse']);Q=op(reps[rid]) if side else reps[rid];need(Q==r['quotient_rows'],'physical converse convention');T=[]
  for i in range(5):
   for u in range(4):T.append(sum(1<<(4*i+v) for v in range(u+1,4))|sum(15<<(4*j) for j in range(5) if Q[i]>>j&1))
  need(T==r['old_rows'],'full physical twenty-name TT4 inflation')
  for query,masks,value in zip(r['block_queries'],r['physical_queries'],r['values_mod16']):
   need(masks==[sum(15<<(4*i) for i in range(5) if bit>>i&1) for bit in query],'complete physical block mask');a,c=query;need(value==physical[rid*2048+2*(32*a+c)+side][0],'retained22-name residue independently recovered by small-physical interpolation');count+=1
 need(count==d['physical_controls']['physical_count']==8,'all eight22-name controls')
 return dict(status='PASS',source_sha256=PINS,other_source_pins=d['source_pins'],dependencies=d['dependency_pins'],complete_prime_orbits=groups,literal_small_physical_counts=374784,maximum_literal_order=11,complete_affine_rows=physical,complete_coefficient_banks=allcoeff,complete32_state_vectors=allstates,retained22_name_residues_recovered=count,scope='Independent exact physical interpolation of every local prime5 coarse bank. Full vectors, cubic/quartic coefficients and all formal states retained. No ambient inference here.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
