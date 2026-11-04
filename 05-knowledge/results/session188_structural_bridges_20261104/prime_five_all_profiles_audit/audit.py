from pathlib import Path
from itertools import combinations,permutations,product
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_five_all_profiles'
PINS={'PROOF.md':'7771476a3389d58d05fded981c1b4e4a06bebba981c6d427ca7ccee47f17929a','verify.py':'89ce51a2beb75a61de804f3fa5b43083701e5071ce7c9013a04d8742769679bc','run_bank.cpp':'6835a3af69126d2b59700055255f5d69287fd171d81f0ff7768e4530af74d098','physical.cpp':'86c1954ef7559b6226867706c204af9e19b0d002bbf1f94d92de5a80dcb027ce','normal.json':'f96dae7967cf2c797701da037961b677ec5fa38ce9a9ca27dabc42b2fa6bdd9b','optimized.json':'f96dae7967cf2c797701da037961b677ec5fa38ce9a9ca27dabc42b2fa6bdd9b'}
PAIRS=list(combinations(range(5),2));TYPES='SPEOC'
def need(c,m):
 if not c:raise RuntimeError(m)
def rows(code):
 T=[0]*5
 for b,(i,j) in enumerate(PAIRS):
  u,v=(i,j) if code>>b&1 else(j,i);T[u]|=1<<v
 return T
def classify():
 groups={};ps=list(permutations(range(5)))
 for code in range(1024):
  T=rows(code)
  if any(all(T[v]&M in (0,M) for v in range(5) if not M>>v&1) for M in range(1,31) if 2<=M.bit_count()<5):continue
  canon=min(sum(bool(T[p[i]]>>p[j]&1)<<b for b,(i,j) in enumerate(PAIRS)) for p in ps);groups.setdefault(canon,[]).append(code)
 need(sorted(groups)==[10,41,76] and sorted(map(len,groups.values()))==[24,120,120],'complete prime quotient domain');return groups
def weights(ts):
 u=[dict(S=0,P=2,E=6,O=6,C=2)[t] for t in ts];return [1]+u+[2*int(t in 'OC') for t in ts]+[4*int(u[i]>0 and u[j]>0) for i,j in PAIRS]
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()==h,'dependency '+f)
 groups=classify();reps=[rows(c) for c in sorted(groups)];wire='3\n'+'\n'.join(' '.join(map(str,Q)) for Q in reps)+'\n'
 with tempfile.TemporaryDirectory(prefix='prime5-all-types-independent-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'interpolate.cpp'),'-o',str(exe)],check=True);p=subprocess.run([str(exe)],input=wire.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);bank=[list(map(int,s.split())) for s in p.stdout.decode().splitlines()]
 need(len(bank)==6144 and all(len(v)==21 for v in bank),'complete physical interpolation bank');need(p.stderr.decode().strip()=='physical_counts 129024 max_names 9','physical enumeration accounting')
 summaries=[];allco=[];singleton=[];types=list(product(TYPES,repeat=5))
 for ri,r in enumerate(d['records']):
  need(r['code']==sorted(groups)[ri] and r['rows']==reps[ri] and r['orbit_size']==len(groups[r['code']]),'representative physical convention');values=[]
  for j,(a,c,T,U) in enumerate(r['full_run_bank']):
   need([a,c]==[j//32,j%32] and [T,U]==bank[ri*2048+2*j:ri*2048+2*j+2],'full integer-mod16 channel vectors');values.append([T,U])
  need(len(values)==1024 and all(v%2==0 for T,U in values for v in T[6:11]+U[6:11]),'complete domain and triple-run parity')
  transformed=[]
  for side in (0,1):
   co=[v[side][:] for v in values]
   for bit in range(10):
    for mask in range(1024):
     if mask>>bit&1:
      for j in range(21):co[mask][j]=(co[mask][j]-co[mask^(1<<bit)][j])%16
   need(all(not any(co[k]) for k in range(1024) if (k>>5).bit_count()>2 or (k&31).bit_count()>2),'no hidden high bidegree');transformed.append(co)
  half=[v for v in range(32) if v.bit_count()<=2];cb=[[a,c,transformed[0][32*a+c],transformed[1][32*a+c]] for a in half for c in half]
  need(cb==r['complete_coefficient_bank'],'all complete named coarse cubic/quartic coefficients');allco.append(cb)
  gaps=[[(a-b)%8 for a,b in zip(T,U)] for T,U in values];state_records=[];blind=[]
  for ts,old in zip(types,r['type_classification']):
   w=weights(ts);active=[(j,v) for j,v in enumerate(w) if v];full=[sum(v*row[j] for j,v in active)%8 for row in gaps];first=next((j for j,v in enumerate(full) if v),None);name=''.join(ts)
   rec=dict(types=name,nonzero_queries=sum(bool(v) for v in full),first_separator=None if first is None else[first//32,first%32,full[first]],full_contrast_sha256=hashlib.sha256(bytes(full)).hexdigest());need(rec==old,'each exhaustive type classification and full-vector digest')
   if first is None:blind.append(name)
   # Full vectors are reconstructed losslessly from the complete21-channel bank and these type words.
   state_records.append(rec)
  need(len(state_records)==3125 and blind==([] if ri<2 else['EEEEE'])==r['blind_types'],'complete type coverage');summaries.append(state_records)
  singles=[]
  for M,old in r['singleton_pattern_mod4']:
   need(1<=M<32,'nonempty singleton pattern');g=[(T[0]-U[0]+2*sum(T[1+i]-U[1+i] for i in range(5) if not M>>i&1))%4 for T,U in values];j=next((j for j,v in enumerate(g) if v),None);need(j is not None and old==[j//32,j%32,g[j]] and g[j]==2,'exact singleton-mod4 certificate');singles.append([M,old])
  need([v[0] for v in singles]==list(range(1,32)),'all31 singleton patterns');singleton.append(singles)
 # Rebuild every actual control's labels and masks; independent literal replay is complementary.
 controls=[];can={'S':[0],'P':[2,0],'E':[14,12,8,0],'O':[6,4,0],'C':[2,4,1]}
 for r in d['literal_physical_controls']:
  ri=sorted(groups).index(r['quotient_code']);Q=reps[ri];side=r['side']
  if side:Q=[31^(1<<i)^v for i,v in enumerate(Q)]
  B=[can[t] for t in r['types']];off=[];n=0
  for b in B:off.append(n);n+=len(b)
  T=[]
  for i,b in enumerate(B):T.extend((v<<off[i])|sum(((1<<len(z))-1)<<off[j] for j,z in enumerate(B) if Q[i]>>j&1) for v in b)
  need(T==r['old_rows'],'actual physical graph');masks=[sum(((1<<len(b))-1)<<off[i] for i,b in enumerate(B) if a>>i&1) for a in r['block_query']];need(masks==r['physical_query'],'actual full child masks');unit=pow(3,r['types'].count('C'),8);need(unit==r['odd_normalization_unit'],'odd C3 scalar');a,c=r['block_query'];v=unit*sum(x*y for x,y in zip(bank[ri*2048+2*(32*a+c)+side],weights(r['types'])))%8;need(v==r['predicted_mod8']==r['actual_mod16']%8,'physical-control modular prediction');controls.append([r['quotient_code'],r['types'],side,masks,v])
 need(len(controls)==18,'all physical controls')
 return dict(status='PASS',source_pins=PINS,dependency_pins=d['dependency_pins'],complete_prime_orbits=groups,physical_small_counts=129024,maximum_small_physical_order=9,complete6144_run_vectors=bank,complete_coefficient_banks=allco,complete_type_classifications=summaries,singleton_certificates=singleton,physical_mask_controls=controls,scope='Independent complete small-physical interpolation and exhaustive local classification. Full state vectors are losslessly determined by retained banks, type records and evaluator. Separate physical audit directly recounts the18 large controls; no ambient inference here.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
