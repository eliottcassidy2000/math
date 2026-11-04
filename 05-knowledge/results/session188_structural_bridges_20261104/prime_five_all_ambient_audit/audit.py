from pathlib import Path
from itertools import combinations,product
from math import prod
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_five_all_ambient'
PINS={'PROOF.md':'1dad821eff721e77ede69a3b0b8a1b62d0a3b0dce7395b883e89ddfa5ad9cea4','verify.py':'ec86992a83ee3a0fdf8aca4f7122139e1d514e4d261f9e9b6a71300aacd29fff','covers.cpp':'434cf002661419b1335658d46e6a8dc3bd7f5c28920be4e367d4e9d55557aca7','physical.cpp':'86c1954ef7559b6226867706c204af9e19b0d002bbf1f94d92de5a80dcb027ce','normal.json':'cdec4c53587ef40f5e33c31b83ee5b83a8522cf082ed6dfb6a068f4e3808fbea','optimized.json':'cdec4c53587ef40f5e33c31b83ee5b83a8522cf082ed6dfb6a068f4e3808fbea'}
REPS=[[20,1,2,7,14],[24,17,3,6,12]];TYPE='PEOC';PROFILES=[[0,1,2,0],[0,1,6,4],[0,1,6,6],[0,1,2,2]]
def need(c,m):
 if not c:raise RuntimeError(m)
def sub(M):
 s=M
 while True:
  yield s
  if not s:break
  s=(s-1)&M
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','status')
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()==h,'dependency '+f)
 shapes=[[1]*5]
 for i in range(5):
  for k in (2,3):s=[1]*5;s[i]=k;shapes.append(s)
 for i,j in combinations(range(5),2):
  for k,l in product((2,3),repeat=2):s=[1]*5;s[i]=k;s[j]=l;shapes.append(s)
 need(shapes==d['run_shape_order'] and len(shapes)==51 and len(set(map(tuple,shapes)))==51,'exhaustive relevant run multiplicities')
 wire='2\n'+'\n'.join(' '.join(map(str,Q)) for Q in REPS)+'\n'
 with tempfile.TemporaryDirectory(prefix='prime5-all-ambient-independent-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'interpolate.cpp'),'-o',str(exe)],check=True);p=subprocess.run([str(exe)],input=wire.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);bank=[list(map(int,s.split())) for s in p.stdout.decode().splitlines()]
 need(len(bank)==1024 and all(len(v)==153 for v in bank),'complete independent cover kernel');need(p.stderr.decode().strip()=='physical_profiles 52224 max_names 11','physical computation accounting')
 half=[a for a in range(32) if a.bit_count()<=2];queries=list(product(half,repeat=2));states=[];coefficients=[];look={}
 for ri,r in enumerate(d['records']):
  need(r['quotient_rows']==REPS[ri] and r['quotient_code']==[10,76][ri],'physical representative coordinate system');values={}
  for j,(a,c,T,U) in enumerate(r['complete_kernel_bank']):
   need((a,c)==queries[j] and [T,U]==bank[ri*512+2*j:ri*512+2*j+2],'every full153-channel basis vector');values[a,c]=[T,U];look[r['quotient_code'],a,c]=[T,U]
  need(len(values)==256,'complete Boolean basis')
  # All bidegree(<=2,<=2) coordinates, including cubic/quartic ones and zeros.
  co=[]
  for a,c in queries:
   coord=[]
   for side in (0,1):
    v=[0]*153
    for x,y in product(sub(a),sub(c)):
     sign=(-1)**((a^x).bit_count()+(c^y).bit_count())
     for j,z in enumerate(values[x,y][side]):v[j]=(v[j]+sign*z)%8
    coord.append(v)
   co.append([a,c,*coord])
  for a,c in queries:
   for side in (0,1):
    got=[0]*153
    for x,y,T,U in co:
     if x&a==x and y&c==y:
      for j,z in enumerate((T,U)[side]):got[j]=(got[j]+z)%8
    need(got==values[a,c][side],'full coefficient/basis reconstruction')
  coefficients.append(co);rr=[];given=iter(r['all512_odd_type_certificates'])
  for ids in product(range(4),repeat=5):
   if sum(i>=2 for i in ids)%2==0:continue
   w=[prod(PROFILES[ids[i]][s[i]] for i in range(5))%8 for s in shapes];active=[(j,z) for j,z in enumerate(w) if z];full=[];first=None
   for a,c in queries:
    T,U=values[a,c];v=[sum(z*(T[3*j+h]-U[3*j+h]) for j,z in active)%8 for h in range(3)];full.extend(v)
    if first is None and v==[4,0,0]:first=[a,c]
   name=''.join(TYPE[i] for i in ids);old=next(given);need(first is not None and old==dict(types=name,query=first,all256_profile_contrasts_sha256=hashlib.sha256(bytes(full)).hexdigest()),'exhaustive odd domain and selected complete-profile separator');rr.append(dict(types=name,query=first,complete256_profile_contrasts_hex=bytes(full).hex()))
  need(len(rr)==512 and next(given,None) is None,'all512 assignments exactly once');states.append(rr)
 # Re-evaluate the singleton and Q41 branches directly from the pinned local bank.
 localpath=ROOT.parent/'prime_five_all_profiles/normal.json';need(hashlib.sha256(localpath.read_bytes()).hexdigest()=='f96dae7967cf2c797701da037961b677ec5fa38ce9a9ca27dabc42b2fa6bdd9b','local complete source bank pin');local=json.loads(localpath.read_text());singleton=[]
 for r in local['records']:
  rows=[]
  for M,query in r['singleton_pattern_mod4']:
   a,c,gap=query;z=r['full_run_bank'][32*a+c];T,U=z[2:];v=(T[0]-U[0]+2*sum(T[1+i]-U[1+i] for i in range(5) if not M>>i&1))%4;need(z[:2]==[a,c] and v==gap==2,'singleton baseline full-response mod4');rows.append([M,query])
  need([z[0] for z in rows]==list(range(1,32)),'all singleton masks');singleton.append(dict(code=r['code'],records=rows))
 need(singleton==d['singleton_pattern_mod4_certificates'],'retained singleton certificates')
 middle=local['records'][1]['full_run_bank'][1];T,U=middle[2:];need(middle[:2]==[0,1] and (T[0]-U[0]+2*sum(T[1+i]-U[1+i] for i in range(5)))%4==2,'Q41 nontrivial full mod4 branch')
 physical=[];can=[[2,0],[14,12,8,0],[6,4,0],[2,4,1]]
 for r in d['physical_ambient_controls']:
  Q=REPS[[10,76].index(r['quotient_code'])];side=r['side']
  if side:Q=[31^(1<<i)^v for i,v in enumerate(Q)]
  ids=[TYPE.index(t) for t in r['types']];B=[can[i] for i in ids];off=[];n=0
  for b in B:off.append(n);n+=len(b)
  T=[]
  for i,b in enumerate(B):T.extend((v<<off[i])|sum(((1<<len(z))-1)<<off[j] for j,z in enumerate(B) if Q[i]>>j&1)|(1<<n) for v in b)
  T += [1<<(n+1),(1<<n)-1];a,c=r['block_query'];masks=[sum(((1<<len(b))-1)<<off[i] for i,b in enumerate(B) if bit>>i&1)|(1<<(n+1)) for bit in (a,c)];need(T==r['old_rows'] and masks==r['query'] and r['outer_rows']==[2,4,1],'physical ambient rows and exact root-incidence masks')
  unit=pow(3,ids.count(3),8);w=[prod(PROFILES[ids[i]][s[i]] for i in range(5))%8 for s in shapes];K=look[r['quotient_code'],a,c][side];v=[unit*sum(z*K[3*j+h] for j,z in enumerate(w))%8 for h in range(3)];need(v==r['core_first_three_mod8'] and unit==r['odd_normalization_unit'],'actual odd scalar/core cover profile');ambient=(3*v[0]+v[1])%8;need(ambient==r['predicted_ambient_mod8']==r['actual_ambient_mod16']%8,'exact C3 exterior expansion prediction');physical.append([r['quotient_code'],r['types'],side,masks,v,ambient])
 need(len(physical)==4 and all((physical[i][-1]-physical[i+1][-1])%8==4 for i in (0,2)),'all actual fixed-ambient controls')
 return dict(status='PASS',source_pins=PINS,dependency_pins=d['dependency_pins'],small_physical_profiles=52224,physical_coordinates=156672,maximum_small_order=11,complete1024_kernel_vectors=bank,all_shape_coordinates=shapes,complete_coefficient_banks=coefficients,all_odd_state_full_contrasts=states,singleton_certificates=singleton,physical_ambient_mask_controls=physical,scope='Independent exact physical ordered-cover interpolation of all basis kernels and full odd-state contrasts, plus exhaustive branch review. Direct large physical recounts are in the complementary prime_five_general_physical_audit. Fixed ambient and source-relative family scope only.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
