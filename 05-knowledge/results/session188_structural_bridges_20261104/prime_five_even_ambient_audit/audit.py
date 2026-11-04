#!/usr/bin/env python3
"""Independent physical-cover interpolation and fixed-ambient audit."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_five_even_ambient'
PINS={'PROOF.md':'1da932d38f9dc57f83defc5ba88cc1736f9805fa827fd258db0ddd3ab03aa05e','verify.py':'cc0917454834af9e6a0c446696e6b965c659d013867d73301c029bfd15831b7e','covers.cpp':'b1667bf213a533271ce813b18309664b13a2adbb2b92b957750ad72ec9db4d02','physical.cpp':'6933dce4c994ed3a02a7461dc5f956157bd05352706ac13ae93db648d1e8da01','physical_controls.py':'44c34f68a2da1ffe87d6f1af6c75186bf1ec35b09b87247ef8f8f68f723b612f','normal.json':'3ca13a3c4945ef2a02ecd478ec38bf0748aff2fbc0cc6dc4273ca04b57f590f9','optimized.json':'3ca13a3c4945ef2a02ecd478ec38bf0748aff2fbc0cc6dc4273ca04b57f590f9'}
def need(c,m):
 if not c:raise RuntimeError(m)
def rows(n,code):
 T=[0]*n
 for b,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>b&1 else(j,i);T[u]|=1<<v
 return T
def slots(T,root,k):
 counts=tuple(k if i==root else 1 for i in range(len(T)))
 @lru_cache(None)
 def f(c,last):
  if not any(c):return 1
  z=0
  for i,v in enumerate(c):
   if v and(last<0 or T[last]>>i&1):
    cc=list(c);cc[i]-=1;z+=f(tuple(cc),i)
  return z
 return f(counts,-1)
def physical_record(Q,query):
 T=[]
 for i in range(5):
  for u in range(4):T.append(sum(1<<(4*i+v) for v in range(u+1,4))|sum(15<<(4*j) for j in range(5) if Q[i]>>j&1)|(1<<20))
 T += [1<<21,(1<<20)-1]
 return T,[sum(15<<(4*i) for i in range(5) if b>>i&1)|(1<<21) for b in query]
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'source pin '+f)
 d=json.loads((SRC/'normal.json').read_text());need(d['status']=='PASS','source status')
 for f,h in d['source_pins'].items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'additional source pin '+f)
 for f,h in d['dependency_pins'].items():need(hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest()==h,'dependency pin '+f)
 with tempfile.TemporaryDirectory(prefix='prime5-cover-independent-') as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'interpolate.cpp'),'-o',str(exe)],check=True)
  p=subprocess.run([str(exe)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);bank=[list(map(int,line.split())) for line in p.stdout.decode().splitlines()]
  need(p.stderr.decode().strip()=='physical_profiles 244 physical_coordinates 732 max_names 11','physical count accounting')
 need(len(bank)==4 and all(len(v)==48 for v in bank),'all independent cover channels')
 records=[]
 for j,query in enumerate([[0,1],[2,18]]):
  T,U=bank[2*j:2*j+2];r=d['ordered_cover_records'][j]
  need(r['query']==query and T==r['source_channels'] and U==r['target_channels'],'every complete affine profile vector')
  base=[(T[k]-U[k])%16 for k in range(3)];derivs=[[(T[3+3*i+k]-U[3+3*i+k])%2 for k in range(3)] for i in range(15)]
  need(base==[8*j,0,0] and derivs==[[int(i in (6,9)),0,0] for i in range(15)],'exact coefficient identity, all q3/q5 channels included')
  need(base==r['base_contrast'] and derivs==r['perturbation_contrasts_mod2'],'retained coefficient matches');records.append((base,derivs))
 # Retain every complete six-coordinate contrast, not merely its digest.
 states=[];hist={}
 for bits in range(32768):
  state=[]
  for base,ds in records:state.extend((base[h]+8*sum(((bits>>i)&1)*ds[i][h] for i in range(15)))%16 for h in range(3))
  need(state in ([0,0,0,8,0,0],[8,0,0,0,0,0]),'complete normalized profile domain')
  hist[tuple(state)]=hist.get(tuple(state),0)+1;states.append(state)
 need(sorted(hist.values())==[16384,16384] and d['formal_states']==32768,'full domain accounting')
 need(hashlib.sha256(bytes(v for s in states for v in s)).hexdigest()==d['complete_state_vector_sha256'],'supplementary state digest')
 # Independently count every collapsed root-run coefficient in small fixed ambients.
 contexts=[]
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   R=rows(n,code)
   for root in range(n):
    b=[slots(R,root,k) for k in range(1,n+1)]
    need(b[0]%2==1 and (8*b[0])%16==8,'odd first ambient coefficient')
    contexts.append([R,root,b])
 need(len(contexts)==285,'all rooted ambient contexts through four names')
 R=[2,4,1];coeff=[slots(R,0,k) for k in range(1,4)];need(coeff==[3,1,0],'exact cyclic exterior run coefficients')
 Q=[24,17,3,6,12];physical=[]
 for side,r in enumerate(d['physical_controls']['records']):
  O=[31^(1<<i)^Q[i] for i in range(5)] if side else Q;T,masks=physical_record(O,[2,18])
  need(r['reverse']==bool(side) and r['old_rows']==T and r['physical_query']==masks and r['outer_quotient']==R and r['root']==0 and r['block_query']==[2,18],'exact24-name physical graph and lifted incidence')
  # The two outside names are uniform to the entire augmented core.
  for mask in masks:need(not(mask>>20&1) and(mask>>21&1),'ports mimic root, no global exterior complement')
  value=sum(coeff[h]*bank[2+side][h] for h in range(3))%16
  need(value==r['value_mod16'],'independent small-physical interpolation plus exact exterior expansion')
  physical.append(dict(old_rows=T,physical_query=masks,value_mod16=value))
 need([r['value_mod16'] for r in physical]==[13,5],'both actual ambient residues')
 return dict(status='PASS',source_pins=PINS,other_source_pins=d['source_pins'],dependency_pins=d['dependency_pins'],physical_profiles=244,physical_profile_coordinates=732,maximum_literal_physical_order=11,complete_affine_cover_vectors=bank,complete32768_profile_contrasts=states,rooted_ambient_controls=contexts,actual24_name_controls_recovered=physical,scope='Special complete-profile modulo16 transfer for regular prime5 and inherited even modulo8 transfer for the other two prime quotients. Physical data recovered by independent small-physical ordered-cover interpolation, not a second direct24-name count. Source-relative original families only; fixed exterior.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
