#!/usr/bin/env python3
"""Complete exact C3 mod16 table by an independent finite run-word automaton."""
from pathlib import Path
from itertools import product
import hashlib,json
from automaton import kernel
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'even_cyclic_mod16'
PINS={'PROOF.md':'d759e62510d0469dd9e864a46117f2bccc02e8256892162664a021821f07a78c','verify.py':'b64c97ad9784657af1ed8478a957deb71b129418d5507cb16d526a6bae8ebee0','literal.cpp':'ab6bc0b847283f2031b3d46c1b9690731f1497c02159152476fae8fbba27b8ec','normal.json':'374b0fe3faf612a7fd99bb8b903741d88b5a7e9155370d8cb1b5e2f8e362be70','optimized.json':'374b0fe3faf612a7fd99bb8b903741d88b5a7e9155370d8cb1b5e2f8e362be70'}
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def main():
 for f,h in PINS.items():need(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h,'frozen source '+f)
 data=json.loads((SRC/'normal.json').read_text());need(data['status']=='PASS','producer status')
 for f,h in data['dependency_pins'].items():need(hashlib.sha256((SRC.parent/f).read_bytes()).hexdigest()==h,'dependency '+f)
 types=[(0,1,b,c,8,d) for b,c,d in product((6,14),(4,12),(0,8))];keys=list(product(range(8),repeat=2));need(list(map(tuple,data['formal_profile_upper_domain']))==types and list(map(tuple,data['block_query_order']))==keys,'entire exact formal domain and physical query order')
 K=[(kernel(a,c,False),kernel(a,c,True)) for a,c in keys];nums=[]
 for q,(F,G),rec in zip(keys,K,data['complete_numerator_polynomials']):
  need(tuple(rec[0])==q,'numerator physical order');D={}
  for e in product(range(6),repeat=3):
   z=0
   for k,c in enumerate((1,-3,3,-1)):
    if min(e)>=k:
     t=tuple(x-k for x in e);z+=c*(F.get(t,0)-G.get(t,0))
   if z:D[e]=z
  need(D=={tuple(z[:3]):z[3] for z in rec[1]},'complete numerator and zero exterior coordinates by run automaton');nums.append(D)
 need(len(data['complete_formal_records'])==512,'complete formal table');hist={};common=set(range(64));exact_records=[];residue_vectors=[]
 for ids,rec in zip(product(range(8),repeat=3),data['complete_formal_records']):
  need(tuple(rec['profile_type_indices'])==ids,'ordered type triple coverage');P=[types[i] for i in ids];pairs=[]
  for kernels in K:pairs.append([sum(c*P[0][e[0]]*P[1][e[1]]*P[2][e[2]] for e,c in d.items()) for d in kernels])
  gaps=[a-b for a,b in pairs];res=[z%16 for z in gaps];need(gaps==rec['integer_functionals'] and res==rec['residues_mod16'],'all64 full exact integer contrasts, not hash comparison');nz=[i for i,z in enumerate(res) if z];need(nz and set(res)<=set((0,8)) and tuple(rec['first_separator'])==keys[nz[0]],'noncollision and first-query convention');need(res[keys.index((1,3))]==8,'fixed actual query on every upper-domain triple');hist[len(nz)]=hist.get(len(nz),0)+1;common&=set(nz);exact_records.append(dict(type_indices=ids,full64_exact_source_target_values=pairs,full64_exact_contrasts=gaps));residue_vectors.append(tuple(res))
 need(hist=={12:64,24:256,28:192} and hist=={int(k):v for k,v in data['nonzero_query_histogram'].items()},'complete residue histogram');universal=[keys[i] for i in sorted(common)];need(universal==[(1,3),(1,5),(2,3),(2,6),(4,5),(4,6)]==list(map(tuple,data['universal_queries'])),'all universal coarse queries')
 # Explicit diagnostic: exactly the three q2 bits affect these C3 contrast vectors.
 lookup={tuple(r['type_indices']):residue_vectors[i] for i,r in enumerate(exact_records)}
 for ids in product(range(8),repeat=3):need(lookup[ids]==lookup[tuple(i&4 for i in ids)],'q3/q5 bits absent from the full C3 contrast vector')
 return dict(status='PASS',source_sha256=PINS,dependency_sha256=data['dependency_pins'],formal_profile_upper_domain=types,queries=keys,formal_triples=512,exact_contrast_coordinates=32768,exact_response_coordinates=65536,distinct_residue_vectors=len(set(residue_vectors)),nonzero_query_histogram=hist,universal_queries=universal,independent_full_numerators=[[[*e,c] for e,c in sorted(p.items())] for p in nums],records=exact_records,scope='Independent complete formal run-table proof and local mathematical review. Actual-profile/literal/named-bank controls are independently reviewed in the separate circuit-worker audit. No ambient assertion follows solely from this local table.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
