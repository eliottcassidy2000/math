#!/usr/bin/env python3
"""Exact short obstruction certificates. Completeness replays the sibling onset solver."""
import argparse,hashlib,heapq,importlib.util,json,math,sys
sys.dont_write_bytecode=True
from pathlib import Path
MASKS=(1187208845135,1453217412461,1487448155501)
PATH43={
MASKS[0]:[1,24,25,39,10,6,43,21,15,34,30,19,17,32,4,5,11,38,26,23,2,14,35,29,20,16,9,40,41,8,28,36,13,12,37,27,22,42,7,18,31,33,3],
MASKS[1]:[1,35,29,20,16,33,31,18,7,42,39,10,15,21,43,38,26,23,13,36,28,8,41,40,9,27,37,12,24,25,11,5,4,32,17,19,6,30,34,2,14,22,3]}
def need(ok,msg):
 if not ok:raise ValueError(msg)
def phase(mask,v):return (mask>>(v-1))&1 if v<=41 else v%2
def adj(mask,n):
 return {v:[u for u in range(1,n+1) if u!=v and math.isqrt(u+v)**2==u+v and phase(mask,u)!=phase(mask,v)] for v in range(1,n+1)}
def valid_path(p,mask):
 n=len(p);need(set(p)==set(range(1,n+1)),'coverage');need(p[0]==1 and p[-1]==(3 if n%2 else 8),'ports')
 need(all(math.isqrt(u+v)**2==u+v for u,v in zip(p,p[1:])),'square edges')
 need(all(phase(mask,v)==(j%2) for j,v in enumerate(p,1)),'fixed phase')
def minimum_component_walk(end):
 # A=0={0,1}, B=1={2,7}, C=2={3,6}, D=3={4,5}.
 # Each intercomponent transition requires a nonempty defect run;
 # crossing AD/BC needs at least two defects.
 group=[0,1,1,0];heap=[(0,0,1)];best={(0,1):0}
 while heap:
  cost,v,seen=heapq.heappop(heap)
  if cost!=best[v,seen]:continue
  if v==end and seen==15:return cost
  for w in range(4):
   if w==v:continue
   new=cost+1+int(group[v]!=group[w]);key=(w,seen|(1<<w))
   if new<best.get(key,10**9):best[key]=new;heapq.heappush(heap,(new,*key))
def run(origin_dir):
 onset_checker=Path(origin_dir)/'check_independent.py'
 need(onset_checker.is_file(),'missing sibling origin-pair checker; supply --origin-dir')
 spec=importlib.util.spec_from_file_location('onset_independent',onset_checker);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 reference=mod.read_reference()
 reference_data=json.loads(reference)
 a=mod.enumerate_by_degrees(41);b=mod.enumerate_by_degrees(42)
 common=sorted(set(a['masks'])&{x&((1<<41)-1) for x in b['masks']})
 need(common==list(MASKS),'complete three-mask catalogue')
 saved=reference_data['nice_pairs'][-1]
 need(saved['n']==41 and sorted(w['odd_mask_common'] for w in saved['witnesses'])==common,'compressed reference catalogue')
 for m,p in PATH43.items():valid_path(p,m)
 # First phase: four degree-two vertices saturate both alternative neighbors of29.
 A=adj(MASKS[0],44)
 required={18:[7,31],42:[7,22],16:[9,20],44:[20,37],29:[7,20,35]}
 need(all(A[v]==ns for v,ns in required.items()),'first obstruction neighborhoods')
 need(all(v not in (1,8) for v in [18,42,16,44,29]),'obstruction vertices internal')
 # Other phases: newly added internal vertex has only one possible incident edge.
 B=adj(MASKS[1],44);C=adj(MASKS[2],43)
 need(B[44]==[37] and C[43]==[38],'degree-one obstructions')
 components=[{0,1},{2,7},{3,6},{4,5}]
 for r in range(8):
  need(any(r in c and (1-r)%8 in c for c in components),'odd-square residue components')
 attachment={}
 for r in range(8):
  residues={(-r)%8,(4-r)%8};grp='AD' if r%4 in (0,3) else 'BC';allowed=components[0]|components[3] if grp=='AD' else components[1]|components[2]
  need(residues<=allowed,'single-defect attachment superpair');attachment[r]=dict(residues=sorted(residues),superpair=grp)
 need(minimum_component_walk(0)==6 and minimum_component_walk(2)==4,'weighted component walk lower bounds')
 defects=[]
 for m in MASKS:
  D=[v for v in range(1,42) if phase(m,v)!=v%2];need(len([v for v in D if v%2])==len([v for v in D if not v%2]),'balanced defects');need(not(set(D)&{1,3,8}),'port defects')
  defects.append(dict(mask=m,defects=D,odd_defects=len([v for v in D if v%2]),even_defects=len([v for v in D if not v%2])))
 hostile={}
 for label,p,mask in [('flipped_mask',PATH43[MASKS[0]],MASKS[0]^2),('wrong_port',PATH43[MASKS[0]][::-1],MASKS[0])]:
  try:valid_path(p,mask);hostile[label]=False
  except ValueError:hostile[label]=True
 need(all(hostile.values()),'hostile controls')
 return dict(status='PASS',completeness_dependency='origin_pairs/check_independent.py',decoded_origin_certificate_sha256=hashlib.sha256(reference).hexdigest(),onset_path_counts=[a['path_count'],b['path_count']],common_phase_masks=common,positive_fixed_phase_43_paths={str(k):v for k,v in PATH43.items()},first_failures=[dict(mask=MASKS[0],n=44,certificate=required),dict(mask=MASKS[1],n=44,neighbor_of_new_vertex=B[44]),dict(mask=MASKS[2],n=43,neighbor_of_new_vertex=C[43])],no_coherent_original_port_family_41_through_44=True,defect_data=defects,attachment_residues=attachment,minimum_weighted_walk_defects=dict(even_A_to_A=6,odd_A_to_C=4),hostile_controls=hostile)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--origin-dir',default=str(Path(__file__).resolve().parent.parent/'session125_origin32_numeration_20261031'/'origin_pairs'),help='directory containing the published origin-pair checker and compressed certificate');ap.add_argument('--output');a=ap.parse_args();s=json.dumps(run(a.origin_dir),indent=2)+'\n'
 if a.output:Path(a.output).write_text(s)
 print(s)
