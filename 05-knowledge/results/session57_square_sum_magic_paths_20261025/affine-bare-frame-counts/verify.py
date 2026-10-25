"""All bare unique-sum frames over a prescribed affine F3-plane."""
import pathlib,json,itertools as it,collections
from field import Field
P=pathlib.Path(__file__).parent

def require(test,message):
 if not test:raise RuntimeError(message)

def key(B,D):return tuple(sorted((tuple(sorted(B)),tuple(sorted(D)))))

def check(r):
 q=r['q'];F=Field(len(r['modulus_low_first'])-1,r['modulus_low_first']);add=F.add;sub=F.sub
 A=set(r['square_plane']);C=r['C'];U={sub(a,C) for a in A};E=range(q)
 require(A<={F.mul(x,x) for x in E},'square-valued witness')
 dirs=sorted({tuple(sorted((0,u,F.neg(u)))) for u in U-{0}})
 isap=lambda B:add(add(B[0],B[1]),B[2])==0
 frames=set()
 for L in dirs:
  cosets={tuple(sorted(add(b,l) for l in L)) for b in E}
  require(len(cosets)==q//3,'line coset count')
  for B in cosets:
   target={sub(a,B[0]) for a in A}
   target_cosets=sorted({tuple(sorted(add(t,l) for l in L)) for t in target})
   require(len(target_cosets)==3,'three target cosets')
   for D in it.product(*target_cosets):
    require({add(b,d) for b in B for d in D}==A,'transversal sums')
    frames.add(key(B,D))
 require(len(frames)==30*q,'all overlapping allowed frames')
 all_types=collections.Counter(sum(isap(T) for T in k) for k in frames)
 require(all_types=={1:24*q,2:6*q},'all type split')
 good={k for k in frames if not set(k[0])&set(k[1])}
 types=collections.Counter(sum(isap(T) for T in k) for k in good)
 require(types=={1:24*(q-9),2:6*(q-9)},'disjoint type split')
 for B,D in good:
  require(len({sub(b,d) for b in B for d in D})==9,'nine distinct differences')
 brute=None
 if q==81:
  tab=[[add(a,b) for b in E] for a in E]
  neg=[F.neg(a) for a in E]
  candidate=[sum(1<<tab[a][neg[b]] for a in A) for b in E]
  full=set();checked=0
  # Enumerate EVERY ambient triple B, then every triple D satisfying all nine
  # cross-sum memberships; unique-sum test is separate, not assumed.
  for B in it.combinations(E,3):
   mask=candidate[B[0]]&candidate[B[1]]&candidate[B[2]]
   if mask.bit_count()<3:continue
   possible=[d for d in E if mask>>d&1]
   for D in it.combinations(possible,3):
    checked+=1
    if len({tab[b][d] for b in B for d in D})==9:full.add(key(B,D))
  require(full==frames,'independent exhaustive all triples census')
  # Gauge action on unordered pairs: explicitly verify thirty free orbits,
  # each having nine nondisjoint and q-9 disjoint members.
  unseen=set(full);orbits=[]
  while unseen:
   B,D=next(iter(unseen))
   orbit={key([tab[b][g] for b in B],[tab[d][neg[g]] for d in D]) for g in E}
   require(len(orbit)==q and orbit<=unseen,'free disjoint gauge orbits')
   require(sum(k in good for k in orbit)==q-9,'nine bad per orbit')
   unseen-=orbit;orbits.append(orbit)
  require(len(orbits)==30,'thirty gauge orbits')
  paths=set();ap_paths=set()
  for B,D in good:
   before=len(paths)
   for b in it.permutations(B):
    for d in it.permutations(D):
     path=tuple(x for pair in zip(b,d) for x in pair);path=min(path,path[::-1])
     recovered=key(path[::2],path[1::2]);require(recovered==key(B,D),'path recovers its frame')
     paths.add(path)
     if isap(B) and isap(D):ap_paths.add(path)
   require(len(paths)==before+36,'no path collisions across frames')
  require(len(paths)==1080*(q-9) and len(ap_paths)==216*(q-9),'unique path counts')
  brute={'ambient_triples_B_examined':q*(q-1)*(q-2)//6,'candidate_ordered_triple_pairs_examined':checked,'all_frames':len(full),'gauge_orbits':len(orbits),'orbit_size':q,'bad_per_orbit':9,'unique_six_vertex_paths':len(paths),'unique_AP_AP_six_vertex_paths':len(ap_paths)}
 return {'q':q,'zero_in_plane':0 in A,'all_frames_before_disjointness':len(frames),'all_disjoint_frames':len(good),'disjoint_AP_AP_frames':types[2],'disjoint_AP_nonAP_frames':types[1],'distinct_unoriented_six_vertex_paths':36*len(good),'distinct_AP_AP_paths':36*types[2],'distinct_unordered_two_path_covers':36*len(good),'brute':brute}
local_paths=[]
for b in it.permutations(range(3)):
 for d in it.permutations(range(3,6)):
  path=tuple(x for pair in zip(b,d) for x in pair)
  mask=sum(1<<(a*3+(z-3)) for a,z in (tuple(sorted(e)) for e in zip(path,path[1:])))
  local_paths.append(mask)
require(len(set(local_paths))==36,'local path masks')
cover_pairs=sum((a|b)==511 for a,b in it.combinations(local_paths,2))
require(cover_pairs==36,'two Hamilton-path covers per skeleton')
rows=[]
for r in json.loads((P/'sources/ap-frame-results.json').read_text())['rows']:
 rows.append(check(r));print('PASS',rows[-1],flush=True)
(P/'results.json').write_text(json.dumps({'status':'PASS','local_Hamilton_paths':36,'local_unordered_two_path_covers':cover_pairs,'rows':rows},indent=2)+'\n')
