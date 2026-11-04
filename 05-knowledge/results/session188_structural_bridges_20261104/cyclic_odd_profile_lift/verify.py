"""Exact odd-total augmented-cover residues and fixed-ambient physical lift."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SOURCE=ROOT.parent/'five_profile_types'
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
@lru_cache(None)
def profile(T):return tuple(u.profile(list(T)))
def phys(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def build(blocks,query,outer,targetblocks=None):
 T,offs=u.blow([2,4,1],blocks);U,_=u.blow([4,1,2],blocks if targetblocks is None else targetblocks);n=len(T);a,c=[phys(mask,blocks,offs) for mask in query];outside=sum(1<<(n+j-1) for j in range(1,len(outer)) if outer[j]&1);T,_=u.blow(outer,[T]+[[0] for _ in outer[1:]]);U,_=u.blow(outer,[U]+[[0] for _ in outer[1:]]);return T,U,[(a|outside,c|outside)]
def run(exe,cases):
 lines=[str(len(cases))]
 for T,U,qs in cases:lines += [f'{len(T)} {len(qs)}',' '.join(map(str,T)),' '.join(map(str,U))]+[f'{a} {c}' for a,c in qs]
 raw=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()));out=[];at=0
 for T,U,qs in cases:out.append([raw[at+4*j:at+4*j+4] for j in range(len(qs))]);at+=4*len(qs)
 need(at==len(raw),'consume all four-coordinate outputs');return out
def modules(T):return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and all(T[v]&M in(0,M) for v in range(len(T)) if not M>>v&1)]
def flip(T,selected):
 R=list(T)
 for M in selected:
  for i in range(len(T)):
   if M>>i&1:R[i]^=M^(1<<i)
 return R

def main():
 source=json.loads((SOURCE/'normal.json').read_text());classes=[r for r in source['canonical_records'] if 'S' not in r['types'] and len(r['source_old_rows'])%2];need(len(classes)==32,'complete odd-total nontrivial type domain');queries=list(product(range(8),repeat=2));records=[];cases=[(r['source_old_rows'],r['target_old_rows'],r['physical_queries']) for r in classes]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/'literal.cpp'),'-o',str(exe)],check=True);values=run(exe,cases);lookup={}
  for r,vv in zip(classes,values):
   need([v[:2] for v in vv]==r['values'],'all frozen canonical Hamilton values reproduced');good=[j for j,(a,b,c,d) in enumerate(vv) if(a-b)%8==4 and(c-d)%8==0];need(good,'every odd-total class has higher-coordinate-zero separator');item={'types':r['types'],'child_rows':r['child_rows'],'source_old_rows':r['source_old_rows'],'target_old_rows':r['target_old_rows'],'physical_queries':r['physical_queries'],'columns':['H_source','H_target','q2_source','q2_target'],'values':vv,'good_query_indices':good,'chosen_block_query':list(queries[good[0]])};records.append(item);lookup[tuple(r['types'])]=item
  need(min(len(r['good_query_indices']) for r in records)==4 and max(len(r['good_query_indices']) for r in records)==20,'complete good-query range')
  extra=[];meta=[]
  odd5=next(r['rows'] for r in source['actual_profile_controls'] if len(r['rows'])==5 and r['odd_Hamilton_cycles']%2 and r['q'][1]%8!=3)
  variants=[(['E','E','O'],[u.tt(4),u.tt(4),u.tt(5)]),(['O','O','O'],[u.tt(5),u.tt(3),u.tt(3)]),(['C','P','E'],[odd5,u.tt(2),[8,1,3,6]]),(['C','C','C'],[odd5,[2,4,1],[2,4,1]])]
  canmap=dict(zip(source['canonical_type_order'],source['canonical_profiles']))
  for types,blocks in variants:
   item=lookup[tuple(types)];query=item['chosen_block_query'];j=queries.index(tuple(query));unit=1
   for typ,B in zip(types,blocks):unit=unit*profile(tuple(B))[1]*pow(canmap[typ][1],-1,8)%8
   outers=[[0],[2,4,1]]+([[6,12,24,17,3]] if sum(map(len,blocks))<=11 else [])
   for outer in outers:
    case=build(blocks,query,outer);idx=len(extra);extra.append(case);meta.append({'types':types,'child_rows':blocks,'source_profiles':[profile(tuple(B)) for B in blocks],'unit':unit,'block_query':query,'canonical_core_values':item['values'][j],'outer_rows':outer,'source_old_rows':case[0],'target_old_rows':case[1],'physical_query':case[2][0],'case':idx})
  # Preserve an unsafe local query and the safe selected query in the same odd type.
  hostile=[];blocks=[u.tt(2),u.tt(2),u.tt(3)];safe=lookup['P','P','O']['chosen_block_query']
  for query in [[1,2],safe]:
   for outer in [[0],[2,4,1]]:
    case=build(blocks,query,outer);idx=len(extra);extra.append(case);hostile.append({'safe':query==safe,'block_query':query,'outer_rows':outer,'physical_query':case[2][0],'case':idx})
  # Every original child-module family in the unequal odd-total PPO source.
  fams=[];mods=[modules(B) for B in blocks];sourceprofiles=[profile(tuple(B)) for B in blocks];dual=[7^safe[1],7^safe[0]]
  for outer in [[0],[2,4,1]]:
   baseline=[]
   for query in [safe,dual]:baseline.append(len(extra));extra.append(build(blocks,query,outer))
   for masks in product(*[range(1<<len(mm)) for mm in mods]):
    target=[];selected=[]
    for B,mm,m in zip(blocks,mods,masks):
     chosen=[s for k,s in enumerate(mm) if m>>k&1];selected.append(chosen);target.append(u.op(flip(B,chosen)))
    targetprofiles=[profile(tuple(B)) for B in target];need(all(b>=a for qa,qb in zip(sourceprofiles,targetprofiles) for a,b in zip(qa,qb)),'root-family child monotonicity');ids=[]
    for query in[safe,dual]:ids.append(len(extra));extra.append(build(blocks,query,outer,target))
    fams.append({'outer_rows':outer,'selected_original_child_modules':selected,'target_child_rows':target,'target_profiles':targetprofiles,'baseline_cases':baseline,'target_cases':ids})
  more=run(exe,extra)
  for r in meta:
   r['values']=more[r.pop('case')][0];a,b,c,d=r['values'];need((a-b)%8==4,'all actual larger-child and ambient separators');
   if len(r['outer_rows'])==1:need((c-d)%8==0 and all((x-r['unit']*y)%8==0 for x,y in zip(r['values'],r['canonical_core_values'])),'actual all-cover unit transfer')
  for r in hostile:r['values']=more[r.pop('case')][0];r['contrast']=[r['values'][0]-r['values'][1],r['values'][2]-r['values'][3]]
  need([r['contrast'][0] for r in hostile[:2]]==[-68,-600] and hostile[0]['contrast'][1]%8==4,'unsafe odd-total query loses its residue');need(all(r['contrast'][0]%8==4 for r in hostile[2:]),'selected query repairs the same type')
  strict=0
  for r in fams:
   base=[more[k][0] for k in r.pop('baseline_cases')];target=[more[k][0] for k in r.pop('target_cases')];d=[v[0]-v[1] for v in base];need(d[0]==-d[1] and d[0]%8==4,'dual baseline direction');chosen=0 if d[0]<0 else 1;need(target[chosen][0]-target[chosen][1]<=d[chosen]<0,'fixed dual member detects every unequal-profile root family');r['baseline_values']=base;r['target_values']=target;r['negative_query_index']=chosen;isstrict=r['target_profiles']!=sourceprofiles;r['profile_strict']=isstrict;strict+=isstrict
  physical_counts=2*sum(len(z) for z in values)+2*sum(len(z) for z in more)
 deps=['five_profile_types/PROOF.md','five_profile_types/verify.py','five_profile_types/normal.json','mod8_transfer/PROOF.md','monotone_dual_queries/PROOF.md','uniform_block_starts/verify.py']
 return {'status':'PASS','odd_total_type_triples':32,'queries_per_triple':64,'canonical_q1_q2_records':records,'larger_child_ambient_controls':meta,'unsafe_and_safe_odd_query_controls':hostile,'unequal_profile_root_family_contexts':len(fams),'profile_strict_family_contexts':strict,'root_family_records':fams,'physical_H_counts':physical_counts,'physical_q2_counts':physical_counts,'source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','literal.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in deps},'scope':'All nontrivial cyclic children except all-even-order>=4 have a higher-cover-zero modulo8 coarse separator in arbitrary fixed ambient, and a dual two-query original-root-family battery. No exact collision or arbitrary-prime quotient theorem.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
