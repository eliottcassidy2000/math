"""Exact cycle-parity collapse and coefficient identities; no separation claim."""
from pathlib import Path
from itertools import product
import hashlib,importlib.util,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('u',ROOT.parent/'uniform_block_starts'/'verify.py');u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def phys(mask,blocks,offs):return sum(((1<<len(B))-1)<<o for i,(B,o) in enumerate(zip(blocks,offs)) if mask>>i&1)
def main():
 h=json.loads((ROOT.parent/'coarse_mod8_affine'/'normal.json').read_text());two=json.loads((ROOT.parent/'coarse_two_cover_affine'/'normal.json').read_text());five=json.loads((ROOT.parent/'five_profile_types'/'normal.json').read_text());graphs=[];index={}
 def add(R):
  key=tuple(R)
  if key not in index:index[key]=len(graphs);graphs.append(R)
  return index[key]
 coeff=[]
 for hr,tr in zip(h['small_model_controls'],two['small_model_controls']):
  need(hr['quotient_rows']==tr['quotient_rows'] and hr['block_queries']==tr['block_queries'],'same independent coefficient domain');Q=hr['quotient_rows'];m=len(Q)
  for hx,tx in zip(hr['extracted_coefficients'],tr['extracted_coefficients']):
   a,c=hx['block_query'];need(tx['block_query']==[a,c],'same fixed query');ids=[add(u.extend(P,a,c)) for P in (Q,u.op(Q))];coeff.append({'name':hr['name'],'quotient_rows':Q,'block_query':[a,c],'cycle_ids':ids,'Hamilton_coefficients':hx['coefficients'],'two_cover_coefficients':tx['coefficients']})
 actual=[];querylist=list(product(range(8),repeat=2))
 for r,tw in zip(five['canonical_records'],two['canonical_controls']):
  need(r['types']==tw['types'],'same canonical order');blocks=r['child_rows'];Q=[2,4,1];unit=pow(3,r['types'].count('C'),8)
  for a,c in [(1,0),(1,2),(1,3)]:
   ix=querylist.index((a,c));p,q=r['physical_queries'][ix];actual.append({'kind':'canonical','types':r['types'],'block_query':[a,c],'child_rows':blocks,'odd_H_product_mod8':unit,'physical_H_q2':tw['values_H_q2_source_target'][ix],'bare_ids':[add(u.extend(P,a,c)) for P in (Q,u.op(Q))],'physical_ids':[add(u.extend(T,p,q)) for T in (r['source_old_rows'],r['target_old_rows'])]})
 for r in two['larger_child_controls']:
  Q=r['quotient_rows'];blocks=r['child_rows'];T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks)
  for (a,c),(p,q),vv in zip(r['block_queries'],r['physical_queries'],r['values_H_q2_source_target']):actual.append({'kind':r['name'],'types':r['types'],'block_query':[a,c],'child_rows':blocks,'odd_H_product_mod8':r['odd_H_product_mod8'],'physical_H_q2':vv,'bare_ids':[add(u.extend(P,a,c)) for P in (Q,u.op(Q))],'physical_ids':[add(u.extend(V,p,q)) for V in (T,U)]})
 lines=[str(len(graphs))]
 for R in graphs:lines += [str(len(R)),' '.join(map(str,R))]
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'cycles';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'cycles.cpp'),'-o',str(exe)],check=True);counts=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(len(counts)==len(graphs),'all complete Hamilton-cycle counts consumed');entries=0
 for r in coeff:
  zs=[counts[i] for i in r.pop('cycle_ids')];r['bare_cycle_counts']=zs;m=len(r['quotient_rows'])
  for z,hh,tt in zip(zs,r['Hamilton_coefficients'],r['two_cover_coefficients']):
   need((tt['B0']-(2*hh['K0']-4-4*(m%2)*z))%8==0,'base formula modulo8')
   for k,b,c in zip(hh['K2'],tt['B2'],tt['B3']):need((b-2*k-2*z)%4==0 and (c-2*z)%4==0,'all double/triple-run coefficients from one cycle parity');entries+=1
 for r in actual:
  bs=[counts[i] for i in r.pop('bare_ids')];ac=[counts[i] for i in r.pop('physical_ids')];r['bare_cycle_counts']=bs;r['physical_cycle_counts']=ac;need(all((a-b)%2==0 for a,b in zip(bs,ac)),'literal physical Hamilton-cycle parity collapse');eps=sum(map(len,r['child_rows']))%2;unit=r['odd_H_product_mod8'];vv=r['physical_H_q2']
  for j,z in enumerate(bs):H,q2=vv[2*j:2*j+2];phi=H*pow(unit,-1,8)%8;psi=q2*pow(unit,-1,8)%8;need((psi-(2*phi-4-4*eps*z))%8==0,'actual normalized two-cover formula including odd scalar');
 # Negative ambient control is deliberately retained.
 Q=[2,4,1];blocks=[u.tt(2),u.tt(2),u.tt(3)];T,offs=u.blow(Q,blocks);U,_=u.blow(u.op(Q),blocks);a,c=phys(1,blocks,offs),phys(2,blocks,offs);cores=[u.extend(V,a,c) for V in (T,U)];profiles=[u.profile(V) for V in cores];oldA,_=u.blow(Q,[T,[0],[0]]);oldB,_=u.blow(Q,[U,[0],[0]]);outside=1<<(len(T)+1);ambient=[u.hp(u.extend(V,a|outside,c|outside)) for V in (oldA,oldB)];source=next(r for r in actual if r['types']==['P','P','O'] and r['block_query']==[1,2]);need([profiles[0][k]-profiles[1][k] for k in (1,2,3)]==[-68,-396,-1116],'literal odd-core profile boundary');need(ambient[0]-ambient[1]==-600 and (source['bare_cycle_counts'][0]-source['bare_cycle_counts'][1])%2==1,'bare cycle mismatch and modular loss')
 deps=['mod8_transfer/PROOF.md','five_profile_types/PROOF.md','five_profile_types/normal.json','coarse_mod8_affine/PROOF.md','coarse_mod8_affine/normal.json','coarse_two_cover_affine/PROOF.md','coarse_two_cover_affine/normal.json','uniform_block_starts/verify.py']
 return {'status':'PASS','source_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','cycles.cpp']},'dependency_pins':{f:hashlib.sha256((ROOT.parent/f).read_bytes()).hexdigest() for f in deps},'distinct_literal_cycle_graphs':len(graphs),'cycle_graphs':[{'rows':R,'Hamilton_cycles':z} for R,z in zip(graphs,counts)],'all_coefficient_controls':coeff,'individual_B_and_C_relations':entries,'physical_cycle_collapse_controls':actual,'odd_total_negative_control':{'child_rows':blocks,'physical_query':[a,c],'bare_cycle_counts':source['bare_cycle_counts'],'full_core_profiles':profiles,'ambient_H_counts':ambient,'ambient_gap':ambient[0]-ambient[1]},'scope':'General child-uniform Hamilton-cycle parity collapse and a bare-quotient test for the odd ambient correction. This does not provide a local separator or prove non-annihilation for arbitrary quotients.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
