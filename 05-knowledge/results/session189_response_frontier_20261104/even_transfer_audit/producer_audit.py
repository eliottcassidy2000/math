"""Cross-check all producer kernels and literal cases, retaining all physical labels."""
from pathlib import Path
import json,hashlib,subprocess,tempfile,math
from kernel_audit import tournament
from physical_audit import profile,substitute,op,augment
D=Path(__file__).resolve().parent;P=D.parent/'even_transfer'
if not P.exists():P=D.parent/'session189-even-transfer'
def need(c,s):
 if not c:raise RuntimeError(s)
def main():
 rows=json.loads((D/'kernel_recount.json').read_text());lines=(P/'kernels.txt').read_text().splitlines();need(len(rows)==len(lines),'all retained kernel rows')
 for r,line in zip(rows,lines):
  n,code,i,k,*v=map(int,line.split());need(r['rows']==tournament(n,code) and r['distinguished']==i and r['multiplicity']==k,'all fixed quotient physical labels');need(v==r['K_1_through_length']+[0]*(9-len(r['K_1_through_length'])),'all independent exact kernels')
 data=json.loads((P/'physical_records.json').read_text());cases=[];labels=[];qcache={}
 for ind,r in enumerate(data['records']+data['complete_outside_domain_bank']):
  blocks=[data['child_blocks'][i] for i in r['child_type_indices']] if 'child_type_indices' in r else r['child_rows'];Q=r['Q'];a,c=r['block_query'];ambient=r.get('ambient',False)
  for side,s in enumerate(r['sides']):
   children=s['target_child_rows'];need(len(children)==len(blocks),'every physical child')
   for C,E in zip(children,blocks):
    for R in [C,E]:
     key=tuple(R)
     if key not in qcache:qcache[key]=profile(R)
    need(qcache[tuple(C)]==qcache[tuple(E)],'corresponding exact complete child profile')
   old,sets=substitute(op(Q) if side else Q,children);need(sets==s['child_sets'],'exact physical child sets')
   A=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);C=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);ext=[]
   if ambient:
    n=len(old);old,_=substitute([2,4,1],[old,[0],[0]]);A|=1<<(n+1);C|=1<<(n+1);ext=[n,n+1]
   need([A,C]==s['physical_query'] and ext==s['exterior_names'] and [len(old),len(old)+1]==s['port_names'],'exact physical masks and ports/exterior')
   physical=augment(old,A,C);need(physical==s['augmented_rows'],'all retained physical arcs')
   need(all(not(physical[u]>>u&1) and all(((physical[u]>>v)&1)+((physical[v]>>u)&1)==1 for v in range(len(physical)) if u!=v) for u in range(len(physical))),'actual physical tournament')
   cases.append(physical);labels.append((ind,side))
 need(math.factorial(18)*math.comb(17,4)<1<<64,'exact first-five arithmetic bound')
 with tempfile.TemporaryDirectory() as tmp:
  binary=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'literal.cpp'),'-o',str(binary)],check=True)
  inp=str(len(cases))+'\n'+''.join(str(len(R))+'\n'+' '.join(map(str,R))+'\n' for R in cases)
  raw=subprocess.check_output([str(binary)],input=inp.encode()).decode()
 values=[list(map(int,line.split())) for line in raw.splitlines()];need(len(values)==len(cases),'all independent literal cases');at=0;records=[]
 for ind,r in enumerate(data['records']+data['complete_outside_domain_bank']):
  vv=values[at:at+2];at+=2;need(all(q==s['ordered_covers'] for q,s in zip(vv,r['sides'])),'all retained exact physical covers');delta=[a-b for a,b in zip(*vv)];need(delta==r['contrast'],'all exact physical contrasts')
  records.append({'record':ind,'block_query':r['block_query'],'ambient':r.get('ambient',False),'source_rows':r['sides'][0]['augmented_rows'],'target_rows':r['sides'][1]['augmented_rows'],'source_q_1_through_5':vv[0],'target_q_1_through_5':vv[1],'contrast':delta})
  if ind<len(data['records']):
   if not r.get('ambient',False):need((delta[1]-2*delta[0])%16==0 and all(x%16==0 for x in delta[2:]),'complete transfer modulus with factorial tail')
   else:need((delta[0]-5*data['records'][r['local_record']]['contrast'][0])%16==0,'same physical cyclic-ambient query multiplier')
 outside=records[len(data['records']):];need(any((r['contrast'][1]-2*r['contrast'][0])%16 or any(x%16 for x in r['contrast'][2:]) for r in outside),'outside-domain boundary is not universalized')
 out=json.dumps(records,separators=(',',':'))+'\n';p=D/'producer_physical_recount.json'
 if p.exists():need(p.read_text()==out,'retained independent physical recount stable')
 else:p.write_text(out)
 return dict(status='PASS',complete_kernel_records=len(rows),padded_exact_kernel_coordinates=9*len(rows),physical_sides=len(cases),exact_physical_cover_coordinates=5*len(cases),complete_child_profiles=len(qcache),maximum_augmented_order=max(map(len,cases)),producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest(),producer_kernel_sha256=hashlib.sha256((P/'kernels.txt').read_bytes()).hexdigest(),producer_physical_sha256=hashlib.sha256((P/'physical_records.json').read_bytes()).hexdigest(),independent_data_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['producer_audit.py','literal.cpp','kernel_audit.py','physical_audit.py']})
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
