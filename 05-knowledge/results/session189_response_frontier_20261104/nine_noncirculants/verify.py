from pathlib import Path
from array import array
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(n):return hashlib.sha256((D/n).read_bytes()).hexdigest()
def valid(T):return all(not(r>>i&1) and all(((r>>j)&1)+((T[j]>>i)&1)==1 for j in range(len(T)) if i!=j) for i,r in enumerate(T))
def modules(T):
 return [s for s in range(1,(1<<len(T))-1) if s.bit_count()>1 and all((T[i]&s) in (0,s) for i in range(len(T)) if not(s>>i&1))]
def main():
 logs={}
 for script in ['anchors.py','affine.py','certificates.py','build_banks.py','physical_mod8.py','residual_checks.py']:logs[script]=json.loads(subprocess.check_output([sys.executable,str(D/script)]))
 from affine import aug,op,solve
 from physical_helpers import build
 anchors=json.loads((D/'anchor_data.json').read_text());data=json.loads((D/'analysis.json').read_text());certs=json.loads((D/'linear_certificates.json').read_text());need(len(anchors)==len(data)==len(certs)==4,'exact four displayed specimens');expected_blind=[[254,511],[511],[511],[54,511]];expected_H=[2251,2239,3027,2479];query_records=[];cube_checks=[]
 for mi,(q,a,cert) in enumerate(zip(data,anchors,certs)):
  Q=q['rows'];need(Q==a['rows']==cert['rows'] and valid(Q) and not modules(Q),'physical aligned prime quotient');need(sum(a['S'])==expected_H[mi] and sum(a['E'])==expected_H[mi],'exact distinct Hamilton counts');need(len(set(r.bit_count() for r in Q))>1,'nonregular hence noncirculant');need(not any(a['p']) and not any(a['sigma']),'both endpoint obstructions');need(a['blind']==expected_blind[mi] and len(a['patterns'])==512,'exact anchor residual scope')
  for V,p in enumerate(a['patterns']):
   need(V==p['V'],'physical singleton-position indexing');d=[(a['d0'][i]+2*sum(a['d2'][i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];t=[(a['S'][i]-a['E'][i]+2*sum(a['start'][i][j] for j in range(9) if V>>j&1))%4 for i in range(9)];need(d==p['first'] and [(2*t[i]+2*d[i])%8 for i in range(9)]==p['second'],'all singleton-sensitive affine anchors')
  need(len(q['queries'])==206,'all selected actual query pairs');equations=[];safe=[]
  for qi,r in enumerate(q['queries']):
   S,U=r['source'],r['target'];need(len(S)==len(U)==56,'complete run and cycle channel vectors');delta=[s-u for s,u in zip(S,U)];base=(delta[0]+2*sum(delta[1:10])+4*sum(delta[19:55]))%8;need(base%4==0 and base==r['baseline'],'affine baseline');need(all(x%2==0 for x in delta[10:19]),'return segment parity');bits=[x%2 for x in delta[1:10]]+[(x//2)%2 for x in delta[10:19]];need(bits==r['linear_bits'],'all profile derivatives');need(r['cycle_safe']==(delta[55]%2==0),'actual cycle parity test');eq=(sum(b<<j for j,b in enumerate(bits)),base//4);equations.append(eq)
   if r['cycle_safe']:safe.append(eq)
   query_records.append(dict(model_index=mi,query_index=qi,query=r['query'],source_augmented_rows=aug(Q,*r['query']),target_augmented_rows=aug(op(Q),*r['query']),port_names=[9,10],source_rooted_cycle_count=S[-1],target_rooted_cycle_count=U[-1]))
  need(solve(equations)==q['all_queries'] and solve(safe)==q['safe_queries'] and q['safe_queries']['inconsistent'],'full safe affine inconsistency')
  ids=cert['safe_certificate']['contradiction_query_indices'];need(len(ids)==[7,5,1,4][mi],'small fixed battery size');chosen=[];row=rhs=0
  for i in ids:
   r=q['queries'][i];need(r['cycle_safe'],'each certificate query cycle safe');eq=(sum(b<<j for j,b in enumerate(r['linear_bits'])),r['baseline']//4);chosen.append(eq);row^=eq[0];rhs^=eq[1]
  need(row==0 and rhs==1,'explicit affine XOR contradiction');zeros=[x for x in range(1<<18) if all(((r&x).bit_count()%2)==b for r,b in chosen)];need(not zeros,'all262144 actual profile-bit assignments checked');cube_checks.append(dict(model_index=mi,assignments_checked=1<<18,blind_assignments=zeros))
  if mi==2:
   r=q['queries'][ids[0]];need(r['query']==[3,132] and r['baseline']==4 and not any(r['linear_bits']) and r['source'][-1]%2==r['target'][-1]%2==0,'single universal nontrivial-profile witness')
  for V in a['blind']:
   if V==511:continue
   r=next(r for r in q['queries'] if r['query']==[1,2]);d=[s-u for s,u in zip(r['source'],r['target'])];need((d[0]+2*sum(d[j+1] for j in range(9) if V>>j&1))%4==2,'complete singleton residual mod4 formula')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'endpoint_check.cpp'),'-o',str(D/'endpoint_check')],check=True);raw=subprocess.check_output([str(D/'endpoint_check')],input=('4\n'+''.join(' '.join(map(str,q['rows']))+'\n' for q in anchors)).encode()).decode();(D/'endpoint_values.txt').write_text(raw)
 for q,line in zip(anchors,raw.splitlines()):need(list(map(int,line.split()))==sum(q['joint_endpoints'],[]),'independent exhaustive joint endpoints')
 banks=json.loads((D/'named_banks.json').read_text());need(len(banks)==8,'complete source and target named banks');full_hashes=[];degrees={d:0 for d in range(5)}
 for b in banks:
  masks=b['physical_masks'];need(masks==[m for m in range(512) if m.bit_count()<=2],'full named bidegree basis');cube=[0]*(1<<18)
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):cube[a+512*c]=b['coefficients'][i][j];degrees[a.bit_count()+c.bit_count()]+=1
  for bit in range(18):
   step=1<<bit
   for block in range(0,1<<18,2*step):
    for j in range(block+step,block+2*step):cube[j]+=cube[j-step]
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):need(cube[a+512*c]==b['basis_values'][i][j],'all named basis values recovered')
  need(all(0<v<=39916800 and v%2 for v in cube),'full reconstructed named cube Hamilton parity and bounds');arr=array('I',cube)
  if sys.byteorder!='little':arr.byteswap()
  full_hashes.append(dict(name=b['name'],side=b['side'],values=1<<18,little_endian_u32_sha256=hashlib.sha256(arr.tobytes()).hexdigest()))
 physical=json.loads((D/'physical_mod8_records.json').read_text());need(len(physical)==1870,'all exact physical interpolation graphs')
 for r in physical:
  query=data[r['model_index']]['queries'][r['query_index']]['query'];expected=build(r['quotient_rows'],r['child_sizes'],*query);need(all(r[k]==v for k,v in expected.items()),'literal interpolation physical names and masks');need(valid(r['augmented_rows']) and r['H']%2==1,'actual interpolation tournament and count')
 residuals=json.loads((D/'residual_records.json').read_text());need(len(residuals)==4,'two local and two ambient residual pairs')
 for r in residuals:
  for s in r['sides']:
   T=s['augmented_rows'];x,y=s['port_names'];aa,cc=s['physical_query'];need(valid(T) and T[x]>>y&1,'actual lifted physical tournament with port edge');old=set(range(len(T)))-{x,y};need(all(((T[i]>>x)&1)==((aa>>i)&1) and ((T[i]>>y)&1)==((cc>>i)&1) for i in old),'actual port masks on every physical nonport name')
   if r['ambient']:
    u,v=s['exterior_names'];root=old-{u,v}|{x,y};need(all(T[i]>>u&1 and T[v]>>i&1 for i in root) and T[u]>>v&1,'both ports inherit fixed cyclic root incidence')
  need(r['contrast']==r['sides'][0]['H']-r['sides'][1]['H'] and r['contrast']%4==2,'literal local/ambient residual gap')
 (D/'query_records.json').write_text(json.dumps(query_records,separators=(',',':'))+'\n');(D/'profile_cube_checks.json').write_text(json.dumps(cube_checks,indent=2)+'\n');(D/'cube_hashes.json').write_text(json.dumps(full_hashes,indent=2)+'\n')
 files=['PROOF.md','anchors.py','anchors.cpp','endpoint_check.cpp','affine.py','runs.cpp','certificates.py','build_banks.py','physical_h.cpp','physical_helpers.py','physical_mod8.py','residual_checks.py','anchor_data.json','anchor_values.txt','endpoint_values.txt','analysis.json','run_values.txt','linear_certificates.json','named_banks.json','named_values.txt','cube_hashes.json','physical_mod8_records.json','physical_mod8_values.txt','residual_records.json','residual_values.txt','query_records.json','profile_cube_checks.json']
 print(json.dumps(dict(status='PASS',scope='Four explicit noncirculant prime9 quotients, arbitrary matching children and fixed ambient; not an order9 census',quotients=4,joint_endpoint_coordinates=324,anchor_primitives=684,singleton_patterns=2048,exact_run_and_cycle_counts=4*206*2*56,certifying_safe_queries=[7,5,1,4],profile_assignments_checked=4*(1<<18),physical_interpolation_graphs=1870,residual_physical_sides=8,complete_named_bank_coordinates=16928,named_degree_counts=degrees,reconstructed_full_named_values=8*(1<<18),sha256={f:sha(f) for f in files}),sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
