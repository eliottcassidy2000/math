#!/usr/bin/env python3
"""Exact producer replay; every check remains active under python -O."""
import collections,csv,gzip,hashlib,json,pathlib,subprocess,tempfile
import prepare as p
import prepare_local as local
P=pathlib.Path(__file__).resolve().parent; B=P.parent/'session189-c3-frontier'

def check(t,m):
 if not t:raise RuntimeError(m)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def tournament(r):
 n=len(r)
 return all(0<=z<1<<n and not(z>>i&1) for i,z in enumerate(r)) and all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))
def bank_input(qs,cs):return str(len(qs))+'\n'+'\n'.join(f'{a} {c}' for a,c in qs)+'\n'+str(len(cs))+'\n'+'\n'.join(str(r['node'])+' '+' '.join(map(str,r['rows'])) for r in cs)+'\n'
def main():
 for name,meta in p.load(P/'DEPENDENCIES.json').items():
  path=P.parent/name;check(path.stat().st_size==meta['bytes'] and sha(path)==meta['sha256'],'dependency '+name)
 cs=p.load(B/'candidates.json')['candidates'];byid={r['node']:r for r in cs};check(len(cs)==2348,'candidate size')
 qs=[(r['a'],r['c']) for r in p.load(P/'battery.json')];check(qs==[(0,1),(1,0),(1,2),(1,3),(3,1),(3,3),(3,5),(0,7),(0,18),(40,256),(4,3),(130,128),(65,9),(3,32)],'battery');check(len(set(qs+[(1<<i,0) for i in range(9)]))==22,'endpoint union')
 check([(r['a'],r['c']) for r in p.load(P/'selected_extras.json')]==qs[12:],'selected extra query provenance')
 initial=p.tab(P/'residues.tsv');extra=p.tab(P/'selected_residues.tsv');check(len(initial)==21132 and len(extra)==117,'bank lengths')
 check([(int(r['node']),int(r['position'])) for r in initial]==[(c['node'],j) for c in cs for j in range(9)],'placement completeness')
 misses=[];chosen=[];hist=collections.Counter();ex={(int(r['node']),int(r['position'])):r for r in extra}
 for r in initial:
  node,j=int(r['node']),int(r['position']);values=[(int(r[f'S{k}']),int(r[f'T{k}'])) for k in range(12)]
  check(all(0<=s<8 and 0<=t<8 and s%2==t%2==1 and (s-t)%8 in (0,4) for s,t in values),'initial residues')
  hits=[(k,s,t) for k,(s,t) in enumerate(values) if s!=t]
  if not hits:
   misses.append(dict(node=node,position=j));r2=ex[node,j];hits=[(12+k,int(r2[f'S{k}']),int(r2[f'T{k}'])) for k in range(2) if r2[f'S{k}']!=r2[f'T{k}']]
  check(bool(hits),'unseparated');k,s,t=hits[0];a,c=qs[k];hist[k]+=1;chosen.append(dict(node=node,position=j,query_index=k,a=a,c=c,physical_a=p.expand(a,j),physical_c=p.expand(c,j),normalized_source=s,normalized_target=t))
 check(misses==p.load(P/'initial_misses.json') and len(misses)==13,'initial misses')
 nodes={r['node'] for r in misses};ec=[r for r in cs if r['node'] in nodes]
 check([(int(r['node']),int(r['position'])) for r in extra]==[(c['node'],j) for c in ec for j in range(9)],'extra completeness')
 check(all(0<=int(r[f'{side}{k}'])<8 and int(r[f'{side}{k}'])%2==1 for r in extra for side in 'ST' for k in range(2)),'extra residue range')
 check([{k:int(v) for k,v in r.items()} for r in p.tab(P/'witnesses.tsv')]==chosen,'witness table')
 inp=bank_input(qs[:12],cs);einp=bank_input(qs[12:],ec)
 check(inp==(B/'bank_input.txt').read_text(),'incoming bank inputs');check(einp==(P/'selected_input.txt').read_text(),'extra bank inputs')
 with tempfile.TemporaryDirectory(prefix='size28-verify-') as tmp:
  exe=pathlib.Path(tmp)/'bank';subprocess.run(['g++','-O3','-std=c++17',str(P/'bank.cpp'),'-o',str(exe)],check=True)
  for text,name in [(inp,'residues.tsv'),(einp,'selected_residues.tsv')]:
   result=subprocess.run([str(exe)],input=text.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);check(result.stdout==(P/name).read_bytes(),'clean recurrence replay '+name)
 records=p.load(P/'four_child_records.json');check(len(records)==64,'four classification size');profile_hist=collections.Counter();family_count=0
 for code,r in enumerate(records):
  rows=local.rows(code);profile=local.profile(rows);modules=local.modules(rows)
  check(r['code']==code and r['rows']==rows and r['profile']==profile and r['degrees']==sorted(z.bit_count() for z in rows),'four source');check(r['nontransitive']==(profile[0]>1) and r['modules']==modules,'four modules')
  expected=[]
  for f in range(1<<len(modules)):
   fam=[m for j,m in enumerate(modules) if f>>j&1];target=local.xor(rows,fam);target_profile=local.profile(target);expected.append(dict(family=fam,rows=target,profile=target_profile))
   if r['nontransitive']:check(profile==target_profile,'four family profile');family_count+=1
  check(r['targets']==expected,'four targets')
  if r['nontransitive']:profile_hist[tuple(profile)]+=1
 check(dict(profile_hist)=={(3,18,36,24):16,(5,22,36,24):24} and family_count==160,'nontransitive classification')
 for typ,(rows,profile) in p.TYPES.items():
  check(local.profile(rows)==profile,'canonical profile');inv=pow(profile[0],-1,8);check([z*inv%8 for z in profile]==[1,6,4,0],'normalized profile')
 controls=p.load(P/'named_controls.json');check(controls==[p.named_control(t) for t in p.TYPES],'complete named controls')
 cubic=quartic=0
 for ctl in controls:
  for side in ['source','converse']:
   coeff=ctl[side]['coefficients'];values=ctl[side]['values'];check(len(coeff)==len(values)==256,'named dimensions')
   check(all(not v for m,v in enumerate(coeff) if m.bit_count()>4),'degree bound control')
   cubic+=sum(v!=0 and m.bit_count()==3 for m,v in enumerate(coeff));quartic+=sum(v!=0 and m.bit_count()==4 for m,v in enumerate(coeff))
   rebuilt=coeff[:]
   for j in range(8):
    for m in range(256):
     if m>>j&1:rebuilt[m]+=rebuilt[m^(1<<j)]
   check(rebuilt==values,'named polynomial evaluation')
 count=0
 with gzip.open(P/'physical_witnesses.jsonl.gz','rt') as f:
  for w in chosen:
   q=byid[w['node']]['rows'];op=[511^z^(1<<i) for i,z in enumerate(q)]
   for typ,(_,profile) in p.TYPES.items():
    line=f.readline();check(bool(line),'physical missing');r=json.loads(line);expected=dict(w,type=typ,profile=profile,source_residue=w['normalized_source']*profile[0]%8,target_residue=w['normalized_target']*profile[0]%8,augmented_source_rows=p.aug(p.inflate(q,w['position'],typ),w['physical_a'],w['physical_c']),augmented_target_rows=p.aug(p.inflate(op,w['position'],typ),w['physical_a'],w['physical_c']))
    check(r==expected,'physical record');check(tournament(r['augmented_source_rows']) and tournament(r['augmented_target_rows']),'physical tournament');check((r['source_residue']-r['target_residue'])%8==4,'physical gap');count+=1
  check(f.read()=='','physical extra')
 names=['PROOF.md','README.md','bank.cpp','battery.json','residues.tsv','selected_residues.tsv','selected_input.txt','selected_extras.json','initial_misses.json','witnesses.tsv','physical_witnesses.jsonl.gz','four_child_records.json','named_controls.json','prepare.py','prepare_local.py','verify.py','DEPENDENCIES.json']
 receipt=dict(status='PASS',scope='Producer finite prime-nine/eight-C3/one-nontransitive-four theorem; source-relative <=28 exclusion uses pinned reviewed inputs and exact one-singleton exterior argument; pending independent review',endpoint_zero_quotients=2348,placements=21132,initial_queries=12,initial_misses=13,additional_queries=2,endpoint_union_queries=22,initial_normalized_residues=2*12*21132,additional_normalized_residues=2*2*117,physical_witness_pairs=count,physical_graphs=2*count,nontransitive_four_labeled=40,nontransitive_four_family_outputs=160,first_query_histogram=[hist[i] for i in range(14)],named_control_response_vectors=6,named_control_coefficients=1536,nonzero_named_cubic_coefficients=cubic,nonzero_named_quartic_coefficients=quartic,sha256={n:sha(P/n) for n in names})
 (P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='sha256'},indent=2))
if __name__=='__main__':main()
