#!/usr/bin/env python3
"""Complete producer bank replay, dependency checks and physical record validation."""
import collections,csv,gzip,hashlib,json,pathlib,subprocess,tempfile
import prepare as p
P=p.P;B=p.B;S=p.S
def check(t,m):
 if not t:raise RuntimeError(m)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def tournament(r):
 n=len(r);return all(0<=z<1<<n and not(z>>i&1) for i,z in enumerate(r)) and all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))
def bank_input(cs,qs):
 lines=[f"{r['node']} {a} {c} "+' '.join(map(str,r['rows'])) for r in cs for a,c in qs]
 return str(len(lines))+'\n'+'\n'.join(lines)+'\n'
def main():
 for name,m in p.load(P/'DEPENDENCIES.json').items():
  f=P.parent/name;check(f.stat().st_size==m['bytes'] and sha(f)==m['sha256'],'dependency '+name)
 cs=p.load(B/'candidates.json')['candidates'];check(len(cs)==2348,'quotient domain');ids=[r['node'] for r in cs];byid={r['node']:r for r in cs};qs,witnesses,miss=p.derive()
 check(len(qs)==16 and qs[14:]==[(0,2),(1,17)],'queries');check(len(set(qs+[(1<<i,0) for i in range(9)]))==24,'endpoint union')
 check(miss==p.load(P/'root_initial_misses.json') and len(miss)==11,'initial misses')
 check(p.load(P/'root_battery.json')==[dict(index=i,a=a,c=c) for i,(a,c) in enumerate(qs)],'battery')
 rs=p.tab(P/'root_residues.tsv');es=p.tab(P/'root_selected_residues.tsv');check(len(rs)==2348*14 and len(es)==22,'bank dimensions')
 check([(int(r['node']),int(r['a']),int(r['c'])) for r in rs]==[(n,a,c) for n in ids for a,c in qs[:14]],'full bank ordering')
 selected_ids=sorted({r['node'] for r in miss});selected_cs=[r for r in cs if r['node'] in selected_ids]
 check([(int(r['node']),int(r['a']),int(r['c'])) for r in es]==[(n,a,c) for n in selected_ids for a,c in qs[14:]],'selected bank ordering')
 for r in rs+es:
  for i in range(46):
   a,b=int(r[f'S{i}']),int(r[f'T{i}']);check(0<=a<8 and 0<=b<8 and a%2==b%2==1 and (a-b)%8 in [0,4],'residue range')
 inherited={int(r['node']):r for r in p.tab(B/'residues.tsv')}
 for index,r in enumerate(rs):
  k=index%14
  if k<12:
   old=inherited[int(r['node'])];check(all(int(r[f'{side}0'])*3%8==int(old[f'{side}{k}']) for side in 'ST'),'inherited baseline')
 text=bank_input(cs,qs[:14]);extra=bank_input(selected_cs,qs[14:]);check(text==(P/'root_bank_input.txt').read_text(),'bank input');check(extra==(P/'root_selected_input.txt').read_text(),'selected input')
 with tempfile.TemporaryDirectory(prefix='size29-verify-') as tmp:
  exe=pathlib.Path(tmp)/'bank';subprocess.run(['g++','-O3','-std=c++17',str(P/'root_bank.cpp'),'-o',str(exe)],check=True)
  for txt,name in [(text,'root_residues.tsv'),(extra,'root_selected_residues.tsv')]:
   result=subprocess.run([str(exe)],input=txt.encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=True);check(result.stdout==(P/name).read_bytes(),'clean recurrence replay '+name)
 # Complete local classification is separately independently verified; establish exact modular inputs here.
 five=p.load(P/'five_local/five_child_records.json');check(len(five)==12 and sum(r['labeled_count'] for r in five if r['nontransitive'])==904,'five classification')
 check({tuple(r['normalized_mod8']) for r in five if r['nontransitive']}=={(1,6,6,0,0),(1,2,2,0,0)},'five profiles')
 four=p.load(S/'four_child_records.json')
 check({tuple(q*pow(r['profile'][0],-1,8)%8 for q in r['profile']) for r in four if r['nontransitive']}=={(1,6,4,0)},'four profiles')
 check(all(r['profile']==t['profile'] for r in four if r['nontransitive'] for t in r['targets']),'four families')
 # The new five-child obstruction is retained rather than incorrectly assuming all profiles match.
 r8=next(r for r in five if r['code']==8);check(r8['profile']==[9,66,186,240,120] and not r8['all_family_profiles_equal'],'five obstruction')
 check(len(witnesses)==126792,'witness count');tsv=p.tab(P/'root_witnesses.tsv');check(len(tsv)==len(witnesses),'witness table size')
 for actual,w in zip(tsv,witnesses):
  restored={k:(json.loads(actual[k]) if isinstance(v,list) else int(actual[k]) if isinstance(v,int) else actual[k]) for k,v in w.items()};check(restored==w,'witness row');check((w['normalized_source']-w['normalized_target'])%8==4,'normalized gap')
 count=0
 with gzip.open(P/'root_physical_witnesses.jsonl.gz','rt') as f:
  for w in witnesses:
   line=f.readline();check(bool(line),'physical missing');r=json.loads(line);expected=p.physical(w,byid[w['node']]['rows']);check(r==expected,'physical record')
   check(tournament(r['augmented_source_rows']) and tournament(r['augmented_target_rows']),'physical tournament');check((r['source_residue']-r['target_residue'])%8==4,'physical gap');count+=1
  check(f.read()=='','physical extra')
 summary=p.load(P/'root_summary.json');check(summary['physical_witness_pairs']==count and summary['endpoint_union_queries']==24,'summary')
 names=p.load(P/'RELEASE_FILES.json');hashes={n:sha(P/n) for n in names if n not in ['receipt.json','receipt.normal.json','receipt.optimized.json','SHA256SUMS']}
 receipt=dict(status='PASS',scope='Producer exact size29 boundary; source-relative >=30 consequence pending independent review; no unrestricted injectivity',endpoint_zero_quotients=2348,initial_query_rows=len(rs),initial_individual_residues=len(rs)*92,selected_query_rows=len(es),selected_individual_residues=len(es)*92,exceptional_five_cases=21132,two_four_cases=84528,inherited_five_cases=21132,initial_misses=len(miss),root_queries=16,endpoint_union_queries=24,physical_witness_pairs=count,physical_graphs=2*count,normalized_gap=4,modulus=8,sha256=hashes)
 (P/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='sha256'},indent=2))
if __name__=='__main__':main()
