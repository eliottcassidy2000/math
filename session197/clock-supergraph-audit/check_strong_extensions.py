#!/usr/bin/env python3
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib,json,subprocess,tempfile,sys,os
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'session197-orbit-wildcard-graph-audit'))
from check_graph import check,h_walk_inclusion_exclusion,literal_h
from check_converse_mass import degree_canonical,relabel
sys.path.insert(0,str(P.parent/'session197-prime-atoms'))
from verify_minimum_ternary import strong
SOURCE=P.parent/'session197-hybrid-constructions'/'clock-supergraph-barrier'
HASHES={'PROOF.md':'0083018e9156f1105b73dbafc4531eb54769c0ed68c75e90f50ea3e1e43411a6','census.cpp':'733e7eb551f3be49f811622cd6b79d30e02d29eb7bbf053c85f449590bcc603e','verify.py':'cf80900fc73fee79aa286f283fa02d413f68a93252113be0321522d0d1c9f10d','VALIDATION.json':'f502eb547e2303662e82c43475c491766801626ff3a28fac81fb884426bfe017'}
SEED=[2,12,145,149,99,143,175,19]

def main():
 for name,h in HASHES.items():check(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==h,('frozen',name))
 with tempfile.TemporaryDirectory(prefix='strong-audit-',dir=P) as tmp:
  exe=Path(tmp)/'search';bf=Path(tmp)/'bank.jsonl'
  flags=['-O0','-D_GLIBCXX_ASSERTIONS'] if os.environ.get('STRONG_AUDIT_DEBUG') else ['-O3']
  subprocess.run(['g++','-std=c++17',*flags,str(P/'check_strong_extensions.cpp'),'-o',str(exe)],check=True)
  run=subprocess.run([str(exe),str(bf)],check=True,text=True,capture_output=True);stats=json.loads(run.stdout)
  bank=[json.loads(line) for line in bf.read_text().splitlines()]
  check(dict(Counter(z['n'] for z in bank))=={8:1,9:27,10:3},'independent complete strong banks')
  target=degree_canonical(SEED);checks=[]
  for z in bank:
   rows=z['rows'];n=len(rows);check(strong(rows),'strong retained graph')
   check(degree_canonical(rows)==tuple(rows),'exact canonical representative')
   H=h_walk_inclusion_exclusion(rows)[0];check(H==z['H'] and H<=541 and H!=541,'independent walk count')
   copies=[list(V) for V in combinations(range(n),8) if degree_canonical(relabel(rows,V))==target]
   check(bool(copies),'actual embedded seed remains')
   checks.append({'order':n,'H':H,'seed_vertex_subsets':copies})
  bankhash=hashlib.sha256(bf.read_bytes()).hexdigest()
 G=[2,64,11,51,39,135,61,95];ports=[0,1,2,5]
 outs=[251,35];values=[]
 for out,want in zip(outs,[539,541]):
  rows=[row|((not(out>>v&1))<<8) for v,row in enumerate(G)]+[out]
  check(strong(rows),'strong rebased child')
  check(h_walk_inclusion_exclusion(rows)[0]==literal_h(rows)==want,'independent literal and walk child count')
  values.append({'rows':rows,'H':want})
 check(h_walk_inclusion_exclusion(G)[0]==literal_h(G)==179 and strong(G),'independent rebased parent')
 check(degree_canonical(G)!=target,'equal-count parent is a different graph')
 check(h_walk_inclusion_exclusion(SEED)[0]==179,'original seed H')
 plus2=[row|((not(19>>v&1))<<8) for v,row in enumerate(SEED)]+[19]
 check(h_walk_inclusion_exclusion(plus2)[0]==539,'original +2 child')
 check(all(541%d for d in range(2,24)) and 5**11<=541**3<5**12,'prime and exact strong bound')
 check(3*179+2**8-2==791>541,'seed-preserving mixed-extension barrier')
 print(json.dumps({'status':'PASS','source_hashes':HASHES,'independent_strong_generation':stats,'strong_bank_sha256':bankhash,'independent_retained_walk_and_embedding_checks':checks,'rebased_parent':G,'rebased_common_ports':ports,'rebased_children':values,'scope':'Independent mixed-strong extension proof and census, plus full rebasing witness checks; source itself enumerates all seed-fixed extensions.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
