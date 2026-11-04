#!/usr/bin/env python3
"""Reconcile independent census/filter audit against the producer's candidates."""
import csv, hashlib, json, pathlib
ROOT=pathlib.Path(__file__).resolve().parent
BASE=ROOT.parent
def need(condition,message):
    if not condition:
        raise RuntimeError(message)
source=BASE/'session189-c3-frontier/candidates.json'
producer=json.loads(source.read_text())
ind=[]
with (ROOT/'candidates.tsv').open() as f:
    for x in csv.DictReader(f,delimiter='\t'):
        ind.append(dict(node=int(x['node']),code=int(x['code']),
                        rows=list(map(int,x['rows'].split(','))),
                        starts=list(map(int,x['starts'].split(','))),
                        ends=list(map(int,x['ends'].split(',')))))
expected=[{**x,'code':int(x['code'])} for x in producer['candidates']]
need(len({x['node'] for x in expected})==len(expected),'Duplicate candidate node')
need(len({x['code'] for x in expected})==len(expected),'Duplicate candidate code')
need(ind==expected,'Candidate set, order, code, adjacency row, or endpoint array mismatch')
with (ROOT/'summary.tsv').open() as f:
    summary=[{k:int(v) for k,v in x.items()} for x in csv.DictReader(f,delimiter='\t')]
last=summary[-1]
need(last['n']==9,'Missing order 9 summary')
for field,key in [('merged','total_merged_nodes'),('prime_merged','prime_nodes'),('p_nonzero','p_nonzero'),('p_zero','p_zero')]:
    need(last[field]==producer[key],str((field,last[field],producer[key])))
need(len(ind)==last['p_zero'],'Candidate count mismatch')
need(all(x['labeled_mass']==x['expected_mass']==1<<(x['n']*(x['n']-1)//2) for x in summary),'Orbit mass mismatch')
need(all(x['ordinary']==2*x['merged']-x['self_converse'] for x in summary),'Merged/ordinary count mismatch')
files=[BASE/'session189-deletion-incoming/source.c',
       BASE/'session189-deletion-incoming/reproduced.out',source,
       BASE/'session189-c3-frontier/select.cpp']
files += [BASE/f'session189-deletion-incoming/data/tower_{n}.txt' for n in range(2,10)]
files += [ROOT/x for x in ('audit.cpp','summary.tsv','all_rows.tsv','candidates.tsv','run.out','run.log','runtime.json','verify.py','controls.cpp','controls.out','REVIEW.md','verify_compact.py')]
hashes={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
receipt=dict(status='PASS',scope='FINITE-EXACT, all tournament isomorphism/converse classes through order 9',
             independent_methods=['degree-only bucket permutation canonical form',
                                  'exact automorphism count as minimum multiplicity',
                                  'ordinary-class distinctness and full labeled-orbit mass',
                                  'pair-seeded distinguisher closure module test',
                                  'literal directed Hamiltonian path DFS on every prime representative'],
             candidate_rows_checked=len(ind),candidate_endpoint_values_checked=18*len(ind),
             all_prime_endpoint_values_checked=18*last['prime_merged'],
             first_candidate=ind[0],last_candidate=ind[-1],summary=summary,
             limitations=['This audit does not establish any C3 blow-up response residue exclusion or all-order lower bound by itself.',
                          'Original producer code was source-audited but not independently regenerated; exhaustive coverage instead follows from independently checked disjoint orbit mass.'],
             sha256=hashes)
(ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('status','candidate_rows_checked','candidate_endpoint_values_checked','all_prime_endpoint_values_checked')},indent=2))
print(json.dumps(last,indent=2))
