from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
for a,b in [('verification.json','verification_optimized.json'),('carrier.json','carrier_optimized.json')]:
    need((HERE/a).read_bytes()==(HERE/b).read_bytes(),('normal/O mismatch',a))
rank=json.loads((HERE/'carrier_rank.json').read_text())
literal=json.loads((HERE/'rank_literal.json').read_text())
need(rank['rank_mod101']==rank['rank_mod103']==literal['rank_mod101']==literal['rank_mod103']==571,'rank certificate')
lines=(HERE/'rank_contexts.txt').read_text().splitlines()
need(int(lines[0])==571,'bank size')
words=[]
for line in lines[1:]:
    v=list(map(int,line.split()));need(v[0]==len(v)-1 and v[0]<=5,'word length')
    need(all(0<=x<=6 for x in v[1:]),'word ports');words.append(v[1:])
need(words==rank['independent_words_mod101'],'exact discovery bank')
data=json.loads((HERE/'carrier_literal.json').read_text())
need(len(data['contexts'])==89 and data['graphs_counted']==2916352,'literal controls')
need(all(r['sum_chi1_H']==r['sum_chi2_H'] for r in data['contexts']),'mean identity')
sources=[
 '05-knowledge/results/primegame_port_spectrum_20261103/verify.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/port_scout.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/growth_identity.py',
 '05-knowledge/results/primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md',
 '05-knowledge/results/regular_parent_symmetry_20261029/verify.py',
 '05-knowledge/results/session145_endpoint_fence_ramsey_20261103/carrier/POOLED_PARITY_CARRIER.md',
 '05-knowledge/results/session144_jacobian_eth_mrb_20261103/jacobian/equivariant_composition_tower.md',
 '05-knowledge/results/session144_jacobian_eth_mrb_20261103/jacobian_followup/equivariant_gauge_hostile.md',
 '05-knowledge/results/jacobian_eth_prime_wait_20261102/WEIGHTED-CLASSIFICATION.md',
 '05-knowledge/results/session31_mirror_mrb_20261022/tree-moments/TREE_MOMENT_BRIDGE.md',
 '05-knowledge/results/overnight_rank_return_20261029/TWO_RUN.md',
 '05-knowledge/results/weighted_crt_clocks_20261028/CLOCKS.md',
 '05-knowledge/results/collatz_atom_route_memory_20261003.md',
 '05-knowledge/results/collatz_five_mirrors_20260929.md',
 '05-knowledge/results/uniform_prime_transport_20261103/README.md',
]
def pin(p):
    b=p.read_bytes().replace(b'\r\n',b'\n')
    return dict(bytes_lf=len(b),sha256_lf=hashlib.sha256(b).hexdigest())
artifacts={p.name:pin(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='manifest.json'}
deps={name:pin(ROOT/name) for name in sources}
out=dict(source_base='5e1cedf452c7cc1ec447205939a9d6cf4c0e8743',
         normalization='Text bytes with CRLF normalized to LF; portable across Git checkouts',
         artifacts=artifacts,dependencies=deps)
(HERE/'manifest.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf-8')
print(json.dumps(dict(artifacts=len(artifacts),dependencies=len(deps),rank=571,normal_optimized_identical=True)))
