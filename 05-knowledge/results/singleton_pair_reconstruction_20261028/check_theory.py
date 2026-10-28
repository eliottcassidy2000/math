"""Reproduce the sharp obstruction, its matrix proofs, and hostile controls."""
from pathlib import Path
from itertools import combinations
import json
import reconstruct as r

HERE=Path(__file__).resolve().parent

def main():
    witness=r.witness_audit()
    a,b=(r.core.presentation(7,x) for x in (1072,4176))
    four=(a,b,r.core.complement(a),r.core.complement(b))
    r.need([r.oracle.hamilton(t) for t in four]==[127,127,93,93],'independent subset Hamilton DP')
    r.need([len(r.oracle.isomorphisms(t,t,True)) for t in four]==[1]*4,'full automorphism search')
    small=r.blue_dead_audit()
    bank=r.rooted_bank(7)
    r.need(bank==json.loads((HERE/'rooted-seven.json').read_text()),'entire seven-vertex bank')
    dead=sum(sum(r.mixed_ends(r.core.presentation(7,x)) and not r.blue(r.core.presentation(7,x))
                 for x in cell)>=2 for cell in bank['nontrivial'])
    controls=r.old.source('refined_controls','05-knowledge/results/root_ports_closed_carriers_20261028/check_carriers.py')
    hostile=controls.strong_controls()
    ta,tb=(controls.build_hostile(x) for x in (False,True))
    cyclic=lambda t:{v for v in combinations(range(len(t)),3)
                    if all(sum(r.core.edge(t,i,j) for j in v)==1 for i in v)}
    r.need(cyclic(ta)==cyclic(tb) and len(cyclic(ta))==67,'same labelled cyclic triples')
    r.need(r.oracle.hamilton(ta)==r.oracle.hamilton(tb)==462915,'81/27 Hamilton control')
    result=dict(status='FINITE-EXACT',witness=witness,blue_dead_audit=small,
                seven_vertex_summary={k:v for k,v in bank.items() if k!='nontrivial'},
                nonsingleton_cells_with_two_blue_dead_members=dead,hostile_controls=hostile,
                hostile_same_labelled_cyclic_triples=67,hostile_same_Hamilton_paths=462915)
    (HERE/'theory-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(witness=witness['masks'],seven_cells=bank['cells'],blue_dead_collision_cells=dead,
                         contexts=sum(witness['context_levels']),hostile_groups=[x['automorphisms'] for x in hostile[-2:]])))
    return result

if __name__=='__main__':main()
