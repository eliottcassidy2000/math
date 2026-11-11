"""Remove only genuine marked automorphism duplicates from the seed census."""
from pathlib import Path
from itertools import permutations,product,combinations
from collections import defaultdict
from functools import lru_cache
import json
import port_transport as pt

HERE=Path(__file__).resolve().parent


@lru_cache(None)
def automorphisms(rows):
    cells=defaultdict(list)
    for v,row in enumerate(rows):cells[row.bit_count()].append(v)
    groups=list(cells.values());out=[];n=len(rows)
    for choices in product(*(list(permutations(cell)) for cell in groups)):
        p=list(range(n))
        for cell,images in zip(groups,choices):
            for v,w in zip(cell,images):p[v]=w
        if all(bool(rows[v]>>w&1)==bool(rows[p[v]]>>p[w]&1) for v in range(n) for w in range(v)):
            out.append(tuple(p))
    return tuple(out)


def move(mask,p):return sum(1<<p[v] for v in range(len(p)) if mask>>v&1)


def canonical_mark(mask,ports,aut):
    return min((move(mask,p),tuple(sorted(p[v] for v in ports))) for p in aut)


def main():
    birth=json.loads((HERE/'birth_atoms.json').read_text());cones=json.loads((HERE/'cones.json').read_text())
    state_sets=defaultdict(set);atom_sets=defaultdict(set);reach_sets=defaultdict(set)
    core_seeds=defaultdict(set);stable_sets=defaultdict(set)
    for item in birth['states']:
        n=item['n'];code=item['code'];rows=tuple(pt.ib.cc.decode(n,code));aut=automorphisms(rows)
        key=(code,canonical_mark(item['selected'],(),aut)[0])
        state_sets[n].add(key)
        if not item['predecessors']:atom_sets[n].add(key)
        if item['reachable_from_singleton']:reach_sets[n].add(key)
        if len(pt.cp.components(rows))==1:core_seeds[n].add(key)
    for item in cones['records']:
        n=item['n'];aut=automorphisms(tuple(item['rows']))
        for r in range(1,len(item['ports'])+1):
            for ports in combinations(item['ports'],r):
                stable_sets[n,r].add((item['code'],canonical_mark(item['selected'],ports,aut)))
    levels=[]
    for n in range(1,8):
        levels.append(dict(n=n,marked_seed_orbits=len(state_sets[n]),marked_deletion_atom_orbits=len(atom_sets[n]),
                           singleton_reachable_orbits=len(reach_sets[n]),strong_marked_seed_orbits=len(core_seeds[n]),
                           strong_stable_subsets=[len(stable_sets[n,r]) for r in range(1,6)]))
    # No relabelling may alter selected-set status, predecessor existence,
    # or reachability. Verify all equivalent entries agree explicitly.
    fibres=defaultdict(set)
    for item in birth['states']:
        rows=tuple(pt.ib.cc.decode(item['n'],item['code']))
        key=(item['n'],item['code'],canonical_mark(item['selected'],(),automorphisms(rows))[0])
        fibres[key].add((item['U'],bool(item['predecessors']),item['reachable_from_singleton']))
    pt.ib.tc.require(all(len(values)==1 for values in fibres.values()),'quotient invariants')
    report=dict(levels=levels,scope='marked isomorphism classes, not just all marks on ordinary representatives',
                stable_core_subset_totals=[sum(len(stable_sets[n,r]) for n in range(1,8)) for r in range(1,6)])
    (HERE/'marked_orbits.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
