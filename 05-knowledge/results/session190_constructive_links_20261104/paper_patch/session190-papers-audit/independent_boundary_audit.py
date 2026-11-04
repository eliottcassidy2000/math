"""Independent edge-subset/union-find audit; standard library only.
Does not import or invoke boundary_signature_demo.py.
"""
from itertools import combinations
from pathlib import Path
from collections import Counter
import json

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / 'attached_papers' / 'boundary_signature_demo.json'
DATA = json.loads(INPUT.read_text())
T = tuple(DATA['terminals_in_cyclic_order'])

def collect(patch, deleted=None):
    n = patch['order']
    vertices = tuple(v for v in range(n) if v != deleted)
    terminals = tuple(t for t in T if t != deleted)
    internals = tuple(v for v in vertices if v not in terminals)
    edges = tuple(tuple(e) for e in patch['edges'] if deleted not in e)
    states = {}
    counts = Counter()
    # Each component contains a terminal. A nonempty forest has 1..|T| components.
    for edge_count in range(max(0,len(vertices)-len(terminals)), len(vertices)):
        for indices in combinations(range(len(edges)), edge_count):
            deg = [0] * n
            chosen = [edges[i] for i in indices]
            for a,b in chosen:
                deg[a] += 1
                deg[b] += 1
            if any(deg[v] != 2 for v in internals) or any(deg[t] > 2 for t in terminals):
                continue
            parent = list(range(n))
            def find(v):
                while parent[v] != v:
                    v = parent[v]
                return v
            cyclic = False
            for a,b in chosen:
                a,b = find(a),find(b)
                if a == b:
                    cyclic = True
                    break
                parent[a] = b
            if cyclic:
                continue
            endpoints = {}
            for t in terminals:
                if deg[t] == 1:
                    endpoints.setdefault(find(t), []).append(t)
            pairs = tuple(sorted(tuple(sorted(e)) for e in endpoints.values()))
            assert all(len(e) == 2 for e in pairs)
            signature = (tuple(deg[t] for t in terminals), pairs)
            states.setdefault(signature, chosen)
            counts[signature] += 1
    return states, counts

def serialized(states):
    return [{'degree_vector':s[0], 'endpoint_pairs':s[1], 'edges':states[s]} for s in sorted(states)]

all_states = {}
all_counts = {}
for name in ('H','H1'):
    for deleted in (None, *range(DATA[name]['order'])):
        all_states[name,deleted],all_counts[name,deleted] = collect(DATA[name],deleted)

A,B = all_states['H',None],all_states['H1',None]
assert A.keys() == B.keys()
assert len(A) == len(B) == 14
expected = set()
for isolated in range(4):
    for through in range(4):
        if isolated == through:
            continue
        degree = tuple(0 if i==isolated else 2 if i==through else 1 for i in range(4))
        pair = tuple(T[i] for i in range(4) if i not in (isolated,through))
        expected.add((degree,(pair,)))
expected.add(((1,1,1,1),((0,2),(4,6))))
expected.add(((1,1,1,1),((0,6),(2,4))))
assert set(A) == expected
old = {v:all_states['H',v] for v in range(DATA['H']['order']) if v not in T}
new = {v:all_states['H1',v] for v in range(DATA['H1']['order']) if v not in T}
coverage = {str(v):[w for w,s in old.items() if s.keys() <= q.keys()] for v,q in new.items()}
assert coverage == DATA['internal_deletion_coverage_new_to_old']
assert all(coverage.values())
terminal_checks = {str(t):all_states['H',t].keys() <= all_states['H1',t].keys() for t in T}
assert terminal_checks == DATA['terminal_deletion_inclusions']
assert all(terminal_checks.values())
assert {str(v):len(s) for v,s in old.items()} == DATA['old_deletion_signature_counts']
assert {str(v):len(s) for v,s in new.items()} == DATA['new_deletion_signature_counts']

# Check every stored witness directly, independent of enumeration representative choice.
def witness_signature(patch, witness):
    vertices = set(range(patch['order']))
    canonical_edges = {tuple(sorted(e)) for e in patch['edges']}
    chosen = {tuple(sorted(e)) for e in witness}
    assert len(chosen) == len(witness) and chosen <= canonical_edges
    adj = {v:set() for v in vertices}
    for a,b in chosen:
        adj[a].add(b); adj[b].add(a)
    assert all(len(adj[v])==2 for v in vertices-set(T))
    assert all(len(adj[v])<=2 for v in T)
    unseen = set(vertices)
    pairs = []
    while unseen:
        seed = min(unseen)
        component = {seed}; frontier = [seed]
        while frontier:
            u=frontier.pop()
            for v in adj[u]-component:
                component.add(v); frontier.append(v)
        unseen -= component
        if len(component)==1:
            assert seed in T and not adj[seed]
        else:
            endpoints = tuple(sorted(v for v in component if len(adj[v])==1))
            assert len(endpoints)==2
            pairs.append(endpoints)
    return tuple(len(adj[t]) for t in T),tuple(sorted(pairs))
for record in DATA['full_signature_witnesses']:
    s = tuple(record['degree_vector']),tuple(tuple(p) for p in record['endpoint_pairs'])
    for name in ('H','H1'):
        assert witness_signature(DATA[name], record[name+'_edges']) == s

# A two-anchor deletion certificate suffices: one old deletion per bipartition side.
anchor_map = {v:3 if v in (1,3,5,7,9) else 8 for v in new}
assert all(old[w].keys() <= new[v].keys() for v,w in anchor_map.items())
assert not any(all(w in ws for ws in coverage.values()) for w in old)

output = {
    'method':'Independent combinations of edge subsets plus disjoint-set cycle test; no imported demo code.',
    'input':str(INPUT),
    'full_state_counts':{'H':len(A),'H1':len(B)},
    'full_fragment_counts':{'H':sum(all_counts['H',None].values()),'H1':sum(all_counts['H1',None].values())},
    'exact_14_state_characterization':True,
    'stored_full_witnesses_validated':28,
    'terminal_deletion_state_counts':{str(t):{'H':len(all_states['H',t]),'H1':len(all_states['H1',t])} for t in T},
    'terminal_deletion_inclusions':terminal_checks,
    'internal_deletion_coverage':coverage,
    'minimal_old_anchor_count':2,
    'two_anchor_map':anchor_map,
    'complete_patch_states':{name:{'full':serialized(all_states[name,None]),'deletions':{str(v):serialized(all_states[name,v]) for v in range(DATA[name]['order'])}} for name in ('H','H1')}
}
(ROOT/'independent_boundary_audit.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({k:v for k,v in output.items() if k!='complete_patch_states'},indent=2))
