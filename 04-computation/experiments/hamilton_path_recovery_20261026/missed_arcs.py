# Arcs on no Hamilton path inside strong tournaments: structure.  For each such arc u -> v test:
#   (a) a vertex w (w != u, v) with T - w non-strong in which u, v lie in strong components that are NOT consecutive
#       (u's component earlier: the arc then 'skips' a component of T - w);
#   (b) the same with w allowed to be any vertex, and which w (degree, position).
import sys
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import tournaments, ham_paths, is_strong
from ham_recovery import components
def sub(adj, keep):
    idx = {v: i for i, v in enumerate(keep)}
    return tuple(sum(1 << idx[u] for u in keep if (adj[v] >> u) & 1) for v in keep), idx
nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
for n in range(3, nmax + 1):
    cnt = 0; ok_a = 0; examples = []
    for adj in tournaments(n):
        if not is_strong(adj): continue
        HP = ham_paths(adj)
        used = {(a, b) for P in HP for a, b in zip(P, P[1:])}
        for u in range(n):
            for v in range(n):
                if not ((adj[u] >> v) & 1) or (u, v) in used: continue
                cnt += 1
                witnesses = []
                for w in range(n):
                    if w in (u, v): continue
                    keep = [x for x in range(n) if x != w]
                    s, idx = sub(adj, keep)
                    if is_strong(s): continue
                    comps = components(s)
                    ci = {keep[x]: i for i, c in enumerate(comps) for x in c}
                    if ci[v] >= ci[u] + 2: witnesses.append((w, len(comps)))
                if witnesses: ok_a += 1
                if len(examples) < 2: examples.append((adj, (u, v), witnesses, [bin(r).count('1') for r in adj]))
    print(f"n = {n}: {cnt} missed arcs inside strong tournaments; with a skipping cut vertex w: {ok_a}")
    for adj, arc, wit, sc in examples:
        print(f"    rows {adj}, scores {sc}, missed arc {arc[0]}->{arc[1]}, cut vertices w (#components of T-w): {wit}")
