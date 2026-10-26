# Joint two-port response histograms of strong atoms (rooted-lines R1): for a strong class A of order n,
#   F_A(h) = #presentations x of A with [tau x; first, last] = h   (presentations = Hamilton paths modulo Aut(A)),
#   B_A(h) = the blue ones (sigma x = x: the path reversal v_i -> v_(n-1-i) is an anti-automorphism of A).
# Classification: which strong atoms share F_A (exactly, as rooted target classes), share it after forgetting ports
# (target free classes), or share the one-sided data.  Also the converse pairing A <-> A^op (F_(A^op)(h) =
# F_A(rho h), rho[C;s,t] = [C^op;t,s]) and the self-converse atoms.  Usage: atom_responses.py n.
import sys, time
from collections import Counter, defaultdict
from rooted import tournaments, canon_rooted, ham_paths, partner, is_strong, relabel
from canon_nauty import canon_batch
n = int(sys.argv[1]); t0 = time.time()
classes = [a for a in tournaments(n) if is_strong(a)]
def aut_size(adj):
    # |Aut| = number of Hamilton paths mapping to the same presentation; count via canonical rooted forms of paths
    return None
resp = {}; blue = {}; tc = {}
for a in classes:
    P = ham_paths(a)
    tg = [canon_rooted(partner(a, p), p[0], p[-1]) for p in P]
    # |Aut(a)| = H / (#distinct presentations); presentations = paths modulo Aut: two paths are equivalent iff the
    # relabelling p -> 0..n-1 gives the same presentation adjacency
    pres = Counter()
    for p in P:
        perm = [0] * n
        for i, v in enumerate(p): perm[v] = i
        pres[relabel(a, perm)] += 1
    aut = len(P) // len(pres)
    assert all(c == aut for c in pres.values())
    tc[a] = len(pres)
    F = Counter(); Bc = Counter()
    for p, h in zip(P, tg):
        perm = [0] * n
        for i, v in enumerate(p): perm[v] = i
        x = relabel(a, perm)                         # the presentation
        # blue: sigma x = x  <=>  i -> n-1-i maps x to its converse
        conv = tuple(sum(1 << j for j in range(n) if (x[j] >> i) & 1) for i in range(n))
        sig = relabel(conv, [n - 1 - i for i in range(n)])
        F[h] += 1; Bc[h] += (sig == x)
    resp[a] = frozenset((h, F[h] // aut) for h in F)
    blue[a] = frozenset((h, Bc[h] // aut) for h in Bc if Bc[h])
free_t = {}
allh = sorted({h for a in classes for h, c in resp[a]})
fc = dict(zip(allh, canon_batch(allh)))
print(f"n = {n}: {len(classes)} strong atoms; presentations tc = H/|Aut| from {min(tc.values())} to {max(tc.values())}; "
      f"{len(allh)} distinct joint rooted targets  [{time.time() - t0:.0f}s]")
groups = defaultdict(list)
for a in classes: groups[resp[a]].append(a)
print(f"  distinct joint response histograms F_A: {len(groups)} (atoms sharing a histogram: {[len(g) for g in groups.values() if len(g) > 1]})")
coarse = defaultdict(list)
for a in classes: coarse[frozenset(Counter({fc[h]: c for h, c in resp[a]}).items())].append(a)
print(f"  after forgetting ports (target free classes with multiplicity): {len(coarse)} distinct")
withblue = defaultdict(list)
for a in classes: withblue[(resp[a], blue[a])].append(a)
print(f"  with the blue refinement B_A: {len(withblue)} distinct")
nblue = sum(1 for a in classes if blue[a])
print(f"  atoms with blue presentations: {nblue}")
# targets: strong or not
st = Counter(is_strong(h) for h in allh)
print(f"  joint targets strong / nonstrong: {st[True]} / {st[False]}")
