"""a8_thm448d.py (AUDIT) - CD_5 vs the THM-448 tower T_31 and its converse: vertex profiles (out- and in-triple
intersections), and H of the out-neighbourhood tournaments (15 vertices)."""
import sys, os, collections, itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import a8_thm448 as m
from aud_common import ham_c
out = []
def P(s):
    print(s); out.append(s); sys.stdout.flush()
def prof(T):
    n = len(T)
    res = collections.Counter()
    for v in range(n):
        c = collections.Counter()
        for a, b in itertools.combinations([u for u in range(n) if u != v], 2):
            c[len(T[v] & T[a] & T[b])] += 1
        res[tuple(sorted(c.items()))] += 1
    return sorted(res.items())
def link_H(T):
    """multiset of H(out-neighbourhood) and H(in-neighbourhood) over the vertices."""
    n = len(T)
    subs = []
    for v in range(n):
        O = sorted(T[v]); I = sorted(u for u in range(n) if u != v and u not in T[v])
        for S in (O, I):
            idx = {u: i for i, u in enumerate(S)}
            subs.append([set(idx[w] for w in T[u] if w in idx) for u in S])
    hs = ham_c(subs)
    return sorted(collections.Counter(hs[0::2]).items()), sorted(collections.Counter(hs[1::2]).items())
Tc = m.cd_tournament(5)
Tt = m.core(m.tower_matrix(5))
Tv = m.converse(Tt)
for name, T in (("CD_5", Tc), ("THM-448 T_31", Tt), ("converse of T_31", Tv)):
    lo, li = link_H(T)
    P("%s: H of out-neighbourhoods -> count %s ; H of in-neighbourhoods -> count %s" % (name, lo, li))
for name, T in (("CD_5", Tc), ("THM-448 T_31", Tt), ("converse of T_31", Tv)):
    p = prof(T)
    P("%s: %d vertex-profile classes, sizes %s" % (name, len(p), [c for _, c in p]))
    for pr, c in p:
        P("     %d vertices: %s" % (c, pr))
open(os.path.join(m.HERE, "a8_thm448d.out"), "w").write("\n".join(out) + "\n")
