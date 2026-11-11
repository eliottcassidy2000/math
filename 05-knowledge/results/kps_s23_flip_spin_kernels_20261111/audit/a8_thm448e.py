"""a8_thm448e.py (AUDIT) - isomorphism tests CD_n vs the THM-448 tower with networkx VF2 (third-party, independent
of the audit's own backtracking)."""
import sys, os, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import a8_thm448 as m
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
out = []
def P(s):
    print(s); out.append(s); sys.stdout.flush()
def dg(T):
    G = nx.DiGraph()
    G.add_nodes_from(range(len(T)))
    for u, s in enumerate(T):
        for v in s:
            G.add_edge(u, v)
    return G
for k in (3, 4, 5):
    Tc = m.cd_tournament(k)
    Tt = m.core(m.tower_matrix(k))
    Tv = m.converse(Tt)
    a = DiGraphMatcher(dg(Tc), dg(Tt)).is_isomorphic()
    b = DiGraphMatcher(dg(Tc), dg(Tv)).is_isomorphic()
    c = DiGraphMatcher(dg(Tt), dg(Tv)).is_isomorphic()
    d = DiGraphMatcher(dg(Tc), dg(m.converse(Tc))).is_isomorphic()
    P("n = %d: VF2: CD_n ~ tower: %s ; CD_n ~ converse of tower: %s ; tower self-converse: %s ; CD_n self-converse: %s" % (k, a, b, c, d))
open(os.path.join(m.HERE, "a8_thm448e.out"), "w").write("\n".join(out) + "\n")
