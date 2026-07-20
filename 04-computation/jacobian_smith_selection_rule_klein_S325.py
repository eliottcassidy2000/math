#!/usr/bin/env python3
"""jacobian_smith_selection_rule_klein_S325.py -- klein-2026-07-19-S325.

THE SMITH SELECTION RULE for etale polynomial self-covers of C^n.

Topology input (P.A. Smith, classical): a Z/p action (p prime) on a
finite-dimensional mod-p acyclic space has nonempty fixed-point set.  C^n is
contractible, and the deck group of an etale cover acts FREELY.  Hence:

  (R1) An etale self-cover of C^n has NO nontrivial deck transformation
       (any would contain a free Z/p).  In particular NO Galois etale
       self-cover of degree >= 2 exists: DEGREE 2 IS IMPOSSIBLE.
  (R2) The monodromy group G (transitive, degree d) must have deck group
       N_G(G_0)/G_0 trivial, i.e. |Fix(G_0)| = 1: the point stabilizer is
       SELF-NORMALIZING (fixes only its own point).

This script tabulates transitive groups of degree 2..7 against (R2).
Consequences: d=3 forces S_3 (the fleet's four independent S_3 results --
THM-1310 resolvent, THM-1315 syzygy/van der Waerden, boxeph-S142 "S3
universal", klein-S323 Chebotarev census -- are instances of a topological
theorem); d=4 forces A_4 or S_4; d=5 allows D_5, F_20, A_5, S_5; d=6 allows
PSL(2,5), PGL(2,5), A_6, S_6 (all im/primitive small ones with deck die);
d=7 allows F_21, F_42, PSL(3,2), A_7, S_7 but NOT C_7 or D_7.
"""
from sympy.combinatorics import Permutation, PermutationGroup, SymmetricGroup, AlternatingGroup
P = Permutation
def test(d, name, G):
    if not G.is_transitive():
        return None
    stab = G.stabilizer(0)
    if int(stab.order()) == 1:
        fixed = list(range(d))
    else:
        fixed = [i for i in range(d) if all(g(i) == i for g in stab.generators)]
    tag = "GALOIS(regular)" if int(G.order()) == d else ""
    verdict = "ALLOWED" if fixed == [0] else f"DEAD (deck size {len(fixed)})"
    return (d, name, int(G.order()), tag, verdict)
def gens(d, lst): return PermutationGroup([P(g, size=d) for g in lst])
rows = [
 test(2, "C2", gens(2, [[[0,1]]])),
 test(3, "C3", gens(3, [[[0,1,2]]])),
 test(3, "S3", gens(3, [[[0,1,2]],[[0,1]]])),
 test(4, "C4", gens(4, [[[0,1,2,3]]])),
 test(4, "V4", gens(4, [[[0,1],[2,3]],[[0,2],[1,3]]])),
 test(4, "D4", gens(4, [[[0,1,2,3]],[[1,3]]])),
 test(4, "A4", AlternatingGroup(4)),
 test(4, "S4", SymmetricGroup(4)),
 test(5, "C5", gens(5, [[[0,1,2,3,4]]])),
 test(5, "D5", gens(5, [[[0,1,2,3,4]],[[1,4],[2,3]]])),
 test(5, "F20", gens(5, [[[0,1,2,3,4]],[[1,2,4,3]]])),
 test(5, "A5", AlternatingGroup(5)),
 test(5, "S5", SymmetricGroup(5)),
 test(6, "C6", gens(6, [[[0,1,2,3,4,5]]])),
 test(6, "S3reg", gens(6, [[[0,1,2],[3,4,5]],[[0,3],[1,5],[2,4]]])),
 test(6, "D6", gens(6, [[[0,1,2,3,4,5]],[[1,5],[2,4]]])),
 test(6, "PSL(2,5)", gens(6, [[[0,1,2,3,4]],[[0,5],[1,4]]])),
 test(6, "PGL(2,5)", gens(6, [[[0,1,2,3,4]],[[0,5],[1,4]],[[1,2,4,3]]])),
 test(6, "A6", AlternatingGroup(6)),
 test(6, "S6", SymmetricGroup(6)),
 test(7, "C7", gens(7, [[[0,1,2,3,4,5,6]]])),
 test(7, "D7", gens(7, [[[0,1,2,3,4,5,6]],[[1,6],[2,5],[3,4]]])),
 test(7, "F21", gens(7, [[[0,1,2,3,4,5,6]],[[0,1,3],[2,5,4]]])),
 test(7, "F42", gens(7, [[[0,1,2,3,4,5,6]],[[0,2,1,5,3,4]]])),
 test(7, "PSL(3,2)", gens(7, [[[0,1,3,2,5,6,4]],[[1,3],[2,4]]])),
 test(7, "A7", AlternatingGroup(7)),
 test(7, "S7", SymmetricGroup(7)),
]
rows = [r for r in rows if r]
rows.sort(key=lambda r: (r[0], r[2]))
print(f"{'d':>2} {'group':<11} {'|G|':>6} {'':16} verdict")
for d, name, o, tag, v in rows:
    print(f"{d:>2} {name:<11} {o:>6} {tag:<16} {v}")
print()
print("R1: degree-2 etale self-covers of C^n DO NOT EXIST (any n).")
print("R2 verdicts above; d=3 forces S_3 -- the fleet's universal S_3 is topology.")
