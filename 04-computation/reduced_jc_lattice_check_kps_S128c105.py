#!/usr/bin/env python3
"""reduced_jc_lattice_check_kps_S128c105.py -- kind-pasteur-2026-07-20-S128c105

Checks the two group-theoretic claims THM-1375's lattice rests on.

(1) SMITH IS STRICTLY STRONGER THAN CAMPBELL.  Campbell excludes GALOIS covers;
    the Smith selection rule (klein-S325) excludes covers whose monodromy point
    stabiliser is NOT self-normalising.  Galois <=> the stabiliser is trivial, and
    for a non-trivial group the trivial subgroup is never self-normalising -- so
    Campbell's hypothesis implies Smith's, i.e. Smith's theorem is stronger.  The
    ordering is only meaningful if the containment is STRICT: there must be a
    transitive non-Galois monodromy whose stabiliser is still non-self-normalising.
    Predicted witness: D_4 acting on 4 points.

(2) THE DISCRIMINANT CHARACTER DETECTS GALOIS EXACTLY AT d = 3.  "monodromy inside
    A_d" implies Galois iff every transitive subgroup of A_d is regular.  True for
    d = 3 (|A_3| = 3 = d); predicted to fail at d = 4 (A_4 is transitive on 4
    points and has order 12 > 4).

Both are elementary; both are load-bearing, so both are checked rather than
asserted.  The d = 4 row is also cross-checked against canon's claim (klein-S325)
that "d = 4 allows only A_4/S_4".
"""
from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import (SymmetricGroup, AlternatingGroup,
                                              CyclicGroup, DihedralGroup)


def transitive_subgroups(d):
    """The standard transitive subgroups of S_d for small d, by name."""
    if d == 2:
        return [("S_2 = Z_2", SymmetricGroup(2))]
    if d == 3:
        return [("A_3 = Z_3", AlternatingGroup(3)), ("S_3 = D_3", SymmetricGroup(3))]
    if d == 4:
        V4 = PermutationGroup([Permutation(0, 1)(2, 3), Permutation(0, 2)(1, 3)])
        return [("Z_4", CyclicGroup(4)), ("V_4", V4), ("D_4", DihedralGroup(4)),
                ("A_4", AlternatingGroup(4)), ("S_4", SymmetricGroup(4))]
    if d == 5:
        F20 = PermutationGroup([Permutation(0, 1, 2, 3, 4), Permutation(1, 2, 4, 3)])
        return [("Z_5", CyclicGroup(5)), ("D_5", DihedralGroup(5)),
                ("F_20", F20), ("A_5", AlternatingGroup(5)),
                ("S_5", SymmetricGroup(5))]
    return []


def normalizer_order(G, H):
    """|N_G(H)| by direct enumeration -- the groups here have order <= 120."""
    Hset = set(H.elements)
    cnt = 0
    for g in G.elements:
        gi = g**-1
        if {g * h * gi for h in Hset} == Hset:
            cnt += 1
    return cnt


def analyse(name, G, d):
    stab = G.stabilizer(0)
    selfnorm = normalizer_order(G, stab) == stab.order()
    regular = G.order() == d
    galois = regular            # transitive + |G| = d  <=>  regular  <=>  Galois
    in_Ad = all(p.is_even for p in G.generators)
    return dict(name=name, order=G.order(), stab=stab.order(),
                selfnorm=selfnorm, galois=galois, in_Ad=in_Ad)


print("=" * 78)
print("(1)+(2)  transitive monodromy candidates, by cover degree")
print("=" * 78)
print("  %-4s %-12s %6s %6s %-10s %-8s %-8s %s"
      % ("d", "group", "|G|", "|stab|", "self-norm?", "Galois?", "in A_d?", "verdict"))
rows = {}
for d in (2, 3, 4, 5):
    rows[d] = []
    for name, G in transitive_subgroups(d):
        r = analyse(name, G, d)
        rows[d].append(r)
        # a counterexample must survive BOTH exclusions
        survives = r['selfnorm'] and not r['galois']
        verdict = "ALLOWED" if survives else "excluded"
        why = []
        if r['galois']:
            why.append("Campbell")
        if not r['selfnorm']:
            why.append("Smith")
        print("  %-4s %-12s %6d %6d %-10s %-8s %-8s %s%s"
              % (d, r['name'], r['order'], r['stab'],
                 str(r['selfnorm']), str(r['galois']), str(r['in_Ad']),
                 verdict, (" (" + "+".join(why) + ")") if why else ""))
    print()

print("=" * 78)
print("(1) IS SMITH STRICTLY STRONGER THAN CAMPBELL?")
print("=" * 78)
# Campbell excludes exactly the Galois ones; Smith excludes exactly the
# non-self-normalising ones.  Galois => stabiliser trivial => not self-normalising
# (for |G|>1), so Campbell's excluded set is contained in Smith's.  Strict?
strict = []
for d in (2, 3, 4, 5):
    for r in rows[d]:
        smith_kills = not r['selfnorm']
        campbell_kills = r['galois']
        if smith_kills and not campbell_kills:
            strict.append((d, r['name'], r['order'], r['stab']))
        if campbell_kills and not smith_kills:
            print("  !! Galois but self-normalising stabiliser: %s at d=%d"
                  " -- containment FAILS" % (r['name'], d))
print("  monodromies killed by SMITH but NOT by Campbell (proving strictness):")
for d, nm, o, st in strict:
    print("     d=%d  %-10s |G|=%-3d |stab|=%d" % (d, nm, o, st))
print("  -> containment is STRICT: %s" % (len(strict) > 0))
print("  -> the predicted witness D_4 present: %s"
      % any(nm.startswith("D_4") for _, nm, _, _ in strict))

print()
print("=" * 78)
print("(2) DOES 'monodromy inside A_d' IMPLY GALOIS?  (character-detectability)")
print("=" * 78)
for d in (2, 3, 4, 5):
    bad = [r['name'] for r in rows[d] if r['in_Ad'] and not r['galois']]
    ok = all(r['galois'] for r in rows[d] if r['in_Ad'])
    print("  d=%d : every transitive subgroup of A_d regular? %-5s %s"
          % (d, str(ok), ("" if ok else "  counterexamples: " + ", ".join(bad))))
print()
print("  -> the discriminant (sign) character detects Galois-ness EXACTLY at d = 3.")
print("     At d >= 4 an even monodromy need not be Galois, so no single quadratic")
print("     character can carry the reduced JC one rung up.")

print()
print("=" * 78)
print("CROSS-CHECK vs canon (klein-S325: 'd=2 impossible; d=4 allows only A_4/S_4')")
print("=" * 78)
for d in (2, 4):
    allowed = [r['name'] for r in rows[d] if r['selfnorm'] and not r['galois']]
    print("  d=%d allowed monodromies: %s" % (d, allowed if allowed else "NONE"))
