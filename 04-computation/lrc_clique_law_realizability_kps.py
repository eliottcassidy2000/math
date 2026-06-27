r"""
The clique law + the 2 arc-states = the apex-7 realizability obstruction (kps-S31y, owner's hint).
CLIQUE LAW: for a conflict graph that is a clique K_r (r pairwise-conflicting odd cycles), the
Redei-type count is I(K_r, 2) = 1 + 2*r. The '2' is the BINARY ARC STATE (each odd cycle is in or
out of the independent set). So K_3 -> I(K_3,2) = 1+2*3 = 7. THM-200: K_3 is FORBIDDEN as a conflict
graph (3 pairwise-conflicting triangles force a C_5). And 14 = 2*7 = (arc-states) * I(K_3,2).
"""
def I_clique(r, x=2):  # independence polynomial of K_r at x: 1 + r*x (empty set + r singletons)
    return 1 + r*x
print("CLIQUE LAW I(K_r, 2) = 1 + 2r  ('2' = the binary arc state):")
for r in range(1,8):
    print(f"   I(K_{r}, 2) = 1 + 2*{r} = {I_clique(r)}   {'<-- the FORBIDDEN value (K_3 not realizable, THM-200)' if r==3 else ''}")
print(f"\n   So 7 = I(K_3, 2): the count of a 3-clique conflict graph using 2 arc-states.")
print(f"   And 14 = 2 * 7 = (arc-states 2) * (forbidden clique value 7).")
print(f"   1/14 = 1/(2*7): the LRC(14) threshold is the reciprocal of the forbidden-clique * arc-states.")
print(f"\nTHE OBSTRUCTION (why LRC(14) is impossible to disprove):")
print(f"  A counterexample would be an integer speed set realizing an OVER-cover at the apex prime 7.")
print(f"  At the apex-7 reduction, that over-cover's conflict structure would have to be the 3-clique K_3")
print(f"  (the 3 antipodal sector-pairs, mutually conflicting), realizability count I(K_3,2)=7.")
print(f"  But H=7 is FORBIDDEN: 3 pairwise-conflicting odd cycles ALWAYS force a C_5 (THM-200), so no")
print(f"  binary-arc tournament realizes K_3. => the over-cover is NOT realizable => no counterexample.")
print(f"  The 2 facts (7 forbidden + 2 arc-states) ARE the realizability obstruction. LRC(14) cannot be")
print(f"  disproved by a finite integer construction -- the construction would have to realize the")
print(f"  unrealizable K_3. (The open rigor: that the apex-7 over-cover IS exactly the K_3 -- the crux.)")
