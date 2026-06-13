#!/usr/bin/env python3
"""
ck_class_ordering.py -- kind-pasteur-2026-03-13-S60

Check: do the 4 H-classes at p=11 have CONSISTENT c_k ordering across all k?
I.e., does Class A > Class B for ALL k?

Previous data:
H=95095: c3=55, c5=594, c7=3960, c9=11055, c11=5505  (N=21169) [CLASS A]
H=93467: c3=55, c5=550, c7=3586, c9=10197, c11=5153  (N=19541) [CLASS B]
H=93027: c3=55, c5=484, c7=3399, c9=9350,  c11=5109  (N=18397) [CLASS C]
H=92411: c3=55, c5=572, c7=3729, c9=10274, c11=4999  (N=19629) [CLASS D]

Ordering by c5: A > D > B > C
Ordering by c7: A > D > B > C
Ordering by c9: A > D > B > C
Ordering by c11: A > B > C > D  <-- DIFFERENT!
Ordering by N: A > D > B > C
Ordering by H: A > B > C > D

Wait -- B > C in c_k for k=5,7,9 but C > D in c11!
And D > B in c5,c7,c9 but B > D in H!

This is EXACTLY the trade-off: D has more cycles of lengths 5,7,9 but
FEWER Hamiltonian cycles (c11) than B. And c11 contributes to alpha_2
in a way that HELPS B's H.

Let me compute: how much does each c_k contribute to H via the OCF?
"""

# The classes (from previous computation):
classes = {
    'A': {'c': {3: 55, 5: 594, 7: 3960, 9: 11055, 11: 5505}, 'H': 95095,
           'N': 21169, 'a2': 10879, 'a3': 1155},
    'B': {'c': {3: 55, 5: 550, 7: 3586, 9: 10197, 11: 5153}, 'H': 93467,
           'N': 19541, 'a2': 11220, 'a3': 1188},
    'C': {'c': {3: 55, 5: 484, 7: 3399, 9: 9350, 11: 5109}, 'H': 93027,
           'N': 18397, 'a2': 11110, 'a3': 1474},
    'D': {'c': {3: 55, 5: 572, 7: 3729, 9: 10274, 11: 4999}, 'H': 92411,
           'N': 19629, 'a2': 10912, 'a3': 1188},
}

print("PER-LENGTH ORDERING:")
for k in [3, 5, 7, 9, 11]:
    ranking = sorted(classes.keys(), key=lambda x: classes[x]['c'][k], reverse=True)
    vals = [f"{cl}:{classes[cl]['c'][k]}" for cl in ranking]
    print(f"  c_{k}: {' > '.join(vals)}")

print(f"\n  N:  ", end="")
ranking = sorted(classes.keys(), key=lambda x: classes[x]['N'], reverse=True)
print(" > ".join(f"{cl}:{classes[cl]['N']}" for cl in ranking))

print(f"  H:  ", end="")
ranking = sorted(classes.keys(), key=lambda x: classes[x]['H'], reverse=True)
print(" > ".join(f"{cl}:{classes[cl]['H']}" for cl in ranking))

print(f"  a2: ", end="")
ranking = sorted(classes.keys(), key=lambda x: classes[x]['a2'], reverse=True)
print(" > ".join(f"{cl}:{classes[cl]['a2']}" for cl in ranking))

print(f"  a3: ", end="")
ranking = sorted(classes.keys(), key=lambda x: classes[x]['a3'], reverse=True)
print(" > ".join(f"{cl}:{classes[cl]['a3']}" for cl in ranking))

# The KEY insight: D has more cycles than B (N_D=19629 > N_B=19541)
# but D has LOWER H than B (H_D=92411 < H_B=93467).
# The ONLY reason: the alpha terms.
# 2*N: D gets +2*(19629-19541) = +176
# 4*a2: D gets 4*(10912-11220) = -1232
# 8*a3: D gets 8*(1188-1188) = 0
# Net: +176 - 1232 = -1056 = H_D - H_B = 92411 - 93467 = -1056 CHECK

print(f"\n\nDETAILED B vs D COMPARISON:")
dN = classes['D']['N'] - classes['B']['N']
da2 = classes['D']['a2'] - classes['B']['a2']
da3 = classes['D']['a3'] - classes['B']['a3']
dH = classes['D']['H'] - classes['B']['H']
print(f"  dN  = {dN:>+6},  2*dN  = {2*dN:>+6}")
print(f"  da2 = {da2:>+6},  4*da2 = {4*da2:>+6}")
print(f"  da3 = {da3:>+6},  8*da3 = {8*da3:>+6}")
print(f"  Total dH = {2*dN + 4*da2 + 8*da3} (actual: {dH})")
print(f"\n  D has MORE cycles (+{dN}) but FEWER disjoint pairs ({da2})")
print(f"  The disjoint pair DEFICIT costs 4*({abs(da2)}) = {4*abs(da2)}")
print(f"  vs the cycle SURPLUS benefit 2*{dN} = {2*dN}")
print(f"  Net: cycle surplus is OVERWHELMED by disjoint pair change")
print(f"  This is the ONLY case at p=11 where more cycles doesn't help!")

# Why? D has c11=4999 vs B's c11=5153. The Hamiltonian cycles are
# the MOST overlapping (every pair of Ham cycles shares ALL vertices),
# so they contribute 0 to alpha_2. Losing Ham cycles REDUCES N
# without changing alpha_2, so the "free" N from Ham cycles matters.

# Actually wait: D has MORE total N but fewer Ham cycles.
# D's extra cycles come from c5, c7, c9 (not c11).
# These shorter cycles OVERLAP more with other cycles?
# No: shorter cycles overlap LESS (smaller vertex sets).

# Let me compute the PER-LENGTH contribution to alpha_2
# alpha_2 = sum of disj(k1, k2) where disj counts disjoint pairs
# with weight n(V1)*n(V2)

# From earlier data:
disj_B = {(3,3): 517, (3,5): 3740, (3,7): 4862, (5,5): 2101}
disj_D = {(3,3): 506, (3,5): 3652, (3,7): 4840, (5,5): 1914}

print(f"\n\nPER-LENGTH DISJOINT PAIR COMPARISON (B vs D):")
for kpair in [(3,3), (3,5), (3,7), (5,5)]:
    db = disj_B.get(kpair, 0)
    dd = disj_D.get(kpair, 0)
    diff = dd - db
    print(f"  disj({kpair[0]},{kpair[1]}): B={db}, D={dd}, D-B={diff:>+5}")

print(f"\n  Total a2: B={sum(disj_B.values())}, D={sum(disj_D.values())}, "
      f"D-B={sum(disj_D.values())-sum(disj_B.values()):>+5}")

# Interesting: D has FEWER disjoint pairs for (3,3), (3,5), (3,7), AND (5,5)
# So D has fewer disjoint pairs overall. But 4*da2 = 4*(-308) = -1232
# which is a BIG reduction.
# Wait: da2 = 10912 - 11220 = -308. So D has FEWER disjoint pairs.
# 4*da2 = 4*(-308) = -1232. This REDUCES H_D relative to H_B.
# And 2*dN = 2*88 = +176.
# So the fewer disjoint pairs HURT more than the extra cycles help.

# BUT WAIT: fewer disjoint pairs means MORE overlaps.
# More overlaps with fewer disjoint pairs should mean LESS alpha_2.
# Lower alpha_2 means LOWER H (since H = 1+2N+4*a2+8*a3).
# So D's lower alpha_2 HURTS it, even though it has more N!

# This is COUNTER-INTUITIVE to the earlier "Paley minimizes alpha_2" narrative.
# Let me re-examine...

# At the TOP: A (Paley) has alpha_2 = 10879 (LOWEST) and N=21169 (HIGHEST).
# A's low alpha_2 "hurts" its H... but its HUGE N advantage compensates.
# For B vs D: D has lower alpha_2 AND higher N, yet LOWER H.
# This seems contradictory!

# NO: Let me recheck the actual H computation.
# H_B = 1 + 2*19541 + 4*11220 + 8*1188 = 1 + 39082 + 44880 + 9504 = 93467
# H_D = 1 + 2*19629 + 4*10912 + 8*1188 = 1 + 39258 + 43648 + 9504 = 92411
# Delta: +176 - 1232 + 0 = -1056. Correct.

# So B has HIGHER alpha_2 (11220 vs 10912) AND LOWER N (19541 vs 19629).
# The HIGHER alpha_2 contributes +4*(11220-10912) = +1232 to B's H.
# The LOWER N costs -2*(19541-19629) = -176.
# Net: B gains 1232-176 = 1056 over D.

# KEY INSIGHT: alpha_2 (disjoint pairs) HELPS H, not hurts it!
# H = 1 + 2N + 4*alpha_2 + 8*alpha_3
# Each term is ADDED. So MORE disjoint pairs = HIGHER H!

# Then WHY does Paley (with LOWEST alpha_2) maximize H?
# Because Paley's N advantage (+1628 over B) gives +3256, which overcomes
# alpha_2 disadvantage of 4*(10879-11220) = -1364.
# Net: 3256 - 1364 = 1892... no, +3256-1364 + (8*(1155-1188)) = 3256-1364-264 = 1628. CHECK.

print(f"\n\n{'='*70}")
print(f"  REVISED UNDERSTANDING OF H-MAXIMIZATION")
print(f"{'='*70}")
print(f"  H = 1 + 2*N + 4*alpha_2 + 8*alpha_3")
print(f"  ALL terms are POSITIVE. More of anything helps.")
print(f"")
print(f"  But N and alpha_2 are ANTI-CORRELATED:")
print(f"  More cycles => more overlaps => fewer disjoint pairs => lower alpha_2")
print(f"")
print(f"  The question: does the 2*dN gain from more cycles")
print(f"  overcome the 4*|d(alpha_2)| loss from fewer disjoint pairs?")
print(f"")
print(f"  ANSWER: For Paley (A) vs others, YES (always).")
print(f"  But for D vs B (within the non-Paley classes), NO.")
print(f"  D has more cycles but it doesn't compensate.")
print(f"")
print(f"  Why does Paley's N advantage SCALE better?")
print(f"  Paley's extra cycles come overwhelmingly from c9 (+858 over B)")
print(f"  and c7 (+374 over B), which are LARGE cycles.")
print(f"  Large cycles overlap with EVERYTHING, contributing 0 to alpha_2.")
print(f"  So Paley's N increase is 'cheap' — it doesn't reduce alpha_2 much.")

# Let's verify: what fraction of Paley's N advantage comes from each k?
print(f"\n  Paley's N advantage by cycle length:")
for k in [3, 5, 7, 9, 11]:
    dc = classes['A']['c'][k] - classes['B']['c'][k]
    print(f"    dc_{k} = {dc:>+5} ({100*dc/(classes['A']['N']-classes['B']['N']):>+6.1f}% of dN)")

print(f"\n  D's N advantage over B by cycle length:")
for k in [3, 5, 7, 9, 11]:
    dc = classes['D']['c'][k] - classes['B']['c'][k]
    print(f"    dc_{k} = {dc:>+5} ({100*dc/(classes['D']['N']-classes['B']['N']):>+6.1f}% of dN)")

# D's extra cycles come from c5(+22), c7(+143), c9(+77), but it LOSES c11(-154).
# The c11 loss is significant: Hamiltonian cycles overlap with everything,
# so losing them reduces N without affecting alpha_2.
# The c5/c7 gains create MORE disjoint pairs (smaller cycles can be disjoint),
# but the alpha_2 STILL drops because the total effect is complex.

print(f"\n\nDONE.")
