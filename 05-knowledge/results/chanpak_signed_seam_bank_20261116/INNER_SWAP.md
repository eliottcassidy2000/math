# A literal inner-spacer +2 from a factor swap

Keep the exact two-core seam defined in [seam_cube.md](seam_cube.md), but use a growing nonempty transitive middle spacer TT_t and fixed left/right guard lengths. Put X=2^a,Y=2^b. With t≥1 no Hamilton path can use a seam contact from A to B around the head without visiting the separating vertices, so the previous exact law factors as

```
H = (q1 X−L1)(q2 Y+E2)·2^t
    +h1h2 + h2(XE1−J1) − h1(YL2+J2).
```

Both factors and both local residuals can move simultaneously while preserving an exact step of two.

Set a=b=1. Fix strong A of order 7 (lexicographic code 2096861, h1=55,q1=860) and strong B of order 6 (code 149, h2=31,q2=328). Keep every internal arc fixed. Change their masks (110,61)→(86,53). Only four head contacts flip. The slope changes factors, not product:

```
(1107,880)=(27·41,20·44) → (1188,820)=(27·44,20·41).
```

Thus the primitive factors 41 and 44 change sides and the slope stays 974160. Simultaneously the left low coordinate XE1−J1 goes 846→548 and the signed right low coordinate −YL2−J2 goes −702→−534. The readout difference is

`31(−298)+55(168)=−9238+9240=2.`

The output pair is therefore **974160·2^t−10679, 974160·2^t−10677** at every t≥1. The first actual witnesses are 1,937,641 and 1,937,643. This gives a fixed-pair exact +2 repair with increasing distance between its active nontransitive cores: it does not need the earlier adjacency correction, and the exact slope is preserved without preserving either factor. A seeded finite search found this pair but no three-output lane; no impossibility claim follows.

Two relevant restrictions:
1. A constant difference on a fixed slope fibre is divisible by gcd(h1,h2), so the earlier adjacent 13/13 cores cannot give a flat unit for this growing-inner-spacer schedule.
2. If one also insists on equal slope for infinitely many synchronously expanding outer guards X=x·2^k,Y=y·2^k, a fixed pair of states must satisfy q1 x ΔE2=q2 y ΔL1 and Δ(L1E2)=0. The first equality forces same-direction movement of L1,E2; a nonzero such move strictly changes their nonnegative product. Thus both L1,E2 must individually freeze on this doubly growing schedule. The successful swap uses fixed outer guards.

`python3 check_inner_swap.py` independently derives the slot coordinates, checks all algebraic invariants, and counts Hamilton paths on both actual tournaments at t=1,2.
