# A unit shield with one exposed arc, and its exact module calculus

Let the old ordinal spine be x → y → Z → z, where x,y,z are single
vertices and Z is any **nonempty** tournament. Add two apices A,B.
The digit at a spine vertex is 2[A beats vertex]+[B beats vertex]:
the word is 12 [0 on the entire block Z] 2. The apices' mutual arc is
the only switch. Both resulting tournaments are strong. For A→B,
a vertex of Z reaches A and B directly, A reaches B and then x,
and B reaches A via x. From x or y or z the apices are reachable,
and B→x→y→Z→z reaches all spine vertices. For B→A, A reaches B
via z, and the other routes still work.

For k≥1 let P_k=P_Z(k) be the count of ordered k-slot path covers
of Z, allowing empties. Write h=H(Z)=P_1. Then

    H(A→B) = 4 P_3 + 9 P_2 − 3 h,
    H(B→A) = 4 P_3 + 9 P_2 − h.

In particular, for any t≥1 take Z=TT_t. Then P_k=k^t, and we
get the unbounded consecutive pair

    4·3^t+9·2^t−3,  4·3^t+9·2^t−1.

The first is 27,29. There is only one arc reversal, and the factor
multiplying +2 is exactly one. For a nontransitive Z the exact gap
is 2H(Z), not 2. This construction realizes a unit switch on
sparse rays, not at every absolute target.

## Proof by positive odd-cycle packs

All cycle packs not involving A or B contribute h. For each apex
cycle configuration, vertices in unused Z give a factor equal to
their own Hamilton count. Each chosen odd cycle has coefficient 2.
If several slots partition Z, weighting a partition by the product
of Hamilton counts is exactly the ordered-slot polynomial. In
particular a parity sign on one slot evaluates it by assigning
each odd cycle or singleton its slot sign. Thus the signed sum on
two slots with signs -,+ is zero (Z nonempty); on three slots with
signs -,+,+ it is h.

- One cycle through A only: its spine run starts y, then visits an
  odd-size subset of Z, giving P_2/2 by the two-slot identity.
- One cycle through B only: start x, choose a subset of Z, and choose
  y and/or z so the length of the spine run is even. Exactly two
  options for y,z fit each parity. This gives 2P_2.
- A and B in distinct cycles: assign Z to the two runs and the
  complement. The A run starts at y and uses an odd-size set in Z;
  the B run starts at x and has a forced optional z for parity.
  This gives (P_3−h)/2 by the signed three-slot identity.
- A and B in the same cycle without using the mutual arc: B's run
  starts x. If A's run starts y, partition Z into the two runs and
  the complement; choosing z in A's run is forced by total parity.
  This gives P_3. If A's run starts z, B uses a nonempty subset of
  Z and y is forced by parity. This gives P_2−h.
- A and B in the same cycle using A→B: B's run starts x. For each
  nonempty subset of Z, y is forced by parity; the empty subset
  permits x alone. Hence this contribution before the factor 2 is
  P_2. For B→A, A starts at y and each subset admits a unique choice
  of whether to end at z, again P_2. There is also the triangle
  ending at z with no y or Z, contributing h.

Assembling the single-cycle coefficients 2 and the double-cycle
coefficient 4 gives both formulas. Every item retains Z's full
slot counts rather than treating its size or factor h as enough
to predict the absolute baseline. Run `python check_shell.py`
for a literal Hamilton-path DP check in both orientations for
every labelled Z on 1–4 vertices.
