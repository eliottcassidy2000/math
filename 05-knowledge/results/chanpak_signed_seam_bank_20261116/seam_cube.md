# Two fixed strong cores: exact seam cancellation and a contact obstruction

Let the tape be TT_a → A → TT_c → B → TT_b, with outer guards a,b≥1, and let x and t be the head and sentinel from the [computational tape session](https://github.com/eliottcassidy2000/math/blob/6d38e8ceb0ccc8a061270e52507820bd4ecf53e6/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md). The head reverses its contacts against sets V in A and W in B. For core i in {1,2}, write h_i,q_i,E_i,L_i,J_i as in that session, and T_i,S_i for its Hamilton path counts ending and starting in the mask. Here z=1 if c=0, zero otherwise; x₀=2^a,y₀=2^b,g=2^c. The entire tournament is strong because every old vertex beats t, t beats the head, and the head reaches the tape via the left guard.

The exact count is
```
H = h1 h2 + g q1 q2 x₀ y₀ + x₀ E* - y₀ L* - J*,
E* = g q1 E2 + h2 E1,
L* = h1 L2 + g q2 L1,
J* = h2 J1 + h1 J2 + g E2 L1 + z T1 S2.
```
The preceding tape identity specialized to these two masks gives it directly. For the nonlocal term, with a predecessor in B and successor in A, the intermediate transitive vertices may be in either slot, so the obstruction is g E2 L1. With predecessor in A and successor in B, the intermediate vertices can belong to neither slot; only adjacency allows this seam contact, and gives T1 S2. Within a single exposed core the other core is traversed once, explaining h2 J1 and h1 J2. These are exact cancellations on one graph, not two separately scaled scalar claims.

The seeded October archive sampling (seed stream “seam-refill-16”, eight Markdown paths) brought up the inverse-trace-lattice [independent review](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/inverse_trace_lattice_20261020/review/README.md) and the [operation-fibre/Farey connections](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/adaptive_budgets_smooth_interfaces_20261020/CONNECTIONS.md). The latter gives a useful precise analogy: sum and product jointly decode an unordered pair but each alone loses distinct fibres. Here the count seeks a deliberately lossy projection on (E*,L*,J*); reconstructing or freezing their full joint state excludes valid repairs. This is an analogy about coordinate choice, not a claimed connection between the open conjectures.

## Exact two-core witness with both mask choices moving

Use strong cores of order five, lexicographic forward arc codes A=855 and B=221, no intervening vertices, left guard a=2, and freely grow b≥1. Keep every internal arc fixed. Change the masks as follows:

| Mask in A | Mask in B | E* | L* | J* | additive residual after the dominant term |
|---:|---:|---:|---:|---:|---:|
| 13 | 25 | 5444 | 5722 | 3442 | 18451 |
| 23 | 13 | 5465 | 5722 | 3524 | 18453 |

The two exact count functions are 29222·2^b+18451 and 29222·2^b+18453, so the difference is +2 at every height. The large constant contributions move together: 4E* increases by 84 and J* increases by 82. This flips five head contacts spread across the two cores; the individual core h factors are 9 and 13, so a naive independent +2 switch in either core would come out as +26 or +18. Here neither internal orientation nor ordinary core count changes and both masks genuinely move. First actual counts are 76895 and 76897.

A complementary two-clock witness uses A order 4 code 17 and B order 3 code 2, a=b+1, c=0 and B mask 1 unchanged. Moving A mask 10 to 9 changes (E*,L*,J*) from (132,159,50) to (138,171,48): both height coefficients change but 2E*−L*=105 is fixed, giving 768·4^b+105·2^b−35 to the same expression ending −33. First counts are 3247 and 3249.

## Why a triangle is a poor second register in one-clock design

For a fixed pair and growing b with a,c fixed, the projected slope is fixed precisely when L* is fixed. For a directed-triangle core, a mask of size s has L=3s and every other profile coordinate depends only on s (rotations act on all masks of each fixed size).

If A is the triangle, L*/3=L2+g q2 s. For B, the general bound 0≤L2≤q2−h2<q2 puts these mask-size bands strictly apart. If B is the triangle, a collision instead requires h1 Δs = −4g ΔL1. Since h1 is odd and |Δs|≤3, this is impossible unless Δs=0. Thus in a fixed slope fibre a triangle can never take an arithmetic step of its own in this schedule, even though a triangle can join an active fixed two-clock port or contribute to the background total. This is a structural obstruction for inserting a small second register, not a proof that every such fibre has width at most two.

## Replay and open issue

Run `python3 check_seam_cube.py`. It recomputes all slot/contact counts by subset DP and compares the resulting formula with a separate subset Hamilton path count of the complete output at several heights. Mask and code conventions agree with the linked tape session. The finite searches test bounded interfaces only; no infeasibility theorem for longer lanes follows.

## A three-odd lane that works independently on both clocks

A second one-core example sharpens the previous symmetric-clock three-odd lane: keep the strong order-nine core 53826135213 and move masks 53 → 77 → 141. All three share (h,q,E,L,T,S)=(1521,26412,10422,11997,631,750), while J decreases 3908→3906→3904. Hence they realize, for independently arbitrary a,b≥1,

`1521 + 26412·2^(a+b) + 10422·2^a − 11997·2^b − 3908 + 2j, j=0,1,2.`

The two batches use four and two head-contact reversals respectively; all internal arcs are identical. This does not have the signed movement of the two-core example, but it isolates more repair freedom in an exact slot-profile fibre and avoids any synchronization assumption between the high clocks. Appending another nontransitive block with fixed contacts would multiply this particular unit step by its path count, because all four mixed-interface readings are frozen.

## Same two cores and two independent clocks, with a seam correction

Both A and B can be the **same** strong order-five core, code 41, with both masks changing. No intervening transitive vertices are used. Their ordinary counts are h=13 and q=104. Exact coordinates are:

| mask | E | L | T | S | J |
|---:|---:|---:|---:|---:|---:|
| 9 | 33 | 33 | 5 | 5 | 4 |
| 3 | 41 | 32 | 6 | 5 | 6 |
| 24 | 32 | 41 | 5 | 6 | 6 |

Move masks (9,9) → (3,24) by four contact flips. The aggregate E*=L*=3861 remains the same, even though **both** local E and L readings change in **both** cores. The signed changes in J* are +52 from the independently scaled J terms, −65 from the back-crossing term E2 L1, and +11 from the adjacent seam term T1 S2. Their sum is −2. Thus the Hamilton count rises by exactly two for **independently arbitrary** outer guard lengths; both count functions are
```
169 + 10816·2^(a+b) + 3861·(2^a−2^b) − J*,   J*=1218,1216.
```
At a=b=1 they give 42215 and 42217. The adjacency term here is mathematically essential: the other changes on this fixed move sum to −13 and would not clear the h=13 divisibility. Inserting an unmarked gap gives a separate back-crossing coefficient while removing the adjacent correction, so this same programmed move then no longer has unit width. This is an exact example of why freezing or preparing the mixed separator independently of the core pair can destroy the repair.

Repository lineage: this h=13 seam mechanism was also seen in an earlier local sweep with another labeling/core encoding (code 955). The profile proof here does not assert a new isomorphism type.
