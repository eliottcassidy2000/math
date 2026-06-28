# LRC14 Erdos-870 Tournament-4 Canary/Filler Lens

Anchors: HYP-3149, HYP-3148, HYP-3147, HYP-3146, HYP-3145, HYP-3144, HYP-3143, HYP-3141, HYP-3140, HYP-3138, HYP-3137, HYP-3134,
HYP-3133, HYP-3124, HYP-3118, HYP-3116, HYP-3093, HYP-3097, T1214,
LTI-275, LTT-173, OPEN-Q-108.

The prompt's two n=4 tournament tables turned out to be more than a naming
scheme.  HYP-3143 captures the exact-order subbasis lesson, HYP-3144 captures
the adjacent pair-function scalarization alarm, HYP-3145 captures the broader
filler-core interface, upstream HYP-3146 captures the cover-versus-scaffold
shift-package policy, upstream HYP-3147 captures the n=3 edge-flip /
Worpitzky kernel, and upstream HYP-3148 captures the live-core deletability
audit.  This note is the fixed-path canary-slice refinement of
that policy.  Together they form a tiny exact laboratory for the interface
pattern in the formalized Erdos-870 solution: a low-order source becomes
proof-useful only after deterministic filler coordinates and canary deletion
rules are named.

In the fixed-Hamiltonian-path model, the free coordinates are

```text
a=(0,2), b=(1,3), c=(0,3)
```

over the path `0->1->2->3`.  The class map has a large `S` fiber:

```text
T={E}, +={a}, -={b}, S={c,ab,ac,bc,abc}.
```

This explains why the displayed `E,a,b,c` table feels group-like but is not a
class group.  It is the pairwise shadow of a quotient map that has already
destroyed a coordinate.

The second model fixes `c` unflipped.  That gives partial score sequence
`(0,1,1,2)` and leaves the two-source table on `x=a,y=b`:

```text
E -> T, x -> +, y -> -, xy -> S.
```

So `c` is the exact canary/filler arc.  With `c=0`, the two-coordinate source
is a perfect four-class transversal.  With `c=1`, all four completions collapse
to `S`.  That is the whole lesson: the missing coordinate is not decorative;
it determines whether the quotient is faithful or a bulk collision.

The LRC14 transfer is to stop asking for a raw group law on edge or tile
classes.  The proof-facing object should be

```text
fixed_path_word
+ c_canary_status
+ xy_completion_table
+ S_bulk_fiber_words
+ deletion/restoration_sidecar
+ edge_tip_tail_exit_or_named_debt.
```

That object is small enough to formalize.  It is also exactly the kind of
packet HYP-3141 wants: edge witnesses should carry tail payload, tip payload,
observer orbit, commutator defect, and terminal exit.  HYP-3149 supplies the
smallest possible worked example of why.

The creative next theorem would be a finite packet lemma:

```text
Every legal local edge quotient that looks like a two-coordinate source must
either identify a canary/filler slice making the completion table faithful, or
emit its collision fiber as named restoration debt.
```

This would not prove LRC14 by itself, but it improves the thought process
around HYP-3141/HYP-3140/HYP-3138: two coordinates can be enough only when the
third coordinate has been fixed, audited, and made visible as a canary.  The
Erdos-870 connection is precisely that formal solutions can split a problem
into a low-order source plus finite deterministic fillers without pretending
the fillers are harmless.

Challenged assumption: vertices need not be tournaments, runners, or raw arcs.
For this problem the useful tournament vertices are proof interfaces:
`fixed_c_xy_transversal`, `erdos870_order2_plus_filler_interface`,
`tiling_hamiltonian_path_cube`, `edge_tip_tail_information_packet`,
`S_bulk_collision_fiber`, `fiber_pgf_or_distribution_sidecar`,
`raw_score_sequence_scalar`, and `raw_einheit_group_table_numerology`.
