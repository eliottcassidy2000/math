# HYP-2497 - Pollock's Sierpinski signal is a carry-pair lift, not a single dyadic obstruction

**Status:** OPEN proof program with exact finite scout through `10^6`.
**Source:** codex-2026-06-13.
**Extends:** HYP-2491, HYP-2490, HYP-2489, HYP-2455, HYP-2457, T815.
**Artifacts:** `04-computation/pollock_sierpinski_carry_scout_codex.py`, `05-knowledge/results/pollock_sierpinski_carry_scout_codex.out`, `07-reflections/pollock-sierpinski-carry-and-defect-pairs.md`.
**External anchors:** Pollock's tetrahedral conjecture remains open; MathWorld records the known `241` four-tetrahedra exceptions and largest `343867`: <https://mathworld.wolfram.com/PollocksConjecture.html>. Brady et al. survey generalized tetrahedral sums and the same frontier: <https://www.tandfonline.com/doi/full/10.1080/00029890.2021.1982635>.

## Claim

The Sierpinski/Waring intuition for Pollock's tetrahedral conjecture should be
lifted from single residues to carry-pair residues.

More precisely:

1. A pure dyadic local obstruction is probably the wrong proof object. In the
   scan, the single tetrahedral atom set `{Te_k mod 2^e}` is surjective for
   every `1 <= e <= 12`.
2. The dyadic signal reappears after imposing HYP-2491's lifted predicate:
   both `r` and `r+tri(k)` must lie in `D_4`. The observed tail pair classes
   stabilize while the dyadic pair universe grows like `4^e`.
3. The Waring/Sierpinski "just below the next atom" mechanism appears as a
   carry-window phenomenon for many four-defects, but not as a complete
   explanation of all defects.

Thus the live proof route is:

```text
Pollock = finite stencil certificate
          + no-long triangular self-correlation of D_4
          + dyadic/carry-pair lifts as pruning certificates.
```

This is a guardrail as much as a conjecture: do not hunt for a missing
single residue class modulo powers of two. Hunt for an incompatible lifted
pair/carry ledger.

## Computed Evidence

The scout reproduces the HYP-2491 frontier:

```text
limit = 1000000
tetrahedral_atoms = 180
D4_defects = 241
max_defect = 343867
five_term_misses = 0
triangular_defect_pairs = 601
last_pair_by_k = (825, 3142, 343867)
```

The triangular pair tail thins sharply:

```text
k >= 0:   601 pairs
k >= 50:  332 pairs
k >= 100: 168 pairs
k >= 300:  28 pairs
k >= 500:   1 pair
```

The last pair is still

```text
3142 -> 343867 = 3142 + tri(825).
```

## Dyadic Anti-Obstruction

For `1 <= e <= 12`, scanning `k < 4*2^e+16` gives every residue class modulo
`2^e` as a single tetrahedral number:

```text
e= 1: 2/2 atom residues
e= 2: 4/4 atom residues
...
e=12: 4096/4096 atom residues
```

This strongly suggests the lemma:

```text
For every e >= 1, {C(k+2,3) mod 2^e : k in Z} = Z/2^e Z.
```

If true, it explains why dyadic residue classes alone cannot prove Pollock.
The Sierpinski gasket is present in Lucas parity:

```text
Te_k = C(k+2,3) is odd iff k == 1 mod 4,
```

but that visible parity fractal does not create a missing tetrahedral residue.
One atom already fills the dyadic residue circle.

## Pair-Residue Compression

The lifted relation is different. For defect pairs with `k>=100`, the number
of observed pair classes modulo `2^e` stays bounded by `168` from `e=8`
through `e=12`, while the full pair universe has size `2^(2e)`.

```text
e= 3: tail pair classes  48 / 64,       compression  0.415
e= 4: tail pair classes 108 / 256,      compression  1.245
e= 5: tail pair classes 149 / 1024,     compression  2.781
e= 6: tail pair classes 163 / 4096,     compression  4.651
e= 7: tail pair classes 166 / 16384,    compression  6.625
e= 8: tail pair classes 168 / 65536,    compression  8.608
e= 9: tail pair classes 168 / 262144,   compression 10.608
e=10: tail pair classes 168 / 1048576,  compression 12.608
e=11: tail pair classes 168 / 4194304,  compression 14.608
e=12: tail pair classes 168 / 16777216, compression 16.608
```

This is not yet a proof, because it is only the observed finite defect set.
But it gives the right carrier: a hypothetical long pair must carry a coherent
dyadic address at both endpoints and at the triangular gap.

## Tournament Analysis

Vertices:

```text
dyadic levels e = 3,4,...,12
```

Pairwise observable:

```text
compression(e) = log2((2^e)^2 / observed tail pair classes mod 2^e)
```

Gauge:

```text
e -> f iff compression(e) > compression(f)
```

Tie Hamiltonian path:

```text
3 -> 4 -> ... -> 12
```

Result:

```text
scores: {3:0.415, 4:1.245, 5:2.781, 6:4.651, 7:6.625,
         8:8.608, 9:10.608, 10:12.608, 11:14.608, 12:16.608}
outdegrees: {3:0, 4:1, 5:2, 6:3, 7:4, 8:5, 9:6, 10:7, 11:8, 12:9}
directed_3_cycles: 0
SCC sizes: ten singleton SCCs
Hamiltonian paths: 1
champion path: 12 > 11 > 10 > 9 > 8 > 7 > 6 > 5 > 4 > 3
```

The tournament is fully transitive. Higher dyadic resolution dominates
because it separates the pair-address universe without producing new observed
tail pair classes after `e=8`.

## Carry Windows

The "just below a power" lower-bound mechanism familiar from Waring/Sierpinski
does have an analogue: many four-defects sit close below a tetrahedral number.
Let `gap_up(d)` be the distance from `d` to the next tetrahedral number.

```text
gap_up <=     1:   2 / 241
gap_up <=     5:  13 / 241
gap_up <=    10:  17 / 241
gap_up <=    50:  55 / 241
gap_up <=   100:  85 / 241
gap_up <=   500: 182 / 241
gap_up <=  1000: 218 / 241
gap_up <=  5000: 240 / 241
gap_up <= 10000: 241 / 241
```

The tail close-to-ceiling defects include:

```text
(37798, upper index 60, gap 22)
(39693, upper index 61, gap 18)
(62158, upper index 71, gap 38)
(64752, upper index 72, gap 72)
(134038, upper index 92, gap 6)
```

But the largest defect is not an ultra-near miss:

```text
343867 lies 5637 below Te_127.
```

So carry windows explain a real chunk of the exceptional geometry but cannot
be the whole proof by themselves.

## Proof Program

The next lemmas are ordered by leverage:

1. **2-adic surjectivity lemma.** Prove `{Te_k mod 2^e}` is all of `Z/2^eZ`
   for every `e`. This closes the false route and gives a clean Sierpinski
   anti-obstruction theorem.
2. **Pair-address necessary conditions.** For a hypothetical pair
   `r,r+tri(k) in D_4`, derive congruence/carry constraints on
   `(r mod 2^e, k mod 2^e)` that survive all four-sum convolution lifts.
3. **No-long-pair theorem.** Prove the HYP-2491 target:
   no such pair exists for `k>825`.
4. **Finite certificate.** Keep the width-3 shell stencil through `k<=825`
   as a compact verification object.

The creative analogy is:

```text
Sierpinski/Waring lower bound:
  a number just below the next kth power forces many summands.

Pollock tetrahedral lift:
  a number in a tetrahedral shell is hard only if two boundary totals,
  separated by the triangular shell gap, both have no four-atom convolution lift.
```

In the repo's boundary-lift language, single residues are scalar shadows.
The real Pollock obstruction is the lifted two-endpoint carry ledger.
