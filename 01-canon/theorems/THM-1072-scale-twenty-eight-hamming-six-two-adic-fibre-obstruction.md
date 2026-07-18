---
id: THM-1072
title: Scale-twenty-eight Hamming-six two-adic fibre obstruction
status: PROVED STRUCTURAL + DUAL INDEPENDENT FINITE-EXACT — the terminal Z/4-fibre deficit is reproduced by the frozen NumPy primary, a separately structured literal-CRT Python flag certificate, and an independent standard-library C++ referee; all agree on the 3,170 scalar survivors, 19,020 owner obligations, threshold census, and 6,628,500 exact reachable-mask incidences
source: codex-2026-07-17-S66 scale-twenty-eight continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-983, THM-986, THM-988, THM-989, THM-990, THM-992, THM-994]
related: [THM-963, THM-969, THM-981, HYP-6820]
verification:
  - 04-computation/lrc13_scale_twenty_eight_hamming_six_two_adic_fibre_obstruction_codex_c28.py
  - 05-knowledge/results/lrc13_scale_twenty_eight_hamming_six_two_adic_fibre_obstruction_codex_c28.out
  - 04-computation/lrc13_scale_twenty_eight_hamming_six_mod_four_flag_codex_c28.py
  - 05-knowledge/results/lrc13_scale_twenty_eight_hamming_six_mod_four_flag_codex_c28.out
  - 04-computation/lrc13_scale_twenty_eight_hamming_six_two_adic_fibre_referee_codex_c28.cpp
  - 05-knowledge/results/lrc13_scale_twenty_eight_hamming_six_two_adic_fibre_referee_codex_c28.out
---

# THM-1072 — scale twenty-eight has a terminal two-adic fibre deficit

The theorem is:

> **The primitive proper AP-centred common-scale-28 Hamming-six face is
> empty.**

The proof-facing certificate retains only the order-two and order-four thick
fibres of `Z/28 -> Z/4`.  It is already terminal.  The competing
`Z/28 -> Z/7` flag is not terminal, and retaining both flags in the full
`4 x 7` CRT grid does not improve a single terminal bound on the scalar
survivor bank.  Thus this face does **not** need a finer residual carrier.

The original frozen primary checks every labelled scalar row; multiplicative
orbits are telemetry only.  Two independently structured replays now agree
with it, so the theorem is promoted to `PROVED STRUCTURAL + DUAL INDEPENDENT
FINITE-EXACT`.

## 1. Hereditary divisor grammar

At common scale `c=28=4*7`, the effective orders and literal unit counts are

```text
D                 1  2  4  7  14  28
# units            1  1  2  6   6  12.
```

For an order word `(D_1,...,D_6)`, hereditary leave-one-out lcm means

```text
lcm(D_i : i != j) = 28       for every j.
```

This has an exact valuation grammar:

```text
at least two coordinates have nu_2(D_i)=2, and
at least two coordinates have nu_7(D_i)=1.                    (1)
```

Indeed, deletion preserves the full 2-adic maximum exactly when at least two
coordinates provide the factor four, and similarly for seven.  This proves
both directions of (1), without enumeration.

Inclusion-exclusion supplies an independent count.  Unweighted order words
give

```text
6^6 - 16384 - 5103 + 1792 = 26961.                            (2)
```

Weighting each order by its number of units gives

```text
28^6 - 52706752 - 151552 + 16576 = 429048576                 (3)
```

literal state words per support.  Across the `binom(12,6)=924` labelled
supports this represents

```text
924 * 429048576 = 396440884224
```

raw labelled states.  The primary checks the lcm predicate against (1) on all
`6^6` order words and checks both inclusion-exclusion identities.

## 2. Literal masks and the exact scalar layer

For provider label `a`, owner label `b`, effective order `D`, and unit `u`, let
`B` be the unique CRT base

```text
B = D*a (mod 13),       B = u (mod D).
```

The local sheet mask is

```text
M(a,D,u;b)
 = {t in Z/28 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]},       (4)
```

where brackets denote the centred residue.  Its size is independent of `u`.
Writing `r=a/b in F_13^*`, literal interval counting gives

```text
D=1 : 28 at r=1, otherwise 0;
D=2 : 14 on {1,6,7}, otherwise 0;
D=4 :  7 on {1,3,4,6,7,9,10}, otherwise 0;
D=7 :  8 at r=1, otherwise 4;
D=14:  6 at r=1, otherwise 4;
D=28:  5 on {1,6,7}, otherwise 4.                            (5)
```

More generally the count behind (5) is

```text
|M(a,D,u;b)|
 = (28/D) * #{x in (-D,D] : x = D*r (mod 13)}.                (6)
```

The primary constructs every CRT base both algebraically and by literal
search, constructs every mask from (4), and checks (5)-(6) mask by mask.

Summing (5) over providers is a necessary scalar capacity test.  On all

```text
924 * 26961 = 24911964
```

support/order contexts, the number of owners reaching scalar capacity `28`
has histogram

```text
0:120024, 1:5260824, 2:10675332, 3:6969052,
4:1724910, 5:158652, 6:3170.                                  (7)
```

Thus only `3170` contexts, on `206` supports, can possibly work at all six
owners.  None contains order one.  They stand for `9275904` literal unit
words and have `1480` distinct capacity vectors.  All other contexts are
already impossible by the union bound at one or more owners.

## 3. The sound anchor/nonanchor relaxation

We use the general inequality proved in THM-994.  Fix an owner and an anchor
family `A`.  For an anchor-unit tuple write

```text
Q = union_(i in A) M_i(u_i).
```

For every choice of the nonanchor units,

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(i notin A) |M_i(u_i) minus Q|
 <= |Q| + sum_(i notin A) max_e |M_i(e) minus Q|.              (8)
```

Consequently

```text
U_A = max_(anchor-unit tuples)
      [|Q| + sum_(i notin A) max_e |M_i(e) minus Q|]           (9)
```

is a sound upper bound on every literal union.  The independent maxima in
(8)-(9) deliberately forget shared nonanchor units and pair overlaps, so they
can only enlarge the attainable union.  An owner with `U_A<28` is terminal.

## 4. Why the `Z/4` premise is proved

For every effective order `D`, replacing `t` by `t+D` in (4) changes its CRT
argument by `13DB`, which is zero modulo `13D`.  Therefore an order-`D` mask
is the pullback of a subset of `Z/DZ`.

In particular, every order-one, order-two, or order-four mask is periodic
under `t -> t+4`, so it is a valid `Z/28 -> Z/4` anchor.  Order one is absent
from (7)'s survivors, leaving precisely the order-two and order-four anchor
providers.  This proves the periodicity premise of THM-994 rather than
inferring it from the observed masks.

## 5. The terminal two-adic certificate

For every one of the `3170` labelled scalar rows, and for each of its six
owners, the primary retains all order-two/order-four masks exactly and applies
(9) independently to every order-seven/order-fourteen/order-twenty-eight
provider.  The `19020` anchor-union banks have only one, two, or three masks:

```text
bank size:       1:12888, 2:936, 3:5196.                       (10)
```

The resulting owner-bound histogram is

```text
U_4: 18:912, 19:900, 20:1032, 22:2064,
     23:8328, 24:4392, 28:1392.                               (11)
```

Most importantly, the number of owners per labelled context at which this
relaxed upper bound still reaches `28` is

```text
0:2018, 1:912, 2:240.                                         (12)
```

There is no row with three, four, five, or six surviving owner projections.
A global unit word would have to cover all `28` sheets at every one of its six
owners.  By (8), every scalar row instead has at least four owner projections
that are impossible for **all** unit choices.  Combined with the scalar
elimination preceding (7), this proves the emptiness statement.  Section 9
records two independent replays of the implication.

No orbit quotient enters this implication.  Multiplication by `F_13^*` splits
the `3170` labelled contexts into orbit-size histogram

```text
4:2, 6:3, 12:262
```

and the `206` supports into `2:1, 12:17`; the primary nevertheless visits all
`3170` rows and all `19020` owner obligations.

## 6. The competing flags and the exact residual answer

The 7-adic carrier is much weaker.  Retaining order-seven masks gives

```text
U_7 live owners/context:
0:304, 1:408, 2:390, 3:144, 6:1924.                            (13)
```

Thus `1924` rows survive at all six owners under the `Z/7` relaxation.  This
failure is uniform at the owner level: on every one of the `19020` obligations

```text
U_4 <= U_7.                                                    (14)
```

Hence the per-owner competing-flag certificate satisfies
`min(U_4,U_7)=U_4` everywhere.

The primary also retains the union of both anchor families in the `4 x 7` CRT
product.  On the scalar survivor bank it finds the stronger exact identity

```text
U_{4x7} = U_4       on all 19020 owner obligations.            (15)
```

Equation (15) is a finite-bank identity, not a claimed all-parameter algebraic
law.  It says that resolving the transverse order-seven fibres does not lower
even one of these bounds.  There is therefore no residual carrier to close:
the first `Z/4` flag already closes every row.

As a sharpness sidecar, immutable-union DP retains every provider mask.  Across
all `19020` owner obligations it constructs `6628500` reachable masks, with
largest bank `2438`, and obtains

```text
exact maximum union:
18:912, 19:900, 20:1056, 22:2088, 23:8280, 24:4392, 28:1392;

exact feasible owners/context:
0:2018, 1:912, 2:240.                                         (16)
```

The mod-four upper bound exceeds the exact maximum by `0` on `18948` owner
rows, by `1` on `48`, and by `2` on `24`; every discrepancy lies strictly
below threshold.  Thus the deliberately lossy `Z/4` relaxation is
**threshold-exact on this entire bank**:

```text
U_4 >= 28    iff    the literal owner union can equal Z/28.    (17)
```

The DP is not needed for the proof implication in Section 5.  It measures the
information discarded by (8) and supplies a stringent internal cross-check.

## 7. The `4 x 7` toothpick/Kakeya picture

Identify `Z/28` with the CRT grid `Z/4 x Z/7`.

- An order-four mask is a seven-cell thick fibre at fixed mod-four coordinate;
  an order-two mask is the union of two such fibres with fixed parity.
- An order-seven mask is a transverse four-cell fibre (or a union of two).
- An order-fourteen mask is assembled from two-cell parity strips.
- An order-twenty-eight mask consists of four or five point needles.

The surprising structural fact is that the thick 2-adic stalk already bounds
every possible transverse completion.  We may independently slide each
7-adic fibre, two-cell strip, and point needle to its best location outside the
retained stalk and still cannot rescue more than two owners.  This is the
scale-28 analogue of THM-994's nested toothpick carrier, but bifiltration does
not mean symmetry: the 2-adic direction dominates the 7-adic one.

The loss ledger is explicit.  The `Z/4` quotient forgets transverse mod-seven
positions, pair overlaps among relaxed order-seven/order-fourteen/order-28
masks, and their shared unit constraints.  Those losses are sound only because
(8) is an upper relaxation.  The full DP in (16)-(17) shows that none of the
forgotten incidence changes the threshold answer.

## 8. Cayley and tournament audit

The scalar ratio switches already distinguish the two flags.  On the twelve
nonzero labels, direct the Cayley arc `a -> b` when `b/a` lies in the
off-diagonal high-cardinality ratio set.

```text
order 2: arcs 24, symmetric edges 24, reciprocal 0,
         directed triangles 0, one SCC;
order 4: arcs 72, symmetric edges 48, reciprocal 24,
         directed triangles 64, one SCC;
order 7: no off-diagonal high-cardinality arc.                 (18)
```

The symmetric order-two shadow is `K6,6` on quadratic versus nonquadratic
residues with three disjoint four-cycles removed.  The symmetric order-four
shadow is `K4,4,4` on the three cosets of the order-four subgroup

```text
{1,5,8,12}, {2,3,10,11}, {4,6,7,9}.                           (19)
```

These are Cayley relations, not tournaments: missing and reciprocal pairs
matter.  In particular the order-seven scalar switch has no off-diagonal
relation at all, matching the weakness of (13), while the two-adic flags carry
nontrivial label geometry.

For the mandatory tournament sidecar, the useful vertices are **owner proof
obligations**, not runners.  Compare owners lexicographically by

```text
(1[U_4>=28], U_4, exact maximum, scalar capacity, U_7),        (20)
```

let the harder key win, and orient equal keys along coordinate order as the tie
Hamiltonian path.  All `3170` completed tournaments are transitive: score
multiset `{0,1,2,3,4,5}`, no directed cycle, six singleton SCCs, and one
Hamiltonian path.  This tournament is telemetry only.  It ranks obligations
but forgets the absolute threshold, so it cannot prove Section 5.

The challenged-assumption audit considered runners/providers, gaps, fixed
circle sections, section boundaries, wall-crossing events, residues, cover
arcs, Fourier modes, matroid circuits, Fano points, chi-seven colours, and
proof obligations as possible vertices.  Owner obligations with their
absolute bounds preserve the terminal predicate.  Runners preserve the input
but not the ownerwise deficit; isolated residues and boundaries destroy fibre
incidence; cover arcs are exact but uncompressed; Fourier modes lose monotone
set union; matroid circuits remember dependencies but not absolute coverage;
Fano/chi-seven colours and completed tournaments discard the `28`-cell
threshold.  This is why (20) is the honest tournament carrier, and also why
the tournament remains a sidecar rather than the proof object.

## 9. Reproducibility and independent replays

The deterministic primary:

1. reconstructs all CRT bases algebraically and by literal search;
2. checks the hereditary valuation grammar and both inclusion-exclusion counts;
3. checks every literal mask against (5)-(6) and proves its periodicity;
4. scans all `24911964` labelled scalar contexts;
5. evaluates `U_4`, `U_7`, and `U_{4x7}` at every survivor/owner obligation;
6. separately constructs every exact reachable-union bank; and
7. records orbit, Cayley, tournament, carrier-loss, and alternate-vertex
   telemetry without quotienting the proof scan.

Normal and `python -O` runs reproduce the primary output byte-for-byte.  The
independent C++ referee compiles cleanly under Apple clang 17 in optimized,
unoptimized, ASan/UBSan, and libc++ debug-hardening modes; all four executable
outputs are byte-identical, and the clang static analyzer reports no
diagnostic.  Frozen SHA-256 values are

```text
Python primary source  0f200ee79d569d611b299e95397ebdd75e348ff2883327c108cdef7c4a412961
primary output         960ed394fdad4b69d17fd23a9d5617d9bd2de48588c5dc131c67045fcb634010
C++ referee source     d0bec69db7571781704e0dc906c258816e129a8563e66c1406ef894a7a6973b8
C++ referee output     ad00571d662430ac14906673b743fd0327d7e379ad0ab3bb48289bb3aa9eb7f1
```

Internal stream digests are

```text
grammar       0d5fb0d79a2c398857b8ba4d286729872b828a8092c17055ef9d06ba72aee366
CRT bases     c9d6db1a24f3960bf71a9c711cc46ce9abf72c26f27975ab16f9c03852cce7a5
literal masks dd9d2fe84c27d15486cd0d67f1237e9b3233455768ca80994c396e35dee8ac1a
scalar rows   311dc8eccfee377d69f134b4bd3d2f335d98dff6031c276789865aa92fb3f6ea
fibre rows    d4d9362a069e5fe445a00d221907d7e170906e1b7ebf34528590c03711fc2819
exact banks   5561960152be08a1b7d61cb8b23e40e6ea51fa25c18086c117823bd9db06d000
```

The frozen primary and its exact-DP sidecar share mask construction and are
therefore not independent of each other.  The promotion instead rests on two
separately structured referees.

The first referee is a literal-CRT Python flag certificate.  It normalizes the
3,170 rows to 1,585 owner keys, stores saturated capacities of the four
mod-four fibres, and separately reconstructs all exact banks.  It agrees on
the scalar histogram, exact feasible-owner census `0:2018,1:912,2:240`, exact
maximum histogram, all 6,628,500 reachable-mask incidences, and zero threshold
mismatches.  Its mod-seven flag has 6,072 owner false positives, confirming
that the successful quotient is genuinely the two-adic one.  Normal and `-O`
runs are byte-identical.

The second referee is standard-library C++.  It changes language, hash
scheme, mask ordering, row representation, and exact-DP traversal.  It checks
the literal-cover-implies-relaxed-cover implication on all 19,020 labelled
owner obligations and independently reproduces every decisive census.  Clang
optimized/unoptimized and ASan/UBSan builds produce the frozen output
byte-for-byte.  Referee hashes are

```text
Python flag source/output  3482a5a103891342f25448c0cc9ca0613494b04b1a05d57221f97a089f38d1f1
                           d815c9197ec3e6b71436485713e763f80621c9cc7a47c064cc728a485313ffa5
C++ referee source/output  d0bec69db7571781704e0dc906c258816e129a8563e66c1406ef894a7a6973b8
                           ad00571d662430ac14906673b743fd0327d7e379ad0ab3bb48289bb3aa9eb7f1
```

The exact-versus-flag contingency is

```text
(exact false, flag false): 17628,
(exact true,  flag true):   1392,
```

with no off-diagonal entry.  The two tournament serializations are
deliberately only sidecars: the primary ends its key with `U_7`, whereas the
C++ referee ends its key with exact-bank size.  Both are transitive on every
row, but agreement of those nonidentical keys is neither assumed nor used by
the proof.

This theorem concerns only the AP-centred Hamming-six common-scale-28 face.  It
does not address Hamming five, non-AP/deep sheets, the global sporadic branch,
or the full LRC14 theorem.  Scale `29` is prime-excluded by THM-983; the next
composite scale should be selected only after a fresh frontier audit.
