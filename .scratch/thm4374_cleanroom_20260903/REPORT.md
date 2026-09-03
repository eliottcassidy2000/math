# THM-4374 clean-room referee report

**Verdict:** PASS WITH WORDING REPAIRS.  The proposed theorem is correct as an
elementary result relative to THM-4365 and THM-4367, and the exact certificates
pass.  Its scope is one fixed `h=420` odd tail and one metric exit consumer.
LRC(14) remains **OPEN**.

## What was independently checked

The independent route was frozen in `DESIGN.md` before reading the proposed
canonical primary.  Its verifier uses integer rational pairs rather than
`fractions.Fraction`, explicitly checks the return formulae and boundary
phases, and implements a direct inverse decoder for `W_17`.

The clean-room run passed `3,390,840` explicit checks:

```text
half-phases                              23,597
structural (a,kappa mod 14) types         6,106
structural scale points                  13,427
unordered structural scale pairs       281,073
literal global decoder rows              47,194
first split counts       t=1:239,955, t=16:30,224, t=17:10,894
maximum W_H fibres       241, 94 (H=1..15), 49 (H=16), 1 (H=17)
```

The phase cover is exact:

```text
t=0: 3370; t=1,...,15: 1303 each; t=16: 682.
```

Both inactive boundary phases were included.  Their least cofinite
representatives are `(P,rho)=(43823,-3371)` and `(50565,3371)`, and both emit
the baseline.

The single sharp class

```text
(a,b,kappa)=(2,47595,18397),  g=1,15,...,3361
```

has 241 current scales, 94 low scales `g<=1303`, 49 middle scales
`631<=g<=1303`, and singleton `W_17` fibres.  The global minimality hostile

```text
P=253031,258645
```

has equal `W_16`, baseline at times 1 through 16, and distinct time-17 exits

```text
97819/3542910, 99989/3621506.
```

## Universal quantifier audit

The finite structural census is not being used as an unbounded empirical
claim.  Every actual strict-active fibre has

```text
c=ag, P=bg, a+1303b=3371kappa, gkappa=1 (mod 14),
```

with even `a` and `ag<=6740`.  For fixed `(a,b,kappa)`, the tail condition
only removes an initial segment of the single arithmetic progression of
scales.  The word partition depends on `a*g` and the residue of `g` modulo
14, while its only possible common-return collision reduces to
`t*kappa=7b`.  This is impossible at `t=16,17` modulo 7.  Therefore the
finite census of the full structural progressions proves a stronger finite
classification whose restriction contains every actual tail fibre.

The return formulae were checked literally on all 13,427 structural scale
points:

```text
c>=2608: E(bg+2)=gkappa/[14(bg+2)];
c<=1242: E(bg+32)=(gkappa+14)/[14(bg+32)];
c<=2606: E(bg+34)=(gkappa+14)/[14(bg+34)].
```

The exact boundary cuts `c=2606/2608` and `c=1242/1244` were included.

For global injectivity, the rotation intervals cover every phase by time 16.
A delayed first entry is followed by another active output.  Two consecutive
active excesses `d_0,d_1` recover the shifted parameter `Q` by

```text
Q=2(1303+47194*d_1)/(47194(d_0-d_1)).
```

For a current low active phase, the time-16 or time-17 once-wrapped return
recovers `P` by

```text
P=(47194-2t*1303-2t*47194*d_t)/(47194(d_t-d_0)).
```

The denominator can vanish only under the already impossible return equation.
This is a direct inverse for `W_17`, not a finite-window extrapolation.  The
47,194-row, two-period replay is a hostile control for conventions only.

Finally, if an equivalence relation is output-compatible and preserved by
`P -> P+2`, iteration gives equal `W_17`; injectivity forces equal parameters.
The congruence corollary is therefore exact.

## Required wording repairs before promotion

1. Say explicitly that `W_H` has `H+1` outputs at times `0,...,H`.
   “Horizon 17” means 17 shifts and 18 observations.
2. State the `241/94/49/1` maxima on the **strict-active current-output
   locus**, equivalently as the largest fibres inside the sets `G` as the
   admissible triple varies.  On the whole tail the inactive `H=0` baseline
   fibre is infinite.
3. In the `r=0` global proof, insert the dependency step: equality of the
   current strict-active exits first puts the two parameters in the same
   `(kappa,b)` fibre by THM-4367, after which Theorem A applies.
4. Replace “the coarsest output-compatible shift congruence is equality” by
   “equality is the only output-compatible forward-shift congruence.”  This
   avoids an order-convention ambiguity.
5. Define an “exact state” in the finite-sidecar sentence as an
   output-emitting deterministic quotient whose kernel is a forward-shift
   congruence.  Then the infinite inactive baseline fibre plus a finite tag
   gives the stated pigeonhole obstruction.
6. Describe the 281,073-pair job as an exhaustion of **structural scale-pair
   types**, after the universal reduction above, rather than as an
   enumeration of the infinitely many actual `(b,kappa)` fibres.
7. On promotion, use status `PROVED ELEMENTARY RELATIVE TO THM-4365 AND
   THM-4367 + VERIFIED-EXACT + INDEPENDENTLY AUDITED`, and add both theorem
   slugs to `depends_on`.  Do not label the global theorem merely
   `FINITE-EXACT`.
8. Treat the existing scratch verifier as a same-design corroboration, not the
   sole independent referee: it and the proposed primary share the same
   structural enumeration skeleton and the same 922,895-check layout.  The
   present direct-decoder verifier supplies the independent audit layer.

None of these repairs changes a constant, witness, formula, implication, or
proof conclusion.

## Execution and hashes

Run from the repository root:

```text
python3 -B .scratch/thm4374_cleanroom_20260903/verify_thm4374_cleanroom.py
python3 -B .scratch/thm4374_cleanroom_20260903/audit_execution_modes.py
```

The mode audit ran normal, `-O`, `-I`, and `PYTHONHASHSEED=314159265`.  Every
stream used raw LF bytes, all four streams agreed within each target, and the
proposed primary also matched its frozen output.

```text
clean-room source SHA-256
  8cc318c4922d00a6f0d8e8879c41bf7e0fd68521b0bd57858627c1fb2ba7124d
clean-room output SHA-256
  03fb2f11e2212fa44d5ff489e793efb8fb71cdafa252e2b266e0d70e4a6aa1c5
mode-audit source SHA-256
  6489520411da145192a735ba90fe2dbf5e8b12f0930c181ff80d00b54edb351c
mode-audit output SHA-256
  baaf6a16d350009e622949c3c17edf09d77611cc26d24398a2d3f30759902f6f
proposed primary source SHA-256
  d4c73423389af2341cfb382f7c189e7fb0339280187bcaade5e4da4b5ea376ba
proposed primary output SHA-256
  9239480fe4dadeaf635f52188517b481f7285f846424e73f943e7af5b422e40c
```

The packet's pre-existing verifier and frozen output also replayed in all four
modes with their advertised hashes
`d39c8ad3e133ad084835602f2dfc714991d6faa4e89790e009c1436f86d86a39`
and `78a00f403babf36a9bf5e485e4736e61df0987552c95e25cb7f373eb1e7d14b5`.

## Strict non-consequences

This theorem identifies an odd-tail parameter from a delayed word of one
first-exit metric.  It does not observe the full physical binder unless `P`
is reconstructed, does not improve the equality-safe time, does not show
entry into the unresolved seam, does not classify arbitrary speed banks or
speed changes, does not decrement an LRC ledger, and does not prove LRC(14).
