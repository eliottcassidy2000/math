---
id: THM-4125
title: "Arbitrary multi-tail projective packing and sharp residue density"
status: >
  PROVED ELEMENTARY PROJECTIVE PACKING + SHARP NATURAL-THRESHOLD TRANSITION +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. For the fixed eleven-speed body U
  and any finite set C of r distinct positive tail multipliers, the phase
  9/19 has clearance distribution 0, 1/19, 2/19 determined exactly by the
  signed projective support of C modulo 19. If no multiplier is zero modulo
  19, the fixed 1/14-safe density is (18-2m)/19, where m is the signed-support
  size, with sharp maximum 16/19 for every r. At the natural threshold
  1/(12+r), the same count holds for r<=6, but every nonzero parameter class
  is safe for r>=7; equality occurs at r=7 and the density is 18/19. A zero
  multiplier residue collapses the displayed phase. LRC(14) remains open.
source: codex-frontier-synthesis-creative-20260825g
depends_on: []
related:
  - THM-764-covering-small-period-signed-pair-deck-and-q25-refutation
  - THM-4119-infinite-supplier-free-eleven-plus-two-residue-family
  - THM-4121-sharp-projective-tail-multiplier-residue-compiler
  - HYP-7812
script: 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125.py
output: 05-knowledge/results/lrc_arbitrary_multitail_projective_packing_thm4125.out
independent_audit_script: 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_arbitrary_multitail_projective_packing_thm4125_independent_audit.out
script_sha256: 82fef5178d83645ba72aafb3c4aadd0dd72e299d9c42bfbfb2d92eab6c044915
output_sha256: 7b659ae86dd858b592d45ee76022cc5adb7f1e2d88558eda342723371b761822
semantic_sha256: 2b80075c494e496faa4247cbe22d59cb07af20b31cc7c426e82d7be8e5dc7192
independent_audit_script_sha256: 73f85fc4997d4085d8f5d29f5b21dd576925e8cdfc1929d0b02f298e6a236aff
independent_audit_output_sha256: a9a622b7d12006c1e4dd93e49abad5d8d2d5088f9904fc379d2493b0a4f2ab61
independent_semantic_sha256: 2b80075c494e496faa4247cbe22d59cb07af20b31cc7c426e82d7be8e5dc7192
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone Fraction implementation exhausts all 262,143 nonempty
  subsets of the eighteen nonzero residue classes, checks the complete
  0/1/2 clearance law, signed-support histogram, ordinary- and
  signed-distinct optima, and literal aligned families at selected tail counts
  up to 40 across the r=6/7 threshold boundary. Normal, optimized, and
  frozen outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room integer implementation imports no primary code and
  enumerates every residue subset as nine independent four-state signed
  pairs. It reproduces the same clearance counts and semantic ledger, checks
  the natural threshold by integer cross-multiplication, and replays the zero
  residue hostile and aligned families. Normal, optimized, and frozen outputs
  byte-match; no theorem failure occurs.
---

# THM-4125 -- arbitrary multi-tail projective packing and sharp density

**PROVED ELEMENTARY PROJECTIVE PACKING + SHARP NATURAL-THRESHOLD
TRANSITION + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4121 aligns the two bad residue pairs of two tails. The same mechanism
does not saturate when more actual tails are added: the obstruction counts
their **signed projective support**, not their number. This separates two
quantities that coincide in the thirteen-speed case. The fixed `1/14`
benchmark retains sharp density `16/19` for arbitrarily many aligned tails,
whereas the natural lonely-runner threshold drops enough at seven tails to
admit every nonzero parameter class.

## 1. Family and projective support

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22).
```

Let `C={c_1,...,c_r}` be any set of `r>=1` distinct positive integers and,
for an integer `t>22`, define

```text
S_(C,t)=U union {ct:c in C}.                              (1)
```

The set in `(1)` has exactly `11+r` distinct positive speeds: all tails are
larger than `22`, and distinct multipliers give distinct tails. Write

```text
|a|_19=min(a mod 19,19-(a mod 19)),
delta_C(t)=min_(v in S_(C,t)) ||9v/19||.                  (2)
```

The body has exact clearance

```text
min_(u in U)||9u/19||=2/19,                              (3)
```

attained, for example, by `u=4` and `u=15`.

Suppose first that every multiplier is nonzero modulo `19`. In the
nine-element quotient

```text
G=F_19^*/{+1,-1},
P(C)={[c mod 19]:c in C},                 m=|P(C)|,        (4)
```

let `m` be the signed projective support size. Since `9^(-1)=17 mod 19`,
define the parameter blocker set

```text
B_C={0} union union_(c in C){+/-17(c mod 19)^(-1)}.       (5)
```

> **Theorem 1 (complete multi-tail clearance law).** If some `c in C` is
> zero modulo `19`, then `delta_C(t)=0` for every admissible integer `t>22`.
> Otherwise
> `(5)` has cardinality `1+2m`, and
>
> ```text
> delta_C(t)=0       if t=0 mod 19,
> delta_C(t)=1/19    if t mod 19 in B_C minus {0},
> delta_C(t)=2/19    if t mod 19 is outside B_C.          (6)
> ```

### Proof

At phase `9/19`, a tail `ct` has numerator

```text
|9ct|_19.                                                 (7)
```

If `c=0 mod 19`, `(7)` vanishes for every parameter. Suppose all multiplier
residues are nonzero. At `t=0 mod 19`, every tail vanishes. At a nonzero
parameter, `(7)` is one exactly when

```text
t=+/-17c^(-1) mod 19.                                    (8)
```

It is otherwise at least two. The two-element sets in `(8)` agree exactly
when the corresponding multipliers agree up to sign. Inversion and
multiplication by `17` induce a bijection of `G`, so distinct elements of
`P(C)` give disjoint signed blocker pairs. Their union with zero therefore
has size `1+2m`. Combining the tail table with the body equality `(3)` gives
all three lines of `(6)`. **QED.**

## 2. Fixed `1/14` density and all sharp packing regimes

Because

```text
1/19 < 1/14 < 2/19,                                      (9)
```

the fixed `1/14` benchmark admits exactly the last line of `(6)`.

> **Theorem 2 (sharp projective packing).** If every multiplier is nonzero
> modulo `19`, the natural density of parameters `t>22` that are safe at the
> fixed `1/14` benchmark and phase `9/19` is
>
> ```text
> rho_14(C)=(18-2m)/19.                                  (10)
> ```
>
> For every fixed number `r` of distinct actual multipliers, the maximum of
> `(10)` is `16/19`. Equality holds exactly when all multiplier residues lie
> in one signed projective class. Arbitrarily many tails attain it, for
> example
>
> ```text
> C_r={1,20,39,...,1+19(r-1)}.                           (11)
> ```

There are `19-(1+2m)=18-2m` last-line residues in `(6)`, proving `(10)`.
Every nonempty support has `m>=1`, so its numerator is at most sixteen; the
family `(11)` has distinct actual multipliers but `m=1`. Conversely, equality
forces `m=1`. **QED.**

The whole packing spectrum is equally explicit.

- With no residue restriction, every `m` in
  `1<=m<=min(r,9)` occurs, so the possible nonzero-residue densities are
  precisely `(18-2m)/19`. Including a multiplier divisible by `19` gives
  density zero.
- If the ordinary residues are required to be pairwise distinct and
  nonzero, then `r<=18`, each signed class contains at most two of them, and
  `m>=ceil(r/2)`. The sharp maximum number of safe classes is
  `18-2ceil(r/2)`, attained by filling signed pairs. If nineteen ordinary
  residues are required, zero is forced and the displayed phase collapses.
- If the signed projective classes themselves are distinct, then `r<=9`,
  `m=r`, and the safe-class count is exactly `18-2r`.

Thus the support-size table is

| signed support `m` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fixed-safe classes | 16 | 14 | 12 | 10 | 8 | 6 | 4 | 2 | 0 |

The examples `C={1,18}` and `C={1,2}` distinguish ordinary-residue
distinctness from projective distinctness: their support sizes are one and
two, hence their densities are respectively `16/19` and `14/19`.

## 3. The seven-tail natural-threshold transition

The `11+r` moving speeds in `(1)` have natural lonely-runner threshold

```text
lambda_r=1/(12+r).                                       (12)
```

The two positive levels in `(6)` compare with `(12)` by

```text
2/19-lambda_r=(2r+5)/(19(12+r))>0,
1/19-lambda_r=(r-7)/(19(12+r)).                          (13)
```

> **Theorem 3 (sharp natural-threshold phase change).** Suppose no multiplier
> is divisible by `19`. At the phase `9/19`, the number of parameter classes
> meeting the closed natural threshold `(12)` is
>
> ```text
> 18-2m,       1<=r<=6,
> 18,          r>=7.                                    (14)
> ```
>
> At `r=7`, every nonzero bad class meets the threshold with equality
> `1/19=lambda_7`. For `r>=8`, those classes meet it strictly. In all cases
> the zero parameter class fails.

For `r<=6`, the second quantity in `(13)` is negative, so only the `2/19`
classes survive and Theorem 2 gives the first line of `(14)`. At `r=7`, both
positive levels in `(6)` meet the weak lonely-runner inequality. For `r>=8`
both exceed it. Since every nonzero parameter has one of those two levels,
all eighteen nonzero classes survive. **QED.**

The threshold transition is independent of projective diversity. For the
maximally spread control `C={1,2,...,9}`, one has `m=9`: the fixed `1/14`
density is zero, but `r=9` makes the natural-threshold density `18/19`. For
the aligned family `(11)`, the fixed density stays `16/19` for every `r`,
while its natural density jumps from `16/19` to `18/19` at `r=7`.

## 4. Inheritance, hostiles, and information loss

The closest proved mechanism is THM-4121: its family is exactly the
specialization `C={1,c}`, and THM-4119 is `C={1,3}`. The canonical hostile is
`C={19}`, for which a zero-residue tail kills the phase for every `t`. The
least-used relevant sidecar is THM-764's signed unit-pair deck; the present
proof uses the same blocker alphabet but does not require THM-764's covering
hypothesis. The proved mod-19 spread lemma recorded in HYP-7812, after the
`forall`-to-`exists` repair in MISTAKE-186, is dual rather than duplicative:
it gives a global necessary support condition from failure below `2/19`,
whereas this theorem classifies one displayed phase as the parameter moves.

The projective compiler has the following exact information budget.

| item | content |
|---|---|
| source | nonzero multiplier residues `c mod 19` |
| target | signed parameter-blocker pairs in `G` |
| map | `[c] -> [17c^(-1)]` |
| preserved | signed support size and all blocker collisions |
| forgotten | integer lifts, ordering, and residue multiplicity |
| required sidecar | actual tail count `r`, because it sets `(12)` |
| cheapest decisive test | the nineteen-entry clearance table `(6)` |

Equal multiplier residues therefore co-locate at the displayed phase, and
opposite residues are reflected. This does not merge actual speeds: the
multipliers are distinct and `t>22`. Lonely runner permits different runners
to occupy the same phase position, so co-location is not an obstruction.

The fixed-density optimum is specific to `U`, phase `9/19`, and the benchmark
`1/14`. A bad residue in `(6)` only fails that displayed benchmark; it is not
a global nonloneliness claim. The natural-density statement is likewise an
explicit-family certificate, not an arbitrary-core closure. Only the
`C={1,c}` thirteen-speed specialization inherits THM-4121's uniform exclusion
from THM-4112's three explicit suppliers. No supplier statement is asserted
for general `C`, and LRC(14) remains open.

## 5. Exact replay

The primary audit uses exact Fractions, enumerates all `2^18-1=262,143`
nonempty subsets of `F_19^*`, and proves the signed-support histogram

```text
((1,27),(2,324),(3,2268),(4,10206),(5,30618),
 (6,61236),(7,78732),(8,59049),(9,19683)).                (15)
```

The independent audit instead chooses one of four states—absent, positive,
negative, or both—in each of the nine signed pairs. Both implementations
check the complete `0/1/2` table, all residue-distinct optima, the threshold
change by exact integer arithmetic, and literal aligned families with as many
as forty distinct actual tails. Run

```text
python3 -B 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125.py
python3 -B -O 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125.py
python3 -B 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125_independent_audit.py
python3 -B -O 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125_independent_audit.py
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_arbitrary_multitail_projective_packing_thm4125_independent_audit.py
```

All six streams match their frozen outputs byte-for-byte, and the independent
semantic ledgers agree. **QED.**
