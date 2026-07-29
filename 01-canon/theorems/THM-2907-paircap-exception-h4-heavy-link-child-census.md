---
id: THM-2907
title: "Pair-cap-exception H4 heavy-link child census and endpoint repair"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  All 18,290 actual H4 pair flags on
  THM-2901's 52 pair-cap exceptions are decided by an exact heavy-link
  recursion.  Literal deletion closes 18,280 flags; the ten survivors
  reconstruct exactly two endpoint-tight 13-speed families, both lonely
  at t=1/28.  Thus all 52 exception branches are discharged, while every
  exception body still has a separate H3 branch and no whole root closes.
source: codex/h4-exception-children-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
related:
  - THM-592-radius-derivative-structure
  - THM-1115-no-blocking-tradeoff-no-gap
verification:
  - 04-computation/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.out
  - 05-knowledge/results/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.ledger.out
  - 04-computation/lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.out
---

# THM-2907 -- pair-cap-exception H4 heavy-link child census and endpoint repair

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Inherited exception universe

Fix one of the `52` exact pair-cap exceptions in THM-2901.  Retain its
literal carrier `C`, mass `h`, attained allowed singleton maximum `q_1`,
and complete forbidden prefix `P`.  THM-2899 gives

```text
q_1<3h/7.
```

Apply THM-2893 with `(k,s)=(5,4)` and singleton cap `B=q_1`.  Every
hypothetical allowed five-cover contains at least two labels in

```text
H_4(C)={w:c_C(w)>=(h-q_1)/4}.                              (1)
```

The strict displayed inequality makes `(1)` finite under THM-735.  The
locked predecessor census reconstructs all `52` carriers and their actual,
globally sealed cores: the core sizes range from `12` to `44`, their sizes
sum to `1,348`, and they have exactly

```text
sum_C binom(|H_4(C)|,2)=18,290                              (2)
```

unordered pair flags.  These actual memberships, rather than cutoff
binomial universes, are the input to this theorem.

## 2. The binary heavy link

Select an unordered pair `L subset H_4(C)` and put

```text
R=C minus D_L,          h_R=|R|,          theta_L=h_R-q_1. (3)
```

If `L` belongs to a hypothetical five-cover and `T` is its remaining
three labels, then every pair `A subset T` satisfies

```text
U_R(A)>=theta_L.                                           (4)
```

Indeed, `L union A` is a four-subset of the parent cover, hence is heavy
by THM-2893.  The literal link identity gives

```text
U_C(L union A)=U_C(L)+U_R(A)>=h-q_1,
```

which is exactly `(4)`.  Thus the child is not an arbitrary three-cover:
its three vertices span a triangle in the symmetric binary link `(4)`.
The selected pair and every label in `P` remain forbidden.

When

```text
theta_L/2>h_R/7,
```

define the globally sealed core

```text
H_2(R)={w:c_R(w)>=theta_L/2}.                              (5)
```

Every edge satisfying `(4)` meets `(5)`, since
`U_R({x,y})<=c_R(x)+c_R(y)`.  A triangle all of whose edges meet `H_2`
has at least two vertices in `H_2`: otherwise its two outside vertices
would form a missed edge.  It is therefore exhaustive to enumerate an
unordered `H_2` pair, check its union literally, and solve for the final
singleton while checking all three link edges.

## 3. Uniform treatment of the nonfinite-link lane

The `H_2` shortcut is not invoked when its discrepancy gap is nonpositive.
It is also deliberately bypassed when its analytic head exceeds `2,000`
or its actual inspected core exceeds `12`; those two choices affect cost,
not the search universe.  In all three cases the verifier uses this exact
generic recursion:

1. a three-cover of `R` has a label with coverage at least `h_R/3`;
2. after deleting it, the remaining two-cover has a label covering at
   least half of the new residual;
3. after deleting that label, the final label must cover the literal
   singleton residual.

The first two thresholds are strictly above one seventh and therefore
have finite THM-735 cores.  For the last step, if `lambda` is the longest
component of the nonempty residual, a danger comb that contains it must
satisfy

```text
w<=floor(1/(7lambda)).                                     (6)
```

Hence the last search is also finite.  Every equality candidate is checked
again by literal residual subtraction; every full triangle is checked
against all three inequalities `(4)`; and every reconstructed five-set is
checked directly against `C`.  This route handles all `179` nonpositive
`H_2` gaps without assuming a nonexistent finite core.

## 4. Exact census

With the fixed optimization thresholds above, the exact route split is

```text
finite H2 route                                      16,357
generic three-cover route                             1,933
  nonpositive H2 gap                                    179
  analytic H2 head above 2,000                          529
  actual inspected H2 core above 12                   1,225
                                                       -----
all H4 pair flags                                    18,290. (7)
```

Among the analytic heads that are inspected before route selection, the
maximum actual `H_2` size is `49`; the finite-`H_2` route itself uses only
sizes at most `12`.  Across both routes the verifier evaluates `184,429`
pair unions and `1,441,996`
terminal singleton-head labels; the maximum terminal singleton cutoff is
only `120`.  The generic route inspects `2,443` first-core and `628`
second-core vertices, with maximum cutoffs `1,241` and `746`.

Literal deletion leaves

```text
empty child completion                               18,280 flags
nonempty child completion                                10 flags
literal-closed exception bodies                          51
literal-open exception bodies                             1. (8)
```

All ten nonempty flags occur on

```text
E=(2,4,6,8,10,12,14),   rank=1,   apex=22,   P={22},
```

and their selected pairs are

```text
(16,18), (16,20), (16,26), (16,48), (18,20),
(18,26), (18,48), (20,26), (20,48), (26,48).             (9)
```

Exhaustive enumeration gives `16` triangle occurrences across `(9)`, but
after reconstructing and deduplicating the parent five-set there are
exactly two possibilities:

```text
(16,18,20,24,26),          (16,18,20,26,48).             (10)
```

Consequently the only full 13-speed completions are

```text
V_A=(2,4,6,8,10,12,14,16,18,20,22,24,26)
   =2*{1,2,...,13},

V_B=(2,4,6,8,10,12,14,16,18,20,22,26,48)
   =2*{1,2,...,11,13,24}.                                 (11)
```

## 5. Endpoint repair

The empty-residual engine works on interval carriers and therefore forgets
measure-zero equality endpoints.  That quotient loss is decisive in
`(11)`.  The original danger predicate is the strict open condition
`||vt||<1/14` (THM-592), whereas the pinned interval engine deliberately
subtracts closed radius-`1/(14v)` comb pieces.  An independent exact
hostile audit finds, for each family in `(11)`,

```text
strict-danger owners at t=1/28       none
equality-boundary owners              2,26
closed-comb good carrier              empty, mass 0.
```

Thus the closed-comb engine really does erase the candidate point, while
the original non-strict target retains it.  At

```text
t=1/28
```

the circle distances for either doubled family equal the distances of its
primitive labels divided by `14`.  For `V_A/2`, the residues are
`1,...,13` modulo `14`.  For `V_B/2`, they are
`1,...,11,13,24=10 (mod 14)`.  None is zero, so in both cases

```text
min_(v in V_A)||vt||=min_(v in V_B)||vt||=1/14.           (12)
```

The hostile perturbations `t=1/28 plus or minus 10^(-9)` both have
clearance strictly below `1/14`, so this is genuinely an isolated equality
repair, not an unrecorded positive-measure interval.

Thus both literal survivors satisfy the non-strict lonely-runner target.
The second primitive family is the tight non-AP family already recorded
in THM-1115, but `(12)` is a direct modular proof and does not depend on
its exact tightness classification.

It follows from `(8)`--`(12)` that every one of the `52` pair-cap
exception branches is discharged: `51` have no literal completion, and
the last has only the two endpoint-safe completions in `(11)`.

## 6. Composition and scope

Relative to the pinned THM-2901 partition, the branch ledger may now be
read as

```text
direct pair-partition closures                         1,835
remaining finite H3 pair-flag obligations             12,919
H4 pair-cap exceptions discharged                         52. (13)
```

Every one of the `52` exception bodies still has at least one distinct
`H_3` branch in that same pinned ledger.  Therefore this theorem alone
closes no whole seven-body root.  It does not decide any of the `12,919`
`H_3` obligations, close the current all-root residual, or prove LRC(14).

No later root-composition theorem is used in either the proof or the
recomposition statement above.

## 7. Verification

The verifier hash-pins the predecessor H4 membership census, which in turn
pins THM-2901's engine and complete branch ledger.  It uses exact rational
interval arithmetic and explicit guards rather than Python `assert`.
Normal and optimized replays are byte-identical.  The semantic digest is

```text
e8d4119f101ae3ac310fe5ca8a056607390ff4d7aa166cb90168983ea7069356.
```

Canonical artifacts:

```text
04-computation/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.py
SHA-256 88d523ea97235471ecce03c06de5cd1e1ba434ccd41fe0633beadf1017aa8fa3

05-knowledge/results/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.out
SHA-256 1df929e106cd16c094886d59f3702ba9bafa395ee906527fe4592a1552e9b458

05-knowledge/results/lrc14_j6_paircap_exception_h4_link_child_census_codex_20260729.ledger.out
SHA-256 2aad1013bd739d82b407d5632e2cd2e9763c1f9edb2bfa127500420bda6fdc9d

04-computation/lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.py
SHA-256 0c332f6522d0ad77185ba7bcbe25ec27cd68b755f652cefa717f3bc91fb74db2

05-knowledge/results/lrc14_j6_paircap_exception_endpoint_convention_audit_codex_20260729.out
SHA-256 a93a4a724dac6c55806f3358c2f5ab25de8f0261c92906a0161414781a717d20
```
