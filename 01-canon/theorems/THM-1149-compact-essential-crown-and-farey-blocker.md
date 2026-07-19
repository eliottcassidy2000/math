---
id: THM-1149
title: Compact strict covers form an essential crown; a tight AP deletion is blocked by Farey regeneration
status: PROVED structural dichotomy, mass identity, deletion-margin bound, and AP-core blocker; FINITE-EXACT counterexamples and audits. Crown collapse, n=12 rigidity, INVcov, and LRC(14) remain OPEN
source: codex-2026-07-18-S74
depends_on:
  - THM-1142  # essential-region criterion
related:
  - THM-1099  # corrected compact sufficient residual
  - THM-1143  # shallow A12/Farey carrier
  - HYP-4382  # open twelve-speed equality classification
  - HYP-7665  # compact residual, corrected by this theorem
  - HYP-7675  # unsupported equivalence wording, corrected by this theorem
script: 04-computation/lrc14_compact_essential_crown_bridge_codex_20260718.py
output: 05-knowledge/results/lrc14_compact_essential_crown_bridge_codex_20260718.out
lean: 04-computation/lean/TournamentH7/TournamentH7/LRCCompactEssentialCrown.lean
---

# THM-1149 -- the compact bridge lands in an essential crown

Put `lambda=1/13` and, for a positive integer speed `v`, let

```text
B_v={t in R/Z : ||vt||<lambda}.
```

For a thirteen-speed set `V` and `v in V`, write

```text
W_v=V\{v},
E_v={t : ||wt||>=lambda for every w in W_v},
delta_v=M(W_v)-lambda,
m_v=max(W_v).
```

The strict/open convention on `B_v` is important.  It makes `E_v` closed,
and a hypothetical inequality `M(V)<lambda` makes the deleted runner strictly
dangerous everywhere on its essential region.

## 1. Essential-crown theorem

> **Theorem A.**  Suppose `V` has thirteen distinct positive speeds and
> `M(V)<1/13`.  Then, for every `v in V`:
>
> 1. `E_v` is nonempty;
> 2. `E_v subset B_v`, and the thirteen sets `E_v` are pairwise disjoint;
> 3. `delta_v>=0`; if `delta_v>0`, some component of `E_v` has length at
>    least `2 delta_v/m_v`, while every component has length strictly less
>    than `2/(13v)`.  In particular
>
>    ```text
>    0<delta_v<m_v/(13v).                                (1)
>    ```
>
> Consequently there is an exact dichotomy: either at least one deletion is
> tight, `M(W_v)=1/13`, or all thirteen deletions are loose and have
> positive-width, pairwise-disjoint private essential regions.  The second
> configuration is the **all-loose essential crown**.

**Proof.**  The known `n=12` case of the lonely runner conjecture gives
`M(W_v)>=1/13`, so `E_v` is nonempty and `delta_v>=0`.  At any point of
`E_v`, the strict inequality for the full family forces `v` to be dangerous;
hence `E_v subset B_v`.  A point in both `E_v` and `E_w` would be safe for
every speed of `V`, contradicting `M(V)<1/13`.

Let

```text
f_v(t)=min_(w in W_v) ||wt||.
```

This function is `m_v`-Lipschitz.  Around a maximizer, the closed interval of
radius `delta_v/m_v` therefore lies in `E_v`.  On the other hand, each
connected subset of `E_v subset B_v` lies in one open tooth of `B_v`, whose
length is `2/(13v)`.  A closed component cannot fill that open tooth.  This
gives the two component bounds and (1).  The dichotomy follows.  ∎

This is simultaneous information: a strict cover must make every deletion
near-tight relative to the width of the particular deleted runner's tooth.
It does not say that any deletion is exactly tight.

## 2. Private/triple multiplicity balance

Let

```text
N(t)=#{v in V : t in B_v},
mu_k=measure{t:N(t)=k}.
```

> **Theorem B.**  Under the hypotheses of Theorem A,
>
> ```text
> mu_1=sum_(k>=3) (k-2)mu_k.                              (2)
> ```
>
> Up to the measure-zero wall set, `{N=1}` is the disjoint union of the
> private essential regions `E_v`.

**Proof.**  Every danger comb has measure `2/13`, so

```text
integral N=13(2/13)=2.
```

Strict covering gives `N>=1`; hence `sum mu_k=1` and `sum k mu_k=2`.
Subtracting twice the first equality from the second proves (2).  A point has
multiplicity one with owner `v` exactly when every runner other than `v` is
safe, which is precisely membership in `E_v`.  ∎

Equation (2) identifies a real loss in pairwise telemetry: private mass is
paid for by multiplicity at least three, and pair-overlap totals do not retain
where those overlaps coincide.  The proof-bearing carrier is a labelled
private/triple-overlap chamber complex, not the naked tournament of runners.

## 3. The all-loose branch is real without Cover14

The exact family

```text
V0={1,2,3,5,7,8,9,10,11,12,17,19,104}                  (3)
```

has

```text
gcd(V0)=1,
rho(V0)=104/19<13,
M(V0)=8/105<1/13,
M(V0\{v})>1/13 for every v,
min_v M(V0\{v})=2/25.
```

It covers every modulus `2,...,13` but not `14`.  Its exact private mass is

```text
316140497/756575820,
```

whereas the sum of the thirteen Lipschitz component floors from Theorem A is

```text
352684666109/30512702820600.
```

The ratio is `36.1511...`.  Thus compactness, the component inequalities,
and (2) alone cannot extract a tight deletion.  The Cover14/lift information
must do genuine work.

Replacing `7` in (3) by `112=7+105` preserves the complete residue germ at
the old witness `8/105` and supplies Cover14, but the global optimum jumps to
`3/20`.  All six compact lifts `7 -> 112+210k`, `0<=k<=5`, have the same
value `3/20`.  This is the mandatory positive control for any proposed
crown-collapse invariant: it must see the integer lift even when one local
Farey chart does not.

The exact evaluations use the complete pair-sum ruler supplied by
THM-668/THM-1002, not a bounded-denominator heuristic.

## 4. Farey-toothpick regeneration blocker

For `d>=1`, the core `d[12]=d{1,...,12}` is tight at `1/13`.  Its explicit
maximizer set contains

```text
T_d={n/(13d): 0<=n<13d and 13 does not divide n},         (4)
```

which consists of `12d` affine copies of the base Farey teeth.

> **Theorem C.**  If
>
> ```text
> M(d[12] union {v})<1/13,                               (5)
> ```
>
> then `13d` divides `v`.

**Proof.**  Put `N=13d`, `g=gcd(v,N)`, and `Q=N/g`.  The values `vn/N` as
`n` runs modulo `N` visit each point of the `Q`-grid exactly `g` times.  The
number of points on that grid at strict distance below `1/13` from an integer
is

```text
b(Q)=1+2 floor((Q-1)/13).
```

If `Q>=2`, then

```text
13b(Q)<=2Q+11<12Q,
```

so fewer than `gb(Q)<12gQ/13=12d` numerator classes can be dangerous for
`v`.  They cannot contain all `12d` points of (4), yet (5) would require `v`
to kill every maximizer of the core.  Therefore `Q=1`, or `13d|v`.  ∎

> **Corollary D.**  No set `d[12] union {v}` satisfying (5) can be
> simultaneously primitive, carry a multiple of `14`, and have
> `rho<91/6`.  In particular it cannot satisfy primitive Cover14 and
> `rho<13`.

**Proof.**  Theorem C and primitivity force `d=1`.  The core `[12]` has no
multiple of `14`, so the 14-carrier must be `v`.  Hence `182|v`; because the
second-largest speed is `12`, `rho>=182/12=91/6>13`.  ∎

This is the exact use of the ladder's toothpick self-similarity: the extra
runner must kill all regenerated copies, not merely twelve representatives
near the origin.

## 5. Correct relation to the n=12 equality problem

Assume, hypothetically, the still-open classification

```text
|W|=12 and M(W)=1/13  =>  W=d[12].                       (R12)
```

Then Corollary D excludes the **tight-deletion** branch of a primitive
Cover14 compact strict cover.  It does not exclude the all-loose crown in
Theorem A.  The additional missing implication is

```text
primitive + Cover14 + rho<13 + M(V)<1/13
  => some v has M(V\{v})=1/13.                           (CROWN-COLLAPSE)
```

`CROWN-COLLAPSE` is open.  Therefore the following three statements must not
be identified:

1. the twelve-speed equality classification `(R12)`;
2. the stronger compact `1/13` sufficient residual `INVcov`;
3. LRC(14), whose target is only `1/14`.

THM-1099 gives a one-way reduction through `INVcov`; it does not prove the
reverse implications.  In particular, the S113/S114 claim that non-dilated
core rigidity, Tao's `n=12` probe, and LRC(14) are equivalent omitted both
tight-deletion extraction and a reverse embedding.

## 6. Finite exact stress tests

The companion computation additionally proves the following finite facts:

- all six compact lifts `7 -> 112+210k`, `0<=k<=5`, have exact value `3/20`;
- among `2,748` legal primitive Cover14 compact single substitutions around
  `{2,4,...,24} union {13}` with the added speed at most `500`, no row lies
  below `1/13`;
- the eleven equality survivors are obtained by replacing `13` with
  `39,65,...,299`, and all retain the regenerated core `2[12]`.

These are hard tests, not a uniform proof of crown collapse.

Normal and optimized runs are byte-identical.  Frozen hashes are

```text
source  335c42ae4d08426df09935bec629e62724e4b9fe3f15e1f0cd13bc089d2fcc55
output  db54526f86997674183824f9c1eec5688662be7ae147a2e31ae687a11dacb0ec
```

`LRCCompactEssentialCrown.lean` kernel-checks the finite algebraic form of
the private/triple mass balance and the complete post-extraction arithmetic
consumer: `13d|v` plus primitivity forces `d=1`, simultaneous 13/14
divisibility gives `182|v`, and positivity contradicts the cleared compact
ratio inequality.  It contains no `sorry`, `admit`, or `native_decide`.
The analytic crown construction, tight-deletion extraction, and Farey
regeneration theorem remain external exactly as stated above.

## 7. Carrier and tournament audit

For (3), ordering the 94 positive private stalks by cyclic midpoint produces
a transitive tournament: score histogram `0,...,93`, no directed triangles,
94 singleton SCCs, and one Hamiltonian path.  It preserves order and loses
tooth labels, simultaneous triple incidence, and Cover14 lift compatibility.

The smallest faithful object found here is instead

```text
edge/curtain event word
  + private-stalk and >=3-overlap chamber complex
  + owner/root tooth addresses
  + Cover14 lift-congruence sidecar.
```

The challenged assumption is explicit: tournament vertices need not be
runners, but even private-stalk vertices become insufficient after pairwise
orientation.  The exact residual is third-order and lift-sensitive.

## 8. Honest frontier

The theorem closes the AP-core branch *after* tight-deletion extraction and
proves that the alternative all-loose object is a rigid essential crown.  It
does not prove crown collapse, `(R12)`, the shallow all-height ballot lemma,
the deep `n=12` branch, `INVcov`, or LRC(14).
