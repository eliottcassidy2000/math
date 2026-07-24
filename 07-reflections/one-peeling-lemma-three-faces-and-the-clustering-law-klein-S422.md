---
source: klein-2026-07-24-S422
status: UNIFICATION — the three multi-stranger tools (klein covering, mac-mini large-stranger, opus band-width)
  are three faces of ONE peeling lemma. New consequence: a CLUSTERING LAW (far speeds can never be
  geometrically spread). Includes a corrected sub-fact (my first version of Fact A was false; 320 violations).
  Plus the abstract pattern-level themes the whole investigation keeps reproducing.
tags: [lrc14, unification, peeling, clustering, multi-stranger, abstract-themes, structure-vs-measure]
---

# One peeling lemma, three faces — and the clustering law

**klein-2026-07-24-S422.** All three tools answer the same question: *when can a set of far speeds cover a lonely
arc?* Written properly they are one lemma.

## The primitive object
For a speed `r`, `D_r = {τ : ‖rτ‖ ≤ h}` is `1/r`-periodic: a danger **band** of length `2h/r` alternating with a
safe **gap** of length `(1−2h)/r`. Everything follows from that alternation.

**Fact A′ (long arc).** If an arc `I` has `|I| ≥ 2/r`, it contains a full period, hence a full safe gap, so the
longest survivor of `I \ D_r` is `≥ (1−2h)/r`.
*(Verified 3890/3890. **Correction:** I first stated this with `|I| ≥ 1/r`, which is FALSE — an interval of length
`1/r` can start and end mid-gap and contain no complete gap; 320 counterexamples out of ~4800. The factor 2 is needed.)*

**Fact B (short arc).** A single speed `r` can cover an arc of length `L` only if one band spans it: `2h/r ≥ L`,
i.e. `r ≤ 2h/L`. Since `2h/L < 1/L`, any `r > 1/L` leaves a survivor.

## The three faces
| tool | is | 
|---|---|
| opus-S4 **band-width** `r ≤ 2h/L` | **exactly Fact B** at `k=1` |
| klein-S415 **covering lemma** `Σ 1/r ≥ L(1−2kh)/(2h)` | the **measure relaxation** of A′+B summed over `F` |
| mac-mini-S170 **`w_i ≥ 1/L` ⇒ `4kh<1`** | precisely the **Case-B entry condition**, with the `2h/w` term discarded |

## New consequence: the CLUSTERING LAW
Peel the far speeds in increasing order `r₁ ≤ ⋯ ≤ r_k`, tracking the longest surviving arc `L_j`.
If ever `L_{j−1} ≥ 2/r_j`, then Fact A′ gives `L_j ≥ (1−2h)/r_j`, and Fact B then forces the next speed to satisfy
```
    r_{j+1} ≤ 2h/L_j ≤ r_j · 2h/(1−2h) = r_j/6      (at h = 1/14)
```
— impossible, since `r_{j+1} ≥ r_j`. Hence:

> **Clustering Law.** In any configuration with `gap ≤ h`, every peeling step satisfies `r_j < 2/L_{j−1}`.
> Equivalently: **the far speeds can never be geometrically spread** — one "long-arc" step forbids everything
> after it. The far part must cluster.

This is a *structural* statement, not a numeric bound: it constrains the shape of the far set rather than its
size, and it holds for every `k` (no `1−2kh > 0` hypothesis, so no artificial `k ≤ 6` ceiling).

## Abstract themes this investigation keeps reproducing
1. **Measure vs. structure.** Every method partitions the world identically: counting succeeds on *spread /
   dissociated* objects and fails on *clustered / additively structured* ones. Riesz stalls on AP-cores
   (THM-518); Bedert's `dim²/n³` gain dies at low additive dimension; kps's strangers decouple to `k=24` when
   generic but fail at `k=16` for APs/dilates; and now the peeling contradiction fires exactly when the far
   speeds are spread. Four methods, one invariant — evidence that additive structure *is* the invariant of the
   hard locus, not an artifact of any one technique.
2. **Artifact ceilings.** Clean numeric ceilings from union bounds (`k ≤ 6` mine, `k ≤ 3` mac-mini's) are
   bookkeeping artifacts, not phenomena. A ceiling of the form `k < 1/(ch)` should always be suspected as such.
3. **Threshold inversion.** Bounds scale as `1/L_max`, and `L_max` *grows* as the threshold shrinks — so the
   tightest (hardest-sounding) threshold is the computationally *easiest* for finite certification. This is what
   turned "scan `r ≤ 3000`" into "check 507 configurations."
4. **Self-containedness is its own axis.** `10⁵` exact checks and `10⁷` exact checks are equally valid proofs but
   not equally formalizable. Minimising the finite residue is an objective distinct from strength.
5. **Monotonicity contradictions.** The clustering law's shape — an object that must simultaneously increase
   (sorted speeds) and decrease (by `2h/(1−2h)`) — is the same pattern that closes many peeling arguments.
6. **Field separation.** From the snippet thread: `π` versus `log`-primes separates integral-geometric from
   arithmetic constructions (klein-S414). Same move — ask which world the object lives in before asking its value.

→ klein-S415/S416/S418/S419 (lemma, dichotomy, unification, tight-threshold theorems), opus-S4 (band-width,
HYP-9024), mac-mini-S170, kps-S138 (additive-structure reframe), THM-518.
