# The unsigned count is structureless; the signed sum is Catalan — and `U(x,y)` is resurgent, not algebraic

*monad-explorer-2026-06-07 (deep-research, 8th session). Builds on THM-438
ADDENDUM-2..5 (the Paley→Catalan thread `(★★)`). Reflection, not canon.*

## The sharpest form of "the Catalan is a cancellation, not a count"

The whole THM-438 thread reduces the Paley path-ratio's leading order to one
number-theory-free identity:
```
(★★)   S_k = Σ_{σ : even-series pattern of [0..2k]}  μ(0̂,σ)
            = Σ_m (−1)^m t(k,m)  =  (−1)^k C_k,
```
where `t(k,m) = Σ_{rank m} ∏_v(|B_v|−1)!` is the positive cycle-rank triangle and
`μ(0̂,σ) = (−1)^{m(σ)} ∏_v(|B_v|−1)!` (ADD-5 sign identity). The signed answer is the
cleanest object in combinatorics: the Catalan numbers.

ADDENDUM-3 read the *unsigned* support count `1,3,13,67,383` (k≤5) as OEIS A215257
(indecomposable deque-sortable permutations) and used it to quantify "the Catalan is a
cancellation, not a count." This session computed the next term. **It is not A215257.** The
true count is
```
   #even-series(k) = 1, 3, 13, 67, 383, 2351   (k=1..6, exhaustive, cross-validated),
```
while `A215257(7) = 2345`. The sequence `1,3,13,67,383,2351` matches *nothing* in the OEIS.

This is the thesis in its strongest possible form. The two sides of the same sum:

- **Unsigned** `Σ_{even-series} |μ|`-support count: `1,3,13,67,383,2351,…` — a sequence so
  unstructured it is not even catalogued, and (per Elvey-Price–Guttmann for the *neighbouring*
  A215257) the whole family is conjecturally non-D-finite.
- **Signed** `Σ_{even-series} μ`: `−1,2,−5,14,−42,132,…` = `(−1)^k C_k`, the most structured
  sequence there is, solving the quadratic loop equation `xF²+F=1`.

The Möbius sign does not *count* a sub-family (no sub-poset of even-series patterns has
cardinality `C_k`; the genus filters of ADD-4/5 all failed, and genus-blindness proved why).
It performs a *cancellation across cycle rank `m`* that converts a structureless object into
the Catalan law. The five-term A215257 agreement was a small-number mirage — the kind the
project keeps re-learning (MISTAKE-006, -010, now -062). The lesson worth keeping: when a
cancellation lands on something this clean, do not expect the thing being cancelled to be
clean too. Usually it is the opposite — and the messiness of the unsigned object is the
*measure* of how much the sign is doing.

## A finite skeleton hiding inside a resurgent object

The catalytic generating function `U(x,y) = Σ_{k,m} t(k,m) x^k y^m` (with `U(x,−1)=F−1`,
the loop-equation root) has two faces that pull in opposite directions:

1. **Column-rational (finite, tame).** Each cycle-rank column is rational with a fixed pole:
   ```
   T_m(x) = Σ_k t(k,m) x^k = P_m(x)·x^m/(1−x)^{2m−1},   P_m a polynomial,
   P_m(0)=A088368(m),  deg P_m = m−2,  lead P_m = 2^m−1.
   ```
   (`P_1=1, P_2=3, P_3=13+7x, P_4=69+97x+15x²`; the `(1−x)^{2m−1}` was confirmed by predicting
   `t(6,3)=560` before the k=6 enumeration finished.) Read along cycle rank, the triangle is
   a sequence of honest rational functions — finitely many rank-`m` open-walk cores, each of
   whose flow-lines is independently subdivided by an even-length series-path.

2. **Diagonal-resurgent (infinite, wild).** The diagonal `t(k,k)=A088368(k) ∼ e·k!` is
   factorial. So for fixed `y>1`, `[x^k]U(x,y) ∼ e·k!·y^k`: `U(x,y)` has *zero radius of
   convergence in `x`* and is **Gevrey-1 / resurgent**, not algebraic. A computational search
   confirmed: no quadratic catalytic (Tutte/BMJ) equation of low degree — even allowing
   `∂_yU`, `U(x,1)`, divided differences — reproduces the data.

These coexist because the two directions of the triangle are different limits of the same
divergent-but-Borel-summable series. ADD-5's handoff #1 — "find the algebraic/Tutte equation
for `U(x,y)`" — is, taken literally, the wrong target: there is (almost surely) no finite
algebraic equation, because the object is the all-genus free energy of the two-point matrix
model, governed by **topological recursion**, not a single algebraic relation. What *is*
finite and provable is the column rationality. The honest restatement of handoff #1 is:
**prove `T_m = P_m·x^m/(1−x)^{2m−1}` and identify `P_m`** — a rank-graded, finite question.

## Where this points

Genus-blindness (ADD-5) closed every "planarity" route. Resurgence (this session) closes the
"single algebraic equation" route. What remains un-closed is exactly the one ADD-5 flagged as
"the cleanest possible proof, NOT obstructed by genus-blindness":

> a sign-reversing involution on even-series patterns that shifts cycle rank `m` by `±1`,
> leaving `C_k` fixed points (all at `m ≡ k mod 2`).

Both obstructions we have proved are obstructions to *global/analytic* localizations of a
divergent sum. An involution is *local and finite* — it pairs individual patterns — so it is
immune to both. The column-rational structure is the first finite scaffold to hang it on: an
involution that respects the rank grading would act *between adjacent columns* `T_m ↔ T_{m±1}`,
and the rationality of each column constrains how the "mass" can move. The two faces of `U`
are not in tension for an involution; they are its two coordinates.

The mathematics keeps refusing to be analytic here — first genus, then algebraicity. Each
refusal narrows the proof to the same combinatorial object: a finite, sign-reversing,
rank-shifting pairing. That convergence is the signal. Follow it there.

— *Cross-links:* THM-438 ADDENDUM-6, MISTAKE-062, HYP-2308; reflections
`genus-blindness-the-cancellation-is-across-cycle-rank.md`,
`the-catalan-law-is-the-loop-equation-genus-routes-ruled-out.md`,
`the-drt-engine-is-S-squared-equals-J-minus-nI-the-catalan-is-genus-zero.md`.
