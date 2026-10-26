---
source: oracle-2026-06-02-S559o
status: synthesis + new cross-link (the user's A+B=C summand graph = the LRC pinch denominator; addition-shadow=pinch, multiplication-shadow=sieve, they meet on the rationals)
tags:
  - lonely-runner
  - summand-graph
  - pinch
  - sieve
  - addition
  - multiplication
  - parity
  - tournament
---

# The Summand Graph Node IS the LRC Pinch Denominator — Addition and Multiplication as the Two Witness Shadows

The user recalled an additive tournament-like construction on ℕ — *nodes `A,B,C` with
arcs `(A→C),(B→C)` iff `A+B=C`* — and asked to connect it to the repo's even/odd and
addition/multiplication threads. It is a recorded object, and it links **directly** to
the live LRC@14 pinch frontier in a way nobody had written down.

## The construction is the SUMMAND GRAPH (oracle-2026-05-03)

`07-reflections/summand-graph-fermat-zeckendorf.md` defines exactly this: the digraph on
ℕ with `a→n` and `b→n` whenever `a+b=n`, `a≠b`. Known facts:
- node `n`'s **incoming-pair-count** is `⌊(n−1)/2⌋` — **parity-split**: odd `n` →
  `(n−1)/2` pairs; even `n` → `n/2−1` (the midpoint pair `n/2+n/2` is excluded).
- **Section 4 there:** "the summand graph IS the tournament staircase" — node `n` ↔
  `δ_{n−2}`, depth-`k` multiplicity `2^{k−2}` mirrors the tournament binary tree, H-values
  fall on Zeckendorf/Fibonacci sums, and the forbidden `{7,21}` recur as nodes.

## The new bridge: node `C` = the LRC pinch denominator

opus-S557 (the pinch lemma): for a set with `M(S)<1/2`, the loneliness radius is
`M=r/s` where the binding pair `(v_a,v_b)` **straddles** the observer and
`s=(v_a+v_b)/gcd` is the **reduced pair-sum**, attained at the pinch time `t*=m/(v_a+v_b)`.
oracle-S555o: the LRC@14 pinch pigeonhole runs over "the `⌊n/2⌋` pairs summing to `n`."
Those pairs `{(a, n−a)}` are **exactly the incoming pairs of node `n` in the summand
graph.** Hence:

> **Summand-graph node `C` = the LRC pinch denominator `s = v_a+v_b`. Its incoming pairs
> are the pinch-pairs; its pair-count `⌊(C−1)/2⌋` is the number of pinch times at that
> denominator.** The user's additive tournament on ℕ *indexes the LRC pinch times.*

The parity split of the pair-count is the parity of the pinch: an **even** denominator
`C=a+b` loses the midpoint pinch (the apex pair `C/2+C/2`) — the very apex/zero-divisor
that obstructs `n=14=2·7` (opus-S559, HYP-2063); an **odd** `C` keeps all `(C−1)/2`.

## Addition and multiplication = the two witness shadows (codex-S365 ⊕ THM-369)

`natural-operation-digraphs-and-product-sum-s365.md` (codex-2026-05-31) collapses the two
two-input operation systems to simple digraphs:
- **ADDITION shadow:** `a→z, b→z iff a+b=z` ⟹ `x→z iff x<z` = the **transitive
  tournament** (complete, acyclic). This is the **pinch** family: a witness time
  `t=m/(a+b)` indexed by a summand pair.
- **MULTIPLICATION shadow:** `xy=z` ⟹ `x→z iff x|z` = the **divisibility DAG** (sparse
  incomplete tournament). This is the **sieve** (THM-369: `t=1/q` is lonely iff **no speed
  is divisible by `q`** — pure divisibility).

The punchline ties to oracle-S555o, which proved **the rational pinch IS the sieve**
(`t=1/(a+b)` lonely ⟺ `(a+b)∤` every speed). In operation-shadow language:

> **The addition-shadow (pinch) and the multiplication-shadow (sieve) COINCIDE on the
> rationals.** A rational witness `t=p/q` is both "a summand-pair pinch at denominator
> `q`" and "a divisibility test mod `q`." The **fine pinch (`q>n`, S556/S556o) is exactly
> where the addition-shadow outruns the divisibility-shadow** — where the additive
> (summand) structure has witnesses the multiplicative (divisibility) structure cannot
> see. That is the open core, restated in the user's two-operation language.

This is the same `+` vs `×` spine as `multiplication-is-repeated-addition-...-s548`
(the hyperoperation tower `succ→+→×`, `∏c_i=exp(−Σlog c_i)`): the pinch lives on the
`+` rung, the sieve on the `×` rung, and `log` welds them — the rationals are where the
weld is exact.

## Parity, once more (the even/odd thread)

- Summand pair-count parity (above) = pinch-count parity; even `C` drops the apex pair.
- Burnside (canon): even cycle types fix `0` tournaments, all-odd fix `2^{pairs}` — so the
  doubling/`×2` is the "killer," and `A000568` is an all-odd sum.
- The **doubled prime `2q`** is "the parity-completion of the prime additive basis"
  (`doubled-primes-as-the-parity-hinge-...-s546o`); `14=2·7` is our hinge. The apex
  (mult of 14) is the even node whose summand-pinch lacks the midpoint — and the even-fold
  measure reduction (S558o, HYP-2065) is the `×2` acting on the summand structure.

## Verdict / next
- **New cross-link (recorded):** the user's `A+B=C` summand graph indexes the LRC pinch
  times (node `C` = pinch denominator `s=a+b`, incoming pairs = pinch pairs); addition- and
  multiplication-shadows are the pinch and the sieve, equal on the rationals (S555o), split
  by the fine pinch (S556o).
- Concrete next: (1) use summand-graph node structure (incoming-pair lattice of `14`, and
  of fine denominators `q∈(14,Cn]`) to organize the **fine-regime pinch pigeonhole** (the
  S555o handoff) — the additive shadow is the right index set; (2) the apex = even node's
  missing midpoint pair — phrase HYP-2063's zero-divisor purely in summand-graph terms;
  (3) Zeckendorf/Fibonacci summand chains (oracle-2026-05-05) as a *sparse* pinch family.

## Artifacts (prior work surfaced)
```
07-reflections/summand-graph-fermat-zeckendorf.md            (oracle-2026-05-03; the A+B=C graph)
07-reflections/lucas-summand-graph-zeckendorf-geometry.md    (oracle-2026-05-05; Fibonacci/Lucas chains)
07-reflections/natural-operation-digraphs-and-product-sum-s365.md (codex-2026-05-31; +/× shadows)
07-reflections/multiplication-is-repeated-addition-...-s548.md    (oracle-2026-06-01; hyperoperation tower)
07-reflections/lrc-n14-the-exact-moments-pinch-pair-and-r-over-s-radius-s557.md (opus; pinch=r/s)
07-reflections/the-pinch-time-pigeonhole-is-the-sieve-and-the-fine-regime-salvage-s555o.md (pinch=sieve)
```
Related: THM-369 (divisibility sieve = the × shadow), HYP-2063 (apex zero-divisor),
HYP-2065 (even-fold = ×2 on the summand structure), S548 (+/× tower), S546o (parity hinge).
