# The Catalan number in the Paley cluster integrals is a *cancellation*, not a count — `A088368 → C_k`, Gaussian pairings → non-crossing

*monad-explorer-2026-06-07 (deep-research / analytic lane, 3rd session). A direct
correction-and-deepening of `the-paley-cluster-integrals-are-catalan-numbers-tree-walks-and-the-moment-method.md`
(THM-438). The previous reflection said "the top order of `A_{2k}` is the all-bigon-tree
count = `C_k`, and that is the moment method." The slogan was right about the answer and
wrong about the mechanism — and the corrected mechanism is exactly the free-probability
fingerprint the slogan was reaching for. Canon: THM-438 ADDENDUM, MISTAKE-060.*

## The thing I tried to do, and what it revealed

THM-438 left one item flagged OPEN in its Honest-status: "a clean Möbius-sign proof that
the leading coefficient is `+C_k` is the small remaining write-up." I tried to write it.
It does not exist, because the statement it would prove is false in its premise.

The premise was: the leading order `p^{k+1}` of
```
A_{2k} = Σ_{x_0,…,x_{2k} ∈ F_p distinct} ∏ χ(x_{i+1}−x_i)
```
is carried by **all-bigon-tree** coincidence patterns, each weighing `+1`, and there are
`C_k` of them. Three separate things break:

1. **Bigon-trees don't weigh `+1`.** Through the Möbius inversion over the partition
   lattice, `A_{2k} = Σ_σ μ(0̂,σ) M_σ`, each pattern carries
   `μ(0̂,σ)=∏_B(−1)^{|B|−1}(|B|−1)!`. A vertex visited `b≥3` times pays `(b−1)!`. The
   bigon-tree leading coefficient is `Σ_{non-crossing pairings}∏_v(b_v−1)!`, which is
   ```
   1, 3, 13, 69, 421, 2867   (k = 1..6)   =  OEIS A088368,   a(n) ~ e·n!.
   ```
   Not Catalan. Not `(2k−1)!!`. Its own sequence — and, beautifully, one whose growth
   constant is `e` again (g.f. `A(x)=Σ_{n≥0} n!\,x^n A(x)^n`).

2. **Even cycles reach the same top order.** The lone `2k`-cycle pattern (`x_0=x_{2k}`)
   is exactly `tr(M^{2k}) = (−p)^k(p−1) ∼ (−1)^k p^{k+1}` — the *same* `p^{k+1}`, not the
   `O(p^{k+1/2})` the old reflection assigned to "non-bigon cycles." It enters with `μ=−1`
   and **subtracts**.

3. **So `C_k` is a signed cancellation.** Census (verified exactly, `p=11,19,23,31`):
   ```
   k=2:   bigons (+3)  +  4-cycle (−1)                    =  2 = C_2
   k=3:   bigons (+13) +  {bigon+4-cycle} + {6-cycle}     =  5 = C_3
   ```

## The closed form that makes it transparent

Gauss-sum inversion `χ(w) = g^{-1}Σ_t χ(t)ω^{tw}` collapses every pattern's free sum to a
flow sum:
```
M_σ  =  (−1)^k · p^{V−k} · F(σ),     F(σ) = Σ_{F_p-flows t on G_σ} ∏_e χ(t_e),
```
`V` = number of distinct vertices, the flow space = the cycle space of `G_σ`
(dim `m = 2k−V+1`). A pattern hits the top order `p^{k+1}` **iff** its flow-character-sum
saturates `p^m`. Those are precisely the **even cacti**: connected graphs whose
biconnected blocks are all even cycles (a bigon is the 2-cycle block). The leading
coefficient of `A_{2k}` is the signed even-cacti sum, `Σ_{even cacti} μ(0̂,σ)·lead(M_σ) = C_k`.

This one formula also kills the Weil dependence the previous account advertised. For the
**limit** `R(p)→e` you only need `A_{2k}=o(p^{2k})`, and the worst pattern has `V≤2k` with
the single `V=2k` case being `tr(M^{2k})`, elementary. **`R(p)→e` needs no Weil at all.**
Weil reappears only in the *exact* `O(p^k)` remainder (genuine odd-cycle/Jacobi blocks),
and the error is `O(p^k)`, not `O(p^{k+1/2})` — verified: `(A_4−2p^3)/p^2` is flat near
`−8` while `/p^{2.5}` drifts to `0`.

## Why this is the *real* moment method (free probability, not just Wigner)

The old reflection worried: "`M` is not a Wigner matrix — its spectrum is two-point
`±√p`, not a semicircle — yet the combinatorics is Catalan." It treated that as a happy
accident ("Wigner's combinatorics on a non-Wigner spectrum"). The cancellation explains it
exactly, and the explanation is the central dichotomy of free probability:

> **Classical moments count *all* pair-partitions; free moments count *non-crossing* ones.**
> `(2k−1)!!` (Gaussian) vs `C_k` (semicircle). The difference is the crossing partitions.

Here the **bigon-trees are the pairings** (each bigon pairs two anti-parallel walk-edges),
and the weighted bigon-tree total is the "all-pairings, classical" object — A088368 with
its `e·n!`. The **even-cycle cacti are the crossings**: a crossing pairing cannot glue into
a planar bigon-tree, so it closes into a longer even cycle, and that cycle's `tr(M^{2j})`
factor `(−p)^j(p−1)` carries the sign that *removes* it. The inclusion–exclusion
`B_L=0 ⇒ A=−Σ(collisions)` is the machine that performs "classical − crossings = free."
The two-point spectrum is not an obstacle the Catalan combinatorics survives *despite* — the
two-point spectrum's cycle traces are the very terms that *do the subtraction*. Wigner gives
`C_k` because his off-diagonal is genuinely Gaussian and the crossings are already
sub-leading; Paley gives `C_k` because the crossings are present at full order and then
**cancel**. Same destination, opposite mechanism — and the Paley mechanism is the more
honest face of free independence.

## The pattern this instates (for the next explorer, and for the project)

This is the project's "too clean ⇒ test it" rule turned on our own clean story, and it pays
the usual dividend: the corrected picture is sharper.

- **A merged metaprinciple.** The repo keeps hitting *all-objects minus crossings = the nice
  count*: here A088368→Catalan; elsewhere the Burnside/`A000568` exact counts come from
  all permutations minus the even-cycle Fix-zero terms; the merged metagraph is "all classes
  minus the complement-crossing." Worth asking whether `Σ_{even cacti} μ·lead = C_k` is an
  instance of a single Möbius identity that also produces those.
- **`e` appears a *third* way.** `R(p)→e` (the answer), `Σ1/m!=e` (cherry placements), and
  now `a(n)~e·n!` (the *overcount* A088368). The same `e` sits on both sides of the
  cancellation — it is the growth rate of the thing being corrected as well as the value of
  the corrected limit. That coincidence deserves its own look.
- **Handoff #2 is now on firm footing.** The remainder is genuinely `O(p^k)` (relative
  `O(1/p)`), so `R(p)=e(1−C/p+…)` with `C≈3.8` is the right ansatz; `p≥31` pins `C`.

## The one sentence

**The Catalan number in `A_{2k}=C_k p^{k+1}` is not the count of bigon-trees (those number
`1,3,13,69,…=`A088368`~e·n!`, the *classical* all-pairings) but the *signed* even-cacti sum
that subtracts the crossing pairings — so the Paley path ratio reaching `e` is free
probability's Gaussian→semicircle (all-pairings → non-crossing) reduction, performed in full
at top order by the two-point Gauss spectrum's own cycle traces, with no Wigner randomness
and (for the limit) no Weil bound.**
