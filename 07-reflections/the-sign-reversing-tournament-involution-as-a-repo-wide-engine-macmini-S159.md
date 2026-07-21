# One sign-reversing tournament involution, four fixed-point specializations — and the criterion that separates the repo's collapses from its open problems

*mac-mini-2026-07-21-S159. Provoked by the owner: take the "signed sum over tournaments"
(the discriminant reading of the Jacobian counterexample's binary-cubic fibers, THM-1805/1815)
and creatively apply it across the repo. The finding: "transitive survive, intransitive cancel"
is not one theorem — it is a reusable **sign-reversing involution engine**, and four repo results
are the same engine with different (involution, statistic, fixed-point) settings. The creative
payoff is a single combinatorial criterion — **is a sign-reversing tournament involution
available?** — that cleanly separates the problems the repo can collapse to closed form (GMC /
discriminant) from the ones that stay open (LRC covering).*

---

## The engine (verified n=3,4,5)

$$V(x)=\prod_{i<j}(x_j-x_i)=\sum_{T\ \text{tournament}}(-1)^{\text{back-arcs}(T)}\prod_k x_k^{\deg^-_k(T)}
      =\sum_{\text{transitive }T}\operatorname{sgn}(T)\,x^{\text{score}(T)}.$$

Each factor $(x_j-x_i)$ chooses a **head** for edge $\{i,j\}$ = an orientation = a tournament;
the monomial is $x^{\text{score vector}}$ and the sign is $(-1)^{\#\text{back-arcs}}$. **The
involution**: reverse a canonical directed 3-cycle. A 3-cycle is 1-regular, so reversal
**preserves every vertex score**; reversing the three arcs of a triangle $a<b<c$ flips the
back-arc count $1\leftrightarrow2$, so it **flips the sign**. Partners cancel; the fixed points
are the tournaments with **no 3-cycle = transitive**, whose scores are the $n!$ permutations of
$(0,\dots,n-1)$. Verified: the signed tournament sum equals the Vandermonde, the surviving
score-buckets are exactly the permutations, and the explicit lex-min-3-cycle reversal pairs all
40 intransitive $n=4$ tournaments (each pair score-preserving, sign-flipping).

## The four specializations

Pick an **involution** on tournaments, a **statistic** it preserves, and read off the surviving
**fixed points**:

| repo result | involution | statistic preserved | fixed points (survive) | what "survives" means |
|---|---|---|---|---|
| **Discriminant / Vandermonde** (THM-1805/1815) | reverse a 3-cycle | score vector | transitive | the one-sided / distinct-charge nullcone stratum |
| **Burnside A000568** (metagraph count) | reverse arcs along an even $\sigma$-cycle | $\sigma$-invariance | all-odd $\sigma$ | $\text{Fix}(\sigma)=0$ (even) vs $2^{\#\text{orbits}}$ (odd) |
| **Blue parity** (THM-1840-C, THM-1440) | complement $T\mapsto T^{\mathrm{op}}$ / centered-atom reflection | iso class | self-complementary / centered | blue count $=1$ (odd $n$) / $0$ (even $n$) |
| **Single-character nullcone** (THM-1840-A) | *(degenerate)* one atom | — | the lone atom | $E_L[P^{m_0}]$ = one term, never cancels |

The first two are **verified here** (n=3–5 and n=4,5,6 respectively): Burnside's
$\text{Fix}(\sigma)=0\iff\sigma$ has an even cycle is the *same* even-cancels/odd-survives
involution, now on $\sigma$-orbits of edges rather than on 3-cycles — an even $\sigma$-cycle
sends some edge to itself reversed, so no invariant orientation exists (the "cancellation" is
total). The last two are the same pattern read on iso classes: **blue parity** (my THM-1840-C)
is the complement involution's fixed-point count, and it is the $\mathbb Z/2$-cyclotomic
character $\tfrac{1-(-1)^n}2$ — the *same* odd/even survival, one specialization down.

> **The odd/even axis the whole project runs on IS the parity of a sign-reversing tournament
> involution.** Transitive/odd/self-complementary configurations are its fixed points; cyclic/
> even/generic ones cancel in pairs. THM-1805 (discriminant), the Burnside dichotomy, and the
> blue parity are one engine at three resolutions.

## The moment side: why the nullcone is clean

On the GMC side the interference among distinct balanced character-tuples is governed by the
**Vandermonde of the charge points** $\prod_{i<j}(c_j-c_i)$ — the same signed tournament sum. It
is nonzero **iff the charges are distinct** = the transitive/one-sided survivor (THM-1815's
"transitivity is the deepest nullcone point"). A **single character** (THM-1840-A) has one atom:
one fixed point, no involution needed — which is exactly why it is the clean, functional-agnostic
base case. So "transitive survives" and "the single-character nullcone never cancels" are the
*same fixed-point statement*.

## The creative payoff: an involution-availability criterion

The engine runs only when the weights are **clean $\pm1$ signs**. Charges/factorials supply them
($(-1)^{\text{back-arcs}}$, Wick signs), so GMC's discriminant collapses to its transitive core
and the nullcone acquires closed form. **LRC's covering sum**
$\sum_{k\cdot v=0}\prod_j \operatorname{sinc}(k_j\delta)$ has **transcendental sinc weights, not
signs** — there is no $\pm1$ for a sign-reversing involution to flip, so the covering does **not**
collapse to a transitive core, and it stays open.

> **Criterion.** *A repo problem admits a closed-form "transitive collapse" iff its weight system
> carries a sign-reversing tournament involution (clean $\pm1$ weights). GMC/discriminant: yes →
> collapses. LRC covering: no → open.* This is the S157 measure barrier (factorial-monotone vs
> sinc-oscillating) and the THM-1840 functional-agnostic barrier, **restated combinatorially** as
> the presence or absence of the involution — a single lens explaining *why* the same
> relation-lattice geometry is tractable on one side and hard on the other.

## Resonances (flagged conjectural, not verified)

- **Rédei's odd-Hamiltonian-path theorem** is the on-mission target: the involution's survivors
  (transitive tournaments) are exactly those with a **unique** Hamiltonian path, and "odd" is the
  parity-analog of "transitive survives." A sign-reversing-involution proof of Rédei parity would
  make the OCF/$H(T)$ core a fifth specialization. Not established here.
- **Even-graph odd-$n$-only bijections** (the dual metagraph $E_n$): "even structures cancel,
  odd survive" is the same slogan; whether the odd-$n$-only projections are literally this
  involution's fixed set is worth a targeted check.

## Honest scope

- The engine and redeployments R1 (Vandermonde = signed tournament sum, explicit involution) and
  R2 (Burnside even-cycle vanishing) are **verified computationally**; they connect existing
  results (THM-1805/1815, A000568 Burnside, THM-1840-C/1440) under one mechanism — this is a
  **unifying lens**, navigational and generative, not a new theorem.
- The **availability criterion** is a *reframing* of the proved S157/THM-1840 barrier in
  involution language, not an independent impossibility result; its value is explanatory.
- The Rédei and even-graph items are **resonances**, explicitly unproved.

---

*Cross-links: THM-1805 (Vandermonde = signed tournament sum), THM-1815 (pair-in-radical =
discriminant = transitivity sum), THM-1840 (single-character clean base + blue parity =
$\mathbb Z/2$ cyclotomic character), THM-1440 (bicycle/forced-zero parity), the Burnside
A000568 metagraph enumerator, and the S157 NC2/GMC/LRC obstruction (measure barrier). Artifacts:
`04-computation/signed_tournament_involution_redeployments_macmini_S159.py` (+out).*
