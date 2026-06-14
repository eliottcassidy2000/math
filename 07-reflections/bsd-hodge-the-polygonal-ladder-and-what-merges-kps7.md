# BSD, Hodge, the polygonal ladder — what merges, what's metaphor, and why Pollock isn't pentagonal

**Source:** kind-pasteur-2026-06-14-S7. A free-associative dispatch (Fermat
polygonal, Goldbach/Lemoine, the line segment, RH's 1/2, BSD, Hodge). The job: keep
the creative threads, but separate theorem from metaphor honestly — the literature
mandate applies hard here. Facts web-verified; repo bridges pulled.

## The honest ledger first

| claim | status |
|---|---|
| Fermat polygonal: every n is a sum of ≤ k k-gonal numbers (k≥3) | THEOREM (Cauchy 1813; Gauss k=3, Lagrange k=4) |
| 2-gonal numbers ARE the integers; the "line segment"/digon is the degree-1 base | TRUE (`P_2(m)=m`) |
| weak/ternary Goldbach: odd n≥7 is a sum of 3 primes | THEOREM (Helfgott 2013) |
| Lemoine: odd n≥7 is `p + 2q`, p,q prime (a prime + a doubled prime) | OPEN (verified to ~10^10) |
| BSD: `rank E(Q) = ord_{s=1} L(E,s)` | OPEN; proved for rank ≤ 1 (Gross–Zagier–Kolyvagin) |
| Hodge: every Hodge class is a ℚ-combination of algebraic cycle classes | OPEN |
| modularity: every elliptic curve /ℚ is modular | THEOREM (Wiles semistable; BCDT all) |
| `Δ = η^24` is the weight-12 cusp form; its coefficients are Ramanujan `τ` | THEOREM |
| **Lemoine `2q`-doubling ⟺ RH's `1/2`** | **NOT a theorem — metaphor only.** No published rigorous link. The "2" of doubling and the "1/2" of the critical line are not connected. Flagged so we don't confabulate. |

With that fence built, here is what is real.

## The real answer to "how is Pollock different from pentagonal": the DEGREE jump

The `k`-gonal numbers are **degree-2** polynomials in their index; sums of them are
a *quadratic* Waring problem, solved exactly by Cauchy (via sums of three squares),
with the clean bound `k`. Pollock's tetrahedral/octahedral numbers are **degree-3**;
sums of them are a *cubic* Waring problem — the world of Waring-for-cubes, the circle
method and Linnik's ternary-form reduction, with messy bounds. Computed (the
figurate degree ladder, `figurate_degree_ladder_kps7`):

```
DEGREE 1  line segment / 2-gon  P_2(m)=m = the integers          bound 1 (trivial)
DEGREE 2  polygons (quadratic): triangular 3, square 4, pentagonal 5, hexagonal 6,
          heptagonal 7, octagonal 8   -> bound = k EXACTLY (Cauchy/Gauss/Lagrange)
DEGREE 3  polytopes (cubic):   tetrahedral 5, octahedral 7, cube 9  -> MESSY
          (tetrahedral conj 5 / proved 8; icosahedral conj 13 / actually 15)
```

So **Pollock ≠ pentagonal = degree 3 ≠ degree 2.** The geometric *dimension* of the
figure (1 segment, 2 polygon, 3 polytope) is exactly the polynomial *degree* is
exactly the *difficulty class*: degree 2 lives in the clean sums-of-squares world
(exact, formula `k`); degree 3 lives in the circle-method world (asymptotic + finite
check, no formula, conjectures that are *wrong by a finite exceptional set* — Pollock
got 13 for icosahedral; the truth is 15). The user's instinct that the segment is
the base of the ladder is right: the digon is the degree-1 floor, below which
"figurate" degenerates to counting.

This ladder is the same one the repo's OCF invariants climb: `c3, c5` are low-degree
power sums (`tr A^k`, spectral, clean — THM-118/498), `H` is higher (non-spectral,
the conflict-graph layer — THM-499). Degree is the universal difficulty knob; the
"efficiency-becomes-proof" thread (kps4/5) is the same statement on the tournament
side.

## The one real gem in the prime musings: the doubled prime does the work of two

Lemoine (`n = p + 2q`) **implies** ternary Goldbach (`n = p + q + q`): the doubled
prime `2q = q + q` lets a single doubled prime stand in for two of the three primes.
So Lemoine is the *finer* conjecture and proved-ternary-Goldbach is its corollary —
the user's "a doubled prime does the work of two of them" is exactly this, and it is
a clean fact, not a metaphor. The repo already runs this `2q` motif on three fronts
(doubled-primes-as-the-parity-hinge S546o): Lemoine's parity hinge, the Burnside
antipodal involution that kills `A000568` fixed points, and the LRC apex halving
`n → n/2` — and **rank-one LRC configurations are exactly `n = 2q`, q prime.** The
*doubling* (pairing / parity) is structural and real across the repo. The *leap to
RH's 1/2* is not — though the repo's recurring "2-adic seam" (the dyadic valuation,
the band-0 layer, the LRC trivial bound `1/(2n)` vs the target `1/n` = a literal
doubling of the gap) is why the metaphor *feels* right. Resonance, not theorem.

## What actually merges with BSD and Hodge: the special-value lens (and a concrete bridge)

The deep, articulable unification — a *lens*, not a theorem-merge — is **special
values of generating functions encode arithmetic.** The repo's *fugacity* `x` is its
"L-function variable `s`": `I(Ω, x=1)` = Zeckendorf/Fibonacci, `I(Ω, x=2)` = the OCF
`H`, with a Binet interpolation across `x` (eureka-zeckendorf-simplex-cuboid;
fugacity axis `x ∈ {−1,0,1,2,6}`). BSD reads the *rank* off `L(E,s)` at the special
point `s=1`; the OCF reads `H` off `I(Ω,x)` at `x=2`; my LRC singular series
`L(S) = lim D(q,S)/q` (THM-501) is a special-value/density that vanishes exactly at
tightness. Same architecture: *arithmetic/combinatorial data is the value of a
generating function at a distinguished point.*

And there is a **concrete bridge to BSD's world**, which the recon flagged as the
repo's missing modularity link — now supplied. The code discriminant
`P_24 = η^24 = Δ` (THM-489) is not just *a* modular form: its coefficients are the
**Ramanujan `τ` function** (verified: `1, −24, 252, −1472, 4830, …`, multiplicative,
satisfying Deligne's bound `|τ(p)| ≤ 2p^{11/2}` — the proved RH-analog for `Δ`;
`eta24_ramanujan_tau_bridge_kps7`). Its Dirichlet series `L(Δ,s) = Σ τ(n) n^{−s}` is
the **prototype modular L-function** — Euler product, functional equation,
Ramanujan–Petersson — and by modularity it is the *same class of object* as BSD's
elliptic-curve L-functions. So the repo, via its pentagonal hub, already computes
with a genuine modular L-function generator. That is the honest, concrete touchpoint:
not "the LRC implies BSD," but "the repo's `Δ = η^24` lives in the modular-L-function
ecosystem that BSD inhabits."

**Hodge, honestly = realizability.** Hodge asks *which cohomology classes are
algebraic* — realized by an actual subvariety. The repo's structural analog is
*which invariant values are realized by an actual tournament*: the gap-free-spectrum
work (THM-498: which `c5` values occur; THM-501: which deficit densities `L(S)`
occur). Both are "is this numerical class realized by a structured geometric object,
or is there an obstruction?" The forbidden values (`c5=10`, `H∈{7,21}`, the LRC
exceptional set) are the combinatorial analog of a Hodge class that *fails* to be
algebraic. Analogy, clearly labeled — but a productive one: it says the repo's
"which values are realizable" program is a baby Hodge problem for tournament
invariants.

## The synthesis

The user's web — polygons, primes, doubling, RH, BSD, Hodge — is unified by one lens
the repo already runs: **a generating function's special value encodes arithmetic,
and which values are realizable is the hard (Hodge-flavored) question.** What I can
stand behind rigorously: (i) Pollock ≠ pentagonal is the degree-3-vs-2 jump
(computed); (ii) Lemoine ⟹ weak Goldbach via `2q = q+q` (the doubled prime does the
work of two); (iii) `Δ = η^24` generates the Ramanujan-`τ` modular L-function, the
repo's concrete handle on the BSD/modular world (verified). What is metaphor —
labeled as such — : the Lemoine-`2q` ↔ RH-`1/2` link, and "BSD/Hodge about the LRC."
Productive metaphors point at the special-value/realizability lens; they are not
theorems, and saying so is the same discipline as the spectral-reframe *boundary*
(kps5): knowing where a bridge stops is as valuable as knowing it carries.

Cross-links: THM-489 (`Δ=η^24`, now → Ramanujan `τ` L-function),
THM-498/499/501 (the degree/spectral hierarchy = the figurate degree ladder's
tournament avatar; the realizability/Hodge analog), THM-118 (`c_k=tr A^k`, degree),
HYP-614 (Dedekind regulator, the Q(√5) special value),
[[eureka-zeckendorf-simplex-cuboid]] (the fugacity-as-`s` axis),
[[doubled-primes-as-the-parity-hinge-cycles-numbers-and-lrc-channels-s546o]],
[[the-lrc-singular-series-kps6]], [[pollock-as-the-bounded-arity-currency-and-the-cycle-spectrum-onset-kps3]].
