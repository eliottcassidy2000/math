---
source: oracle-2026-05-31-S16
status: exploratory synthesis (long contemplative session on the n=14 frontier)
tags:
  - lonely-runner
  - n14-frontier
  - modulus-sieve
  - pushforward-measure
  - relation-lattice
  - torus-knot
  - view-obstruction
  - even-n
  - representations
---

# Lonely Runner at n = 14 — the First Open Case, and What It Represents

This is a deliberately abstract session. The task was to contemplate *what is
fundamentally going on* in the 14-runner case, from both the proof and the
disproof side, and to find new representations for the objects involved. It
builds on the repo's endpoint-protection / anti-Bohr-boundary program
(S356–S362) but tries to step outside it.

## 0. Why n = 14 is exactly the wall

External status (verified by web search this session):

- **8 runners**: Rosenfeld, computer-assisted (`arXiv:2509.14111`, 2025-09).
- **9, 10 runners**: Rosenfeld + sieve refinement (`arXiv:2511.22427`,
  `arXiv:2512.01912`).
- **11, 12, 13 runners**: Sungkawichai–Trakulthongchai (`arXiv:2604.23906`,
  2026-04-26).

So **n = 14 (k = 13 speeds, threshold 1/14) is the first open case** as of
2026-05. This is not an arbitrarily chosen frontier; it is *the* frontier.

The reason it is the wall is structural. A key engine of the recent proofs is a
polynomial-method theorem: every speed tuple `≡ (1,2,…,k) (mod p)` with
`gcd = 1` satisfies LRC **when `k+1` and `p > k²+k` are both odd primes**. For
n = 14 we have `k+1 = 14 = 2·7`, which is **even and composite**, so this tool
gives *nothing* for the near-initial-segment families at n = 14. The runners'
count being prime is exactly what powered n = 11 (`k+1=11`) and n = 13
(`k+1=13`); n = 12 and n = 14 are the composite holdouts, and n = 12 has now
fallen, leaving n = 14.

**Slogan:** n = 14 is the first case where the runner count is even *and*
composite *and* not yet beaten by an ad hoc argument. Its difficulty is the
difficulty of composite, even `n`.

## 1. The modulus-sieve reading and the "whack-a-mole" lemma

This is the cleanest new handle of the session, and it is elementary.

**Modulus-sieve lemma.** Let `V` be `k = n-1` nonzero integer speeds at
threshold `1/n`. For each integer `m` with `2 ≤ m ≤ n-1`:

> if **no** speed in `V` is divisible by `m`, then `t = 1/m` is lonely,
> because `‖v/m‖ ≥ 1/m > 1/n` for every `v` with `m ∤ v`.

For `m = n`: if no speed is divisible by `n`, then `t = 1/n` is a *boundary*
lonely witness (`‖v/n‖ ≥ 1/n` for all `v`).

**Corollary (sieve necessary condition).** A counterexample at level `n` must
contain, for **every** `m ∈ {2,…,n}`, a speed divisible by `m`. In particular it
must contain a multiple of every prime power `≤ n` and a multiple of `n` itself.

This is the baseline residue-covering used by the sieve proofs, but reading it
as a divisibility *cover* makes the n = 14 picture vivid. Specialize to
`m ∈ {2,…,14}`.

**The initial segment explained.** `V = {1,…,13}` covers every modulus
`m ∈ {2,…,13}` (the speed `m` divides itself) and misses **only** `m = 14`.
Computation (`lonely_runner_n14_frontier_s363.py`) confirms it is tight
(forbidden length 1, max gap 0) and its **six** surviving boundary witnesses are
exactly the units `a/14`, `gcd(a,14)=1`. So:

> the tightness of the extremal example localizes at the *single modulus scale
> it fails to cover*, namely 14, and the witnesses are the units there.

This is a much more satisfying "why" for the unit-boundary skeleton (S359–S362)
than "it's the Dirichlet equality case": the unit skeleton is the *residue of
the one uncovered sieve level*.

**Whack-a-mole.** To turn `{1,…,13}` into a counterexample we must cover `m=14`,
which forces a multiple of 14 into the set. But we have only 13 slots, so we
must drop a small speed `m₀`. Computation:

| edit | uncovered modulus becomes | result |
|------|---------------------------|--------|
| drop 13, add 14 | 13 | positive gap 1/308 opens (`t=1/13` lonely) |
| drop 8, add 14  | 8  | positive gap 2/147 opens (`t=1/8` lonely)  |
| drop 9, add 14  | 9  | positive gap 5/539 opens |
| drop 11, add 14 | 11 | positive gap 3/392 opens |

Covering modulus 14 **re-opens** a coarser modulus scale, and a *macroscopic*
gap appears: loneliness becomes *easier*, not harder. You cannot perturb the
extremal example into a counterexample; you fall off the knife-edge into easy
territory.

The only escape is to cover both 14 and the dropped scale with a *single larger*
speed (e.g. drop 13, add 182 = 2·7·13). Computation shows this works for the
sieve but the larger speed contributes only sparse, tiny intervals, so a
positive gap survives anyway. The extreme version: a single super-speed
`360360 = 2³·3²·5·7·11·13` is divisible by every `m ∈ {2,…,14}` at once, yet
the set `{1,…,12, 360360}` still has a (tiny) positive gap `1/420420`.

**The dilemma, exactly.** A counterexample at n = 14 must satisfy the
divisibility sieve, and there are only two ways to satisfy it:

- **(i) small speeds** — they cover macroscopic forbidden measure, but the
  initial segment misses `m = 14`, and any swap to cover 14 re-opens a coarse
  scale (a big gap);
- **(ii) large speeds** — they cover fine moduli, but their forbidden intervals
  are too sparse to keep the cover full (a tiny gap survives).

13 speeds cannot win both games. This is, I believe, the most honest one-line
description of *what is fundamentally going on* at n = 14.

## 2. The gap-floor hypothesis

The whack-a-mole evidence sharpens into a hypothesis, and — to my surprise — the
data supports a clean *general* form, not just an n = 14 form:

> **HYP (gap floor / no-multiple-of-n at full measure).** Any primitive speed
> set at level `n` whose forbidden union has **full measure**
> (`forbidden_length = 1`, i.e. tight or a cover) contains **no** speed
> divisible by `n`.
>
> n = 14 special case: every primitive 13-set containing a multiple of 14 has
> `max_gap > 0` (it is not even tight).

**Why this would settle n = 14 (and any n).** A full-cover counterexample is a
full-measure forbidden union, and by the modulus-sieve lemma it *must* contain a
speed divisible by `n`. The gap-floor hypothesis says a full-measure union has
*no* speed divisible by `n`. Contradiction. So the hypothesis at level `n`
implies **LRC at level `n`**.

**The evidence is striking.** Every catalogued tight instance avoids multiples
of `n`:

```
n=4: (1,2,3)                                  -- no multiple of 4
n=5: (1,2,3,4), (1,3,4,7)                     -- no multiple of 5
n=6: (1,2,3,4,5), (1,3,4,5,9)                 -- no multiple of 6
n=7: (1,2,3,4,5,6)                            -- no multiple of 7
n=8: (1,..,7), (1,4,5,6,7,11,13), (1,2,3,4,5,7,12) -- no multiple of 8
```

initial-segment *and* sporadic alike. And `lonely_runner_n14_gapfloor_s363.py`
(26019 primitive 13-sets *forced* to contain a multiple of 14: a complete
structured-perturbation sweep of 6019 plus 20000 randomized) found **no** tight
set; the minimum gap stayed strictly positive throughout — exactly
`1/728 ≈ 0.00137` at `{1,2,3,4,5,7,8,9,10,11,12,13,98}`, `98 = 2·7²`.

This is not a free lunch — a measure count alone does not forbid a multiple of
`n`, since `(n-2)·(2/n) ≥ 1 - 2/n` — and it should be checked against the
Goddyn–Wong classification of tight instances. But it is a remarkably sharp,
falsifiable target: **one full-measure primitive set containing a multiple of `n`
refutes it; a proof of it cracks the conjecture level by level.** To me this is
the most promising concrete lead the session produced.

## 3. The pushforward-measure identity: where the difficulty actually lives

Step outside combinatorics. The map `φ(t) = (v₁t,…,v_k t) mod 1` pushes Lebesgue
measure on `R/Z` to a measure `μ_V` on the torus `Tᵏ`. Its Fourier coefficients
are trivial to compute:

```
hat{μ_V}(a) = ∫₀¹ e^{2πi (a·v) t} dt = [ a·v = 0 ],   a ∈ Zᵏ.
```

So **`hat{μ_V}` is the indicator of the relation lattice**
`L(V) = { a ∈ Zᵏ : a·v = 0 }`. Equivalently `μ_V` is Haar measure on the closed
orbit (a circle, since integer speeds always give a *closed* geodesic). Then the
measure of lonely times is an exact sum over the relation lattice:

```
Leb(lonely times) = ∫_{Tᵏ} 1_B dμ_V = Σ_{a ∈ L(V)} ∏ᵢ f_n(aᵢ),
```

where `B = [1/n, 1-1/n]ᵏ` is the central box and
`f_n(0) = 1 - 2/n`, `f_n(a) = -sin(2πa/n)/(πa)` for `a ≠ 0`.
Verified numerically in `lonely_runner_lattice_sum_s363.py` (e.g. speeds
`(1,2,5)` at n=4: exact 0.2 vs lattice sum 0.1975; tight cases tend to 0
slowly).

Three things this representation makes obvious:

1. **The only thing the averaging "sees" about `V` is its relation lattice
   `L(V)`.** Two speed sets with the same `L` (e.g. `V` and `cV`) have identical
   lonely structure. The runners *are* their integer relations.

2. **`Leb(lonely) ≥ 0` is automatic**, so a pure density/Fourier statement can
   never be the whole conjecture. All the content sits in the **measure-zero
   tight stratum**, where the lattice sum is *exactly 0* but a boundary witness
   must still exist. This is a clean explanation of why Horvat–Stoffregen warn
   that classical Diophantine-approximation/density methods cannot settle LRC,
   and why the repo's endpoint-protection program (the measure-zero boundary)
   is the right altitude.

3. The box indicator's 1D coefficients decay only like `1/a` (a discontinuous
   function), which is why the tight stratum is delicate and why a finite-degree
   Fourier certificate (HYP-1812 Fejér/Riesz kernels) is fighting slow decay.
   The natural sharpening is the **Beurling–Selberg box minorant**: LRC holds
   for `L(V)` iff there is a trig polynomial `M ≤ 1_B` with
   `Σ_{a∈L(V)} hat{M}(a) > 0`. This is an explicit LP whose tight value at the
   threshold `1/n` is the conjecture. (Jensen's mixed-threshold Fourier
   expansions, `arXiv:2605.27941`, live in exactly this object.)

## 4. The torus-knot / view-obstruction representation

For *integer* speeds the orbit `t ↦ (v₁t,…,v_k t)` is **always a single closed
geodesic** of period 1 — the `(v₁,…,v_k)` torus knot on `Tᵏ`. The forbidden set
is the union of slabs `Sᵢ = {x : ‖xᵢ‖ < 1/n}` (the coordinate hyperplanes of
`Zᵏ`, thickened by `1/n`). LRC is then:

> every primitive `(v₁,…,v_k)` torus knot with distinct nonzero winding numbers
> pierces the central cube `B` of side `1 - 2/n` centered at the lattice.

This is Cusick's view-obstruction picture. For n = 14 the cube has side
`1 - 2/14 = 6/7`, and the forbidden shell is the slab of thickness `1/14` on each
side of every coordinate hyperplane (total slab width `1/7` around each
integer). The conjecture is the exact threshold where *every* admissible rational
direction's knot dips into the central cube.

**What a runner represents, in this picture:** runner `vᵢ` is the `i`-th
coordinate projection `Tᵏ → R/Z`, and its forbidden slab is the preimage of the
central bad arc under that projection. The speed set is a finite subset of the
*Pontryagin dual* `Ẑ = Z` of the circle — a set of *characters*. LRC is a
statement about subsets of the dual group: can a single point of the total space
avoid all the bad pullbacks at once.

Why the geometry says "expect no counterexample": the central box is a *small*
target — volume `(6/7)¹³ ≈ 0.135`, so the forbidden shell is the **bulk** (≈ 86%
of the torus). Naively that makes the box hard to hit. But the knot is a **long**
curve: it winds `vᵢ` times in coordinate `i`, so its length grows with the
speeds and it gets that many quasi-independent chances to fall into the 14% box.
Only a finely tuned resonance can keep the *entire* curve inside the 86% shell.
The extremal `(1,…,13)` knot is precisely the finest such resonance, and even it
only *touches* `∂B` (at the six unit points), never the interior — that is the
tight stratum again. A counterexample would be a knot that stays strictly in the
open shell for its whole length; the gap-floor data (§2) says forcing the
required `m=14` resonance into the winding vector instead *pushes* the knot into
the box.

## 5. The even/composite anatomy of 14 = 2·7

n = 14 is the first even composite frontier, and its arithmetic is worth
dissecting because the proof for even `n` probably needs tools the odd-prime
method lacks.

- **Antipodal (2-adic) tool.** Because n is even, `t = 1/2` makes
  `‖v/2‖ = 1/2 ≥ 1/14` for every *odd* `v` and `0` for every even `v`. So if all
  speeds are odd, `t=1/2` is instantly lonely. **A counterexample must contain
  an even speed** — automatically supplied by the required multiple of 14, but
  the involution `t ↦ t + 1/2` is a genuine structural symmetry available only
  for even `n`. This is the natural candidate to *replace* the missing odd-prime
  polynomial method: an even-`n` argument should exploit the `Z/2` antipodal
  involution / 2-adic descent.

- **CRT factoring of the unit skeleton.** `Z/14 ≅ Z/2 × Z/7`, and
  `(Z/14)* ≅ 1 × (Z/7)* ≅ Z/6` (cyclic, generated by 3:
  `3,9,13,11,5,1`). The six unit witnesses `a/14` are exactly the residues
  `(1 mod 2, a mod 7)` with `a ∈ (Z/7)*`. **The whole unit skeleton is a mod-7
  phenomenon** sitting on the odd residues mod 14; the 2-part is rigid. In this
  sense *n = 14 is "n = 7 doubled."* Since the 7-runner case is classically
  proven (Barajas–Serra), a promising — if speculative — descent is to factor a
  hypothetical n=14 counterexample through its mod-7 and mod-2 shadows and
  derive a contradiction with the proven n = 7 case. This is the sieve-descent
  HYP-1813 read through the prime factorization of 14.

## 6. Proof vs disproof: honest assessment

**Disproof (counterexample at n = 14):** essentially not expected. A
counterexample is now boxed in by: (a) it must realize divisibility by every
`m ∈ {2,…,14}` (so it contains speeds divisible by 11 and by 13 — "expensive"
large-modulus speeds); (b) total forbidden measure budget is only
`13·(2/14) = 13/7 ≈ 1.857`, so overlaps must absorb `≈0.857` with full coverage
and zero gap; (c) every endpoint must be protected (S359); (d) the gap-floor
search found *no* tight primitive set with a multiple of 14. The knife-edge
structure — the extremal example is an isolated tight point and every
perturbation opens a gap — is strong heuristic evidence against any
counterexample.

**Proof:** the realistic route mirrors the recent papers — bound the minimal
counterexample (Tao), sieve residues (the modulus-sieve lemma, massively
restricting allowed residue patterns), then dispatch the residual
near-initial-segment families. The **open gap specific to n = 14** is that the
odd-prime polynomial method does not apply (14 is even/composite). The new
suggestion of this session: for even `n`, use the **`t ↦ t+1/2` antipodal
involution / 2-adic structure** as the substitute, and factor the residual
analysis through the proven `n = 7` case via the `14 = 2·7` CRT decomposition.
A complete proof is not in reach from this workspace today, but the *shape* of
what is missing is now sharp.

## 7. Concrete next steps

1. **Attack the gap-floor hypothesis** structurally: show any primitive 13-set
   with a multiple of 14 has a positive gap. Even a proof restricted to sets
   whose other 12 speeds lie in `{1,…,K}` for small `K` would be new.
2. **Beurling–Selberg LP at n = 14**: build the box minorant `M` of fixed degree
   and minimize `Σ_{a∈L(V)} hat{M}(a)` over the sieve-admissible lattices; a
   positive certificate for all admissible `L` would be an analytic proof.
3. **CRT-descent experiment**: take candidate near-counterexamples, reduce mod 7
   and mod 2, and test whether the proven n = 7 structure forbids the
   reduction — make HYP-1813 concrete at `14 = 2·7`.
4. **Classify tight 13-sets**: is the initial segment the *only* primitive tight
   set at n = 14, or are there sporadic ones (as at n = 8)? A targeted search in
   the sieve-admissible neighborhood would settle the local picture.

## Appendix: artifacts produced this session

```
04-computation/lonely_runner_n14_frontier_s363.py   (modulus-sieve / whack-a-mole)
04-computation/lonely_runner_n14_gapfloor_s363.py   (gap-floor stress test)
04-computation/lonely_runner_lattice_sum_s363.py    (pushforward identity check)
05-knowledge/results/lonely_runner_n14_frontier_s363.out
05-knowledge/results/lonely_runner_n14_gapfloor_s363.out
05-knowledge/results/lonely_runner_lattice_sum_s363.out
```

## Parallel work (added at merge time)

On pushing this session I found that **codex independently ran sessions
S363–S382 on the same n = 14 frontier** while I worked, and converged on the
*same* central tension: codex **HYP-1839** ("fourteen-runner tightness and
counterexamples require opposite gate behavior") is the same statement as the
gap-floor here (§2), and codex **HYP-1816** (CRT-gate descent) and **HYP-1832**
(torsion/CRT leak dichotomy) overlap the modulus-sieve and `14 = 2·7` CRT lenses.
This is independent confirmation that the *opposite-gate / no-multiple-of-n*
tension is the right way to see n = 14. The distinct contributions of this
session are: (i) the **general** "no-multiple-of-`n` at full measure" form with
its n = 4..8 corroboration; (ii) the **pushforward lattice-sum identity** locating
all difficulty in the measure-zero stratum; (iii) the explicit **Beurling–Selberg
LP** sharpening. My hypotheses are indexed as HYP-1850 / 1841 / 1842, cross-linked
to the codex entries.
