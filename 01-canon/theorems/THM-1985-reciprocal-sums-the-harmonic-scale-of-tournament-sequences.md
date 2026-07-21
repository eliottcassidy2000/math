---
id: THM-1985
title: "RECIPROCAL SUMS AS A HARMONIC-SCALE FACE OF THE POLY/#P TOWER (CORRECTED SUPPORT SEMANTICS)"
status: "PARTIALLY RETRACTED by MISTAKE-209 and repaired by THM-2000/2005. The exact simplex/figurate telescoping identities remain PROVED. Census rows are finite numerical term-multiset prefixes, not proved transcendental constants; global H-spectrum divergence is OPEN; support masses require collision-tax correction."
author: opus-2026-07-20-S447
depends_on: [THM-1920, THM-1930, THM-1970, THM-1975, THM-1980, THM-805, THM-819, THM-1926]
related: [THM-2000, THM-2005, MISTAKE-209, THM-1370-h-spectrum-omits-7-21-all-n.md]
external: "Downey-Ong-Sellers (CMJ), sums of reciprocals of figurate numbers (triangular=2); Abel-Dini theorem; Bertrand series; Erdos-Borwein constant."
cite_by_filename: true
---

# THM-1985 — Reciprocal sums: the harmonic-scale face of the poly/#P tower

> **CORRECTION (MISTAKE-209; THM-2000).**  The exact figurate identities in
> Section 1 survive.  Sections 2--3 originally counted repeated sequence terms
> despite saying “subset,” called unnamed numerical constants transcendental
> without proof, and treated `THM-1370-h-spectrum-omits-7-21-all-n.md`'s conjectural global H-spectrum
> completeness as proved.  Support masses remove repetitions; arithmetic type
> is open unless proved; and H-spectrum reciprocal divergence is conditional.
> Gauss's triangular-number identity also gives an exact product/theta form for
> the labeled-tournament support mass.  THM-2005 additionally corrects the
> cyclic-triangle maximum and gives the full support Dirichlet profile.

Owner: the reciprocal of an integer sequence is a **sub-series of the harmonic series**; study
`Σ 1/a_n` for the repo's sequences (figurate reciprocals, Abel–Dini, Bertrand). This gives the
**reciprocal-sum analogy** for THM-1970/1975's formula/`#P` edge.  It is not
an equivalence with computational complexity.

## 1. Figurate invariant-sizes → rational sums `k/(k−1)`

`char_S`'s degree-`k` coefficient is `C(n,k)`-shaped (THM-1920: the subleading `= C(n,2) = #arcs`),
and

> **`Σ_{n≥k} 1/C(n,k) = k/(k−1)`** (telescoping / hockey-stick), verified `k=2..6`.

So the **tournament's own sizes** have rational reciprocal sums:

| size = char_S coefficient | figurate | `Σ 1/size` |
|---|---|---|
| arc count `C(n,2)` (subleading) = triangular `T_{n−1}` | triangular | **2** |
| # tiles `C(n−1,2)` | triangular | 2 |
| triple-slot count `C(n,3)` (**not** `c3`-max) | tetrahedral | **3/2** |
| var-max `2·C(n,3)` (transitive, THM-1930) | 2·tetrahedral | **3/4** |

**`Σ 1/(arc count) = 2` is the Downey–Ong–Sellers triangular identity realized on the tournament** —
the `char_S` subleading-coefficient series sums to exactly `2`, while the plain harmonic partial sum
already passes `2` by `n=5` (the owner's contrast).  The formula
`k/(k-1)` applies to these pure binomial rows `C(n,k)`; degree alone does not
determine a general invariant's reciprocal mass.

The true maximum-cyclic-triangle row is
`sum_(n>=3)1/M_3(n)=75/4-24 log 2` by THM-2005.

## 2. Counting sequences → fast-converging numerical constants

Super-exponential growth makes the reciprocal sum converge fast; it does **not**
by itself imply irrationality or transcendence.  The original computation
printed **indexed** prefix values, not support masses or proved closed forms:
A000568 was `2.8535...` when started at `n=1` (and is `3.8535...` when
started at `n=0`), while the A038375, A051337, A002854, and A000255 rows were
likewise offset-sensitive finite prefixes.  THM-2005 instead deduplicates
A000568 and proves the support-mass interval
`1.853534132290116317333715265823800054971816577029176...`
`< D(1) <`
`1.853534132290116317333715265823800054971816608078302...`.
The **Cayley–Dickson levels** `n=2^k+1` give `Σ_{k≥1} 1/(2^k+1) = 0.7645` (Erdős–Borwein cousin of
`Σ1/(2^k−1)=1.6067`).  The realizable neighbor support
`{1+2^k:k>=1}={3,5,9,...}` has mass `0.764499...`; the old
`1.2645` termwise row included the impossible even value `H=2`.

## 3. The H-value spectrum → conditional divergence target

`THM-1370-h-spectrum-omits-7-21-all-n.md` proves that `7,21` never occur and that every other odd value through
`609` occurs; it states global completeness as a **conjecture**.  Therefore
linear density and divergence of `Σ1/(H-value)` remain open.  Conditional on
positive lower density, the intended harmonic-edge conclusion follows.  The
independent classical Abel--Dini statement remains valid:
for a divergent `Σa_n` with partial sums `S_n`, `Σ a_n/S_n` diverges but `Σ a_n/S_n^{1+ε}` converges
for every `ε>0`.  Exponent one is the divergent critical member.  The comparison
with **kps THM-1980** ("Rédei parity is the last formula") is heuristic, not an
identification of analytic and complexity parameters.  The
**Bertrand** boundary is `Σ 1/(n \ln n)` (diverges), `Σ 1/(n(\ln n)^α)` converges iff `α>1`.

## The picture

```
  reciprocal sum                sequence                     tower position
  ────────────────────────────────────────────────────────────────────────
  rational k/(k-1)     figurate invariant-SIZES (char_S)    poly, degree k   (deep convergence)
  numerical constant   COUNTING seqs (A000568/38375/…)      the object census
  OPEN                  the global H-VALUE support           conjectural edge
```

This is a heuristic harmonic-scale comparison, not a recovery theorem for the
poly/`#P` tower.  The full support profile and Abel criterion in THM-2000/2005
retain the tail data that one reciprocal scalar erases.

## 4. The harmonic constants are already in the repo — and meet the figurate sizes at the LRC extremal

A full sequence sweep (S447 agent) found the harmonic number and a genuine
tournament zeta function, plus a proposed `gamma` analogy.  The exact LRC
connection is the primitive harmonic measure below:

- **`THM-819` (correcting the old THM-805 generalization):** the LRC
  interval-core measure is the **primitive** harmonic sum
  `2/[q(q+1)] sum_(u<q,(u,q)=1)1/u`, with `q=k+1`.  It equals
  `H_k/C(k+2,2)` exactly when `q` is prime.  At `k=12,q=13`:
  `H₁₂/C(14,2) = 6617/194040` (exact, `= mac-mini THM-736`'s `|G′({1..12})|`). The staircase
  Smith-network resistance is `H_{n−2}` (`THM-805`/HYP-6865). **So the divergent-edge object (the
  harmonic number `H_k`) and the convergent figurate size (the triangular `C(k+2,2)`, `Σ1/·=2` of
  §1) combine into the single most important LRC extremal measure.**
- **The loneliness floor spectrum IS the harmonic series:** the `n`-runner floor is `M=1/n`
  (`CONSTANTS-INDEX`), so `Σ_n M_floor(n) = Σ 1/n` diverges — the LRC floor values *are* the harmonic
  terms. (And the deep-well `M = n/Φ₆(n)`, `Φ₆(n)=n²−n+1` cyclotomic; `Σ1/Φ₆(n)=1.798` converges.)
- **`THM-1926` — the tournament ζ:** `ζ_T(u)=1/det(I−uA)=∏_{\text{primitive cycles }p}(1−u^{ℓ(p)})^{−1}`
  — the Euler product with **cycles as primes**, `≡1` on the acyclic part. The determinant-side
  (poly) generating function of the moment tower.
- **Proposed analogy only (`THM-1970` reflection):** the `char_S -> H` defect
  was compared with a tournament `gamma`; no canonical limiting constant is
  proved.
- **Refinement (agent correction):** only the *pure binomials* `C(n,k)` have
  the proved rational sums `k/(k−1)` here.  The `1+2^d`, `2^{n−2}+1`, and
  Cayley–Dickson `2^k+1` families converge rapidly, but no arithmetic type is
  inferred from that growth (§2).  Staircase 3-cycles
  `k(k−1)` telescope to **1**; `A002088` and `A001764` are **absent** from the repo.

So the reciprocal-sum program closes a loop: the figurate `2` (triangular = tournament size), the
harmonic `H_k` (resistance / the divergent edge), and the `ζ_T` Euler product are three faces of one
structure, and the LRC deep-well — the extremal the whole covering program orbits — is literally
`H_k / (\text{triangular})`, harmonic ÷ figurate.

## Open

1. **Normalize and certify the other counting supports.** Fix canonical
   offsets, remove collisions, and prove tails before attempting arithmetic
   recognition.  THM-2005 completes this for A000568; A038375 and the other
   census rows remain open.
2. **The H-value support.** Prove or refute global completeness outside
   `7,21`; only then ask for a density or logarithmic reciprocal coefficient.
3. **Bertrand placement of the near-boundary sequences.** Does any repo sequence grow like
   `n(\ln n)^α` (the Bertrand scale) rather than polynomially or exponentially? (Candidate: the
   LRC covering-modulus / minimal-witness-denominator sequences.)

## Verification

`04-computation/reciprocal_sums_of_repo_sequences_opus_S447.py` (+ `.out`)
is retained as a superseded exploratory prefix calculator.  Its simplex rows
survive; use the THM-2000/2005 referees for support-normalized claims.
