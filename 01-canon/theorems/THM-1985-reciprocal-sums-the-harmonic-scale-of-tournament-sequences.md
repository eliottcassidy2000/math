---
id: THM-1985
title: "RECIPROCAL SUMS AS A HARMONIC-SCALE FACE OF THE POLY/#P TOWER (CORRECTED SUPPORT SEMANTICS)"
status: "PARTIALLY RETRACTED by MISTAKE-209 and repaired by THM-2000. The exact simplex/figurate telescoping identities remain PROVED. Census rows are finite numerical term-multiset prefixes, not proved transcendental constants; global H-spectrum divergence is OPEN; support masses require collision-tax correction."
author: opus-2026-07-20-S447
depends_on: [THM-819, THM-1920, THM-1930, THM-1970, THM-1975, THM-1980, THM-805, THM-1926]
related: [THM-2000, MISTAKE-209, THM-1370-h-spectrum-omits-7-21-all-n.md]
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
> The row formerly called `c_3`-max is also corrected: `C(n,3)` counts all
> triple slots, not cyclic triples; the true maximum-cyclic sequence has mass
> `75/4-24log2` (THM-2000).
> Gauss's triangular-number identity also gives an exact product/theta form for
> the labeled-tournament support mass.

Owner: the reciprocal of an integer sequence is a **sub-series of the harmonic series**; study
`Σ 1/a_n` for the repo's sequences (figurate reciprocals, Abel–Dini, Bertrand). This gives the
**reciprocal-sum face** of THM-1970/1975's formula/`#P` harmonic edge: a sequence's *growth* (its
place in the poly tower) is exactly its reciprocal sum's *convergence*.

## 1. Figurate invariant-sizes → rational sums `k/(k−1)`

`char_S`'s degree-`k` coefficient is `C(n,k)`-shaped (THM-1920: the subleading `= C(n,2) = #arcs`),
and

> **`Σ_{n≥k} 1/C(n,k) = k/(k−1)`** (telescoping / hockey-stick), verified `k=2..6`.

So the **tournament's own sizes** have rational reciprocal sums:

| size = char_S coefficient | figurate | `Σ 1/size` |
|---|---|---|
| arc count `C(n,2)` (subleading) = triangular `T_{n−1}` | triangular | **2** |
| # tiles `C(n−1,2)` | triangular | 2 |
| all triple slots `C(n,3)` (**not** `c₃`-max) | tetrahedral | **3/2** |
| var-max `2·C(n,3)` (transitive, THM-1930) | 2·tetrahedral | **3/4** |

**`Σ 1/(arc count) = 2` is the Downey–Ong–Sellers triangular identity realized on the tournament** —
the `char_S` subleading-coefficient series sums to exactly `2`, while the plain harmonic partial sum
already passes `2` by `n=5` (the owner's contrast).  The **pure binomial
column** `C(n,k)` has reciprocal sum `k/(k−1)`; degree alone does not force
this value for an arbitrary polynomial invariant.

## 2. Counting sequences → fast-converging numerical constants

Super-exponential growth makes the reciprocal sum converge fast; it does **not**
by itself imply irrationality or transcendence.  The original computation
printed termwise prefixes for A000568, A038375, A051337, A002854, and A000255,
but its offsets were not canonical and its A000568 list contained a wrong
term.  Those decimals are historical diagnostics, not sequence constants;
THM-2000 retains only explicitly ranged, support-deduplicated prefixes.
The **Cayley–Dickson levels** `n=2^k+1` give `Σ_{k≥1} 1/(2^k+1) = 0.7645` (Erdős–Borwein cousin of
`Σ1/(2^k−1)=1.6067`); the `H=1+2^{n−2}` SC-neighbour series sums to `1.2645`.

## 3. The H-value spectrum → conditional divergence target

`THM-1370-h-spectrum-omits-7-21-all-n.md` proves that `7,21` never occur and that every other odd value through
`609` occurs; it states global completeness as a **conjecture**.  Therefore
linear density and divergence of `Σ1/(H-value)` remain open.  Conditional on
positive lower density, the intended harmonic-edge conclusion follows.  The
independent classical Abel--Dini statement remains valid for `a_n>0`:
for a divergent `Σa_n` with partial sums `S_n`, `Σ a_n/S_n` diverges but `Σ a_n/S_n^{1+ε}` converges
for every `ε>0` — **there is no convergence at the exact exponent-one boundary**, the precise analogue of **kps
THM-1980** ("Rédei parity is the last formula"; no poly invariant beats the last bit). The
**Bertrand** boundary is `Σ 1/(n \ln n)` (diverges), `Σ 1/(n(\ln n)^α)` converges iff `α>1`.

## The picture

```
  reciprocal sum                sequence                     tower position
  ────────────────────────────────────────────────────────────────────────
  rational k/(k-1)     figurate invariant-SIZES (char_S)    poly, degree k   (deep convergence)
  numerical constant   COUNTING seqs (A000568/38375/…)      the object census
  OPEN                  the global H-VALUE support           conjectural edge
```

**A sequence's reciprocal sum measures its position on the harmonic scale, and this recovers the
poly/`#P` tower**: the polynomial (figurate) invariant-sizes converge to rationals (formula-land),
the `H`-value set would diverge under the still-open positive-density conjecture, while the counting sequences fill the middle. The Downey–Ong–Sellers
triangular `2` is the tournament's own size; the harmonic edge of THM-1970 is where `H`'s values
live; Abel–Dini/Bertrand are the microscope at the boundary.

## 4. The harmonic constants are already in the repo — and meet the figurate sizes at the LRC extremal

A full sequence sweep (S447 agent) confirms the repo already contains the harmonic number, `γ`, and
a ζ-analogue, and — the key synthesis — **they meet the figurate sizes of §1 exactly at the LRC
flagship extremal**:

- **`THM-819` (correcting the old all-`k` reading of THM-805):** the general
  base measure is the **primitive harmonic sum over a triangular denominator**,
  `2/[(k+1)(k+2)] sum_{u<=k,(u,k+1)=1}1/u`.  It equals
  `H_k/C(k+2,2)` exactly when `k+1` is prime.  Thus the flagship `k=12`
  instance really is
  `H₁₂/C(14,2) = 6617/194040` (exact, `= mac-mini THM-736`'s `|G′({1..12})|`),
  because `13` is prime; the formula is not an all-`k` law. The staircase
  Smith-network resistance is `H_{n−2}` (`THM-805`/HYP-6865). **So the divergent-edge object (the
  harmonic number `H_k`) and the convergent figurate size (the triangular `C(k+2,2)`, `Σ1/·=2` of
  §1) combine into the single most important LRC extremal measure.**
- **The loneliness floor spectrum IS the harmonic series:** the `n`-runner floor is `M=1/n`
  (`CONSTANTS-INDEX`), so `Σ_n M_floor(n) = Σ 1/n` diverges — the LRC floor values *are* the harmonic
  terms. (And the deep-well `M = n/Φ₆(n)`, `Φ₆(n)=n²−n+1` cyclotomic; `Σ1/Φ₆(n)=1.798` converges.)
- **`THM-1926` — the tournament ζ:** `ζ_T(u)=1/det(I−uA)=∏_{\text{primitive cycles }p}(1−u^{ℓ(p)})^{−1}`
  — the Euler product with **cycles as primes**, `≡1` on the acyclic part. The determinant-side
  (poly) generating function of the moment tower.
- **`THM-1970`/reflection** already fixes the `γ`-analogue: the `char_S→H` **defect** is the
  tournament `γ` (the finite part after the poly tower is subtracted).
- **Refinement (agent correction):** only the *pure binomials* `C(n,k)` give rational sums `k/(k−1)`;
  the `1+2^d`, `2^{n−2}+1`, and Cayley–Dickson `2^k+1` families converge to
  **Erdős–Borwein-type numerical constants** (§2); no arithmetic type is
  asserted here.  Staircase 3-cycles
  `k(k−1)` telescope to **1**; `A002088` and `A001764` are **absent** from the repo.

So the reciprocal-sum program closes a loop: the figurate `2` (triangular = tournament size), the
harmonic `H_k` (resistance / the divergent edge), and the `ζ_T` Euler product are three faces of one
structure, and the LRC deep-well — the extremal the whole covering program orbits — is literally
`H_k / (\text{triangular})`, harmonic ÷ figurate.

## Open

1. **Certify the census tails before naming constants.**  The displayed
   decimals are finite historical prefixes (and the S447 A000568 list contains
   a wrong term), so continued-fraction or inverse-symbolic searches are not
   meaningful until canonical offsets, support deduplication, and rigorous
   tail intervals are supplied.
2. **The global H-value support.** First prove enough density to decide
   divergence.  Only then ask whether `Σ 1/(H-value)` has an asymptotic
   `c log x` and compute `c`.
3. **Bertrand placement of the near-boundary sequences.** Does any repo sequence grow like
   `n(\ln n)^α` (the Bertrand scale) rather than polynomially or exponentially? (Candidate: the
   LRC covering-modulus / minimal-witness-denominator sequences.)

## Verification

The S447 script/output are historical and retracted as a canonical referee by
MISTAKE-209/213.  Use
`04-computation/support_harmonic_abel_dini_figurate_surface_thm2000.py` and
THM-2000/2005 for the optimization-safe support, Abel--Dini, figurate,
automatic, and tournament certificates.
