---
id: THM-1840
title: "Two-charge constant-term nonvanishing and the charge-projecting functional corollary"
status: >
  PROVED, WITH THE FUNCTIONAL SCOPE CORRECTED BY THM-2070. The charge-0 coefficient of
  (a u^p + b u^{−q})^m is a single term at the return times (the only (k,m−k) with pk = q(m−k)):
  elementary and machine-verified across (p,q) = (1,1),(1,2),(2,3),(3,5),(2,2). The
  functional corollary holds for charge-PROJECTING F, meaning F(u^k)=0 for k!=0 and F(1)!=0.
  The former weaker phrase "depends only on charge" did not imply this selection rule and is
  retracted; THM-2070 gives an explicit counterexample. The general many-charge case is not
  closed by the binomial formula: its support returns can cancel coefficientwise.
source: klein-2026-07-20-S393 (owner: prove the abstract 'both-signs single-character nullcone non-vanishing' by the functional-agnostic method; think cyclotomic)
depends_on: []  # the binomial/charge-balance proof below is self-contained
related:
  - THM-1830-dvdez-conditional-edifice-toral-route-lrc
  - THM-1745-multistraddle-combines-by-max-and-the-moment-bound-is-a-return-time
  - THM-1770-gmc2-detection-depth-grows-with-radial-degree
  - THM-1810-the-bosonic-permanent-side-is-why-gmc2-is-hard
  - THM-1805-the-vandermonde-is-a-signed-tournament-sum-intransitivity-cancels
  - THM-790-blue-parity-law-proved
  - THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation
---

# THM-1840 — two-charge nonvanishing and charge projection

## The theorem

Let `Λ(u) = a·u^p + b·u^{−q}` with `a,b ≠ 0` and `p,q > 0` — **both-signs** (one positive
character `+p`, one negative `−q`), **single-character** on each sign. Put `g = gcd(p,q)`,
`m₀ = (p+q)/g`.

> **The charge-0 coefficient of `Λ^m` is a single term:**
> `[u^0]Λ^m = C(m, qm/(p+q))·a^{qm/(p+q)}·b^{pm/(p+q)}` when `m₀ | m`, and `0` otherwise.
> It is **uncancellable** (the unique `(k,m−k)` with `pk = q(m−k)`) and **nonzero** at every
> return time `m ∈ m₀·ℤ_{>0}`.

**Charge-projecting corollary.** Let `F` be a linear functional on
`ℂ[u,u^{−1}]` satisfying

```text
F(u^k)=0 for k!=0,                 F(1)!=0.
```

Then `F(h)=F(1) CT(h)` for every Laurent polynomial `h`, and therefore

```text
F(Λ^m) = F(u^0)·[u^0]Λ^m ≠ 0   for every m ∈ m₀·ℤ_{>0}.
```

The functional supplies the scalar `F(1)`; **the charge lattice supplies the non-vanishing.**
So the both-signs single-character nullcone is non-vanishing for every
charge-projecting functional at once — verified on three genuinely different ones:

| `F` | role | `F(Λ^5)`, `Λ = 2u²+3u^{−3}` (`m₀=5`) |
|---|---|---|
| `CT` (constant term) | **TNC / toral** | `720` |
| radial-weighted (`w₀=7`) | **GMC / Gaussian** | `5040` |
| Bonferroni-weighted (`w₀=5`) | **LRC** | `3600` |

Same `[u^0]Λ^5 = 720`, three charge-projecting functionals, all nonzero.

The selection rule is load-bearing. Merely assigning an arbitrary value to
each nonzero charge does not yield the displayed identity: THM-2070 constructs
a diagonal-by-charge functional for which the nonzero-charge terms cancel the
constant term. That earlier weaker definition is withdrawn.

## Cyclotomic content

The charge-0 condition `p·k = q·(m−k)` is **balance on the charge lattice** — equivalently, the
constant-term projection is the average over the `(p+q)`-th roots of unity, and the return times
are the `m` for which the positive and negative characters **realign** modulo `p+q`. The first
return `m₀ = (p+q)/gcd(p,q)` is exactly mac-mini's coprime-pair realignment time (THM-1745) and
the Lonely-Runner two-runner realignment. The both-signs non-vanishing is thus a *cyclotomic*
fact, which is why it transfers across charge-projecting functionals: roots of
unity do not care which scalar such a functional assigns to the constant term.

## The residual is many-charge cancellation — and it is THM-1770

With a **middle** charge present, the extreme-pair charge-0 term is no longer alone:

```text
Λ = u^1 + B·u^0 + u^{−1},   [u^0]Λ² = 2·(1·1) + B²  =  2 + B².
```

The extreme-pair term `2` is **cancellable** by the charge-0 self-term `B²`, and it cancels at
`B = i√2` — which is *exactly* the THM-1770 depth-growth witness
`P = Z + (−i√2 + i√2·|Z|²) + Z̄`. So:

> **Single character (binomial): uncancellable, closed under charge projection at the return
> time. One extra (middle) character: the first cancellation, and it is the witness that makes
> the detection depth exceed the span (THM-1770).**

The many-charge stratum is where several cyclotomic returns collide and can
cancel — the `bosonic/permanent, no sign to cancel` wall of THM-1810. It is
resolved globally by THM-2022, but not by this two-charge lemma.

## The atoms, and the parity lore (owner's threads)

- **The 3-cycle atom.** A directed 3-cycle is the minimal *intransitive* unit (THM-1805, the
  cancelling unit of the Vandermonde) — the tournament face of the both-signs *atom*: a minimal
  balanced (charge-0-reachable) configuration. One atom = a single-character straddle, closed
  here; **two independent atoms** (two coprime-pair straddles with distinct returns) is where the
  collisions of §"residual" begin — the boundary between the charge-projecting zone and the
  GMC-hard zone. (The owner's "`n ≥ 13` admits two 3-cycle atoms" is the tournament/LRC-14
  incarnation of this two-atom threshold; the charge-side statement is that two straddles first
  produce a colliding charge-0 term there. Flagged for its own pass.)
- **Blue parity (THM-790: `1` at odd `n`, `0` at even `n`).** The blue (grid-symmetric =
  `σ`-fixed) count parity is the **fixed-point parity of an involution**, the same shape as the
  charge-0 projection's parity: the constant term is the average over the charge involution
  `q ↦ −q`, whose fixed points (charge 0) carry the surviving weight, exactly as `σ`-fixed
  (blue) tilings carry the SC weight (THM-1415 §1). The `1/0` odd/even alternation is that
  fixed-point count mod 2 — a clean instance of the project's SC/complement-parity lore living on
  the same charge-involution as the non-vanishing here.

## Scope

The single-character (two-charge) both-signs nullcone is now **closed for every
charge-projecting functional simultaneously**, cyclotomically — the rigorous form of THM-1830's common parent. The
many-charge case is the first cancellation stratum, pinpointed by THM-1770 and
framed by THM-1810. This improves the mechanism ledger: TNC, GMC, and LRC
share a proved two-charge base case, while beyond it one must retain the
functional's actual orbit weights. THM-2022 supplies the required whole-face
argument for GMC(2).

*Files: `04-computation/functional_agnostic_nullcone_klein_S393.py` (+ `.out`).*
