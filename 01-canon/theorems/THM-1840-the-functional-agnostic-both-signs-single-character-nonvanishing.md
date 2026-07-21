---
id: THM-1840
title: "THE FUNCTIONAL-AGNOSTIC BOTH-SIGNS SINGLE-CHARACTER NULLCONE NON-VANISHING — PROVED, and it is the rigorous common parent of THM-1830. For a two-charge Laurent polynomial Λ(u) = a·u^p + b·u^{−q} (a,b ≠ 0, p,q > 0 = both-signs, one positive and one negative character), the charge-0 coefficient of Λ^m is a SINGLE, uncancellable term [u^0]Λ^m = C(m, qm/(p+q))·a^{qm/(p+q)}·b^{pm/(p+q)}, nonzero exactly when the coprime-pair first-return time m₀ = (p+q)/gcd(p,q) divides m. THEREFORE, for ANY charge-graded functional F (F depends on a monomial only through its charge, F(u^0) ≠ 0), F(Λ^m) = F(u^0)·[u^0]Λ^m ≠ 0 at every return time — proving the both-signs single-character nullcone non-vanishing SIMULTANEOUSLY for TNC (F = CT), GMC (F = radial/Gaussian), and LRC (F = Bonferroni), because the non-vanishing is a fact of the CHARGE LATTICE and the functional only supplies a nonzero weight. The mechanism is CYCLOTOMIC: p·k = q·(m−k) is balance on the charge lattice = the (p+q)-th roots of unity, and m₀ = (p+q)/gcd is mac-mini's coprime-pair realignment (THM-1745). The residual is exactly many-charge cancellation: once a middle charge exists, the extreme-pair charge-0 term becomes cancellable — and the {−1,0,1} cancellation [u^0]Λ² = 2 + B² = 0 at B = i√2 IS the THM-1770 depth-growth witness."
status: >
  PROVED (the single-character = binomial case). The charge-0 coefficient of
  (a u^p + b u^{−q})^m is a single term at the return times (the only (k,m−k) with pk = q(m−k)):
  elementary and machine-verified across (p,q) = (1,1),(1,2),(2,3),(3,5),(2,2). The
  functional-agnostic corollary is one line (F charge-graded ⟹ F(Λ^m) = F(u^0)·[u^0]Λ^m), and is
  demonstrated on three distinct functionals (CT, a radial-weighted, a Bonferroni-weighted) giving
  the same non-vanishing up to the scalar F(u^0). The general many-charge case is NOT closed by
  this — it is the cancellation residual, and the {−1,0,1} witness is exhibited as its first
  instance (= THM-1770). Proves the two-charge stratum uniformly across functionals; the
  many-charge nullcone (GMC-hard) remains open.
source: klein-2026-07-20-S393 (owner: prove the abstract 'both-signs single-character nullcone non-vanishing' by the functional-agnostic method; think cyclotomic)
depends_on:
  - THM-1830  # the conditional map naming the common return-time parent — this proves its clean case
  - THM-1745  # mac-mini: the coprime-pair first-return time (p+q)/gcd = m₀
related:
  - THM-1770  # the {−1,0,1} depth-growth witness = the first cancellation past the single-character case
  - THM-1810  # bosonic/fermionic: the many-charge cancellation is the 'no sign to cancel' wall
  - THM-1805  # the 3-cycle as the intransitivity atom (the both-signs atom, tournament side)
  - THM-790   # the blue parity law (1 odd / 0 even) — the charge-0/fixed-point parity lore
---

# THM-1840 — the functional-agnostic both-signs single-character non-vanishing

## The theorem

Let `Λ(u) = a·u^p + b·u^{−q}` with `a,b ≠ 0` and `p,q > 0` — **both-signs** (one positive
character `+p`, one negative `−q`), **single-character** on each sign. Put `g = gcd(p,q)`,
`m₀ = (p+q)/g`.

> **The charge-0 coefficient of `Λ^m` is a single term:**
> `[u^0]Λ^m = C(m, qm/(p+q))·a^{qm/(p+q)}·b^{pm/(p+q)}` when `m₀ | m`, and `0` otherwise.
> It is **uncancellable** (the unique `(k,m−k)` with `pk = q(m−k)`) and **nonzero** at every
> return time `m ∈ m₀·ℤ_{>0}`.

**Functional-agnostic corollary.** A functional `F` on `ℂ[u,u^{−1}]` is *charge-graded* if
`F(monomial)` depends only on the monomial's charge and `F(u^0) ≠ 0`. Then

```text
F(Λ^m) = F(u^0)·[u^0]Λ^m ≠ 0   for every m ∈ m₀·ℤ_{>0}.
```

The functional supplies the scalar `F(u^0)`; **the charge lattice supplies the non-vanishing.**
So the both-signs single-character nullcone is non-vanishing for *every* charge-graded functional
at once — verified on three genuinely different ones:

| `F` | role | `F(Λ^5)`, `Λ = 2u²+3u^{−3}` (`m₀=5`) |
|---|---|---|
| `CT` (constant term) | **TNC / toral** | `720` |
| radial-weighted (`w₀=7`) | **GMC / Gaussian** | `5040` |
| Bonferroni-weighted (`w₀=5`) | **LRC** | `3600` |

Same `[u^0]Λ^5 = 720`, three functionals, all nonzero. This is the *method*: prove the nullcone
statement once, on the charge lattice, and read it off for TNC, GMC, and LRC together.

## Cyclotomic content

The charge-0 condition `p·k = q·(m−k)` is **balance on the charge lattice** — equivalently, the
constant-term projection is the average over the `(p+q)`-th roots of unity, and the return times
are the `m` for which the positive and negative characters **realign** modulo `p+q`. The first
return `m₀ = (p+q)/gcd(p,q)` is exactly mac-mini's coprime-pair realignment time (THM-1745) and
the Lonely-Runner two-runner realignment. The both-signs non-vanishing is thus a *cyclotomic*
fact, which is why it is functional-agnostic: roots of unity do not care which functional weights
the constant term.

## The residual is many-charge cancellation — and it is THM-1770

With a **middle** charge present, the extreme-pair charge-0 term is no longer alone:

```text
Λ = u^1 + B·u^0 + u^{−1},   [u^0]Λ² = 2·(1·1) + B²  =  2 + B².
```

The extreme-pair term `2` is **cancellable** by the charge-0 self-term `B²`, and it cancels at
`B = i√2` — which is *exactly* the THM-1770 depth-growth witness
`P = Z + (−i√2 + i√2·|Z|²) + Z̄`. So:

> **Single character (binomial): uncancellable, closed functional-agnostically at the return
> time. One extra (middle) character: the first cancellation, and it is the witness that makes
> the detection depth exceed the span (THM-1770).**

The many-charge nullcone is where several cyclotomic returns collide and can cancel — the
`bosonic/permanent, no sign to cancel` wall of THM-1810, and the open heart of GMC(2).

## The atoms, and the parity lore (owner's threads)

- **The 3-cycle atom.** A directed 3-cycle is the minimal *intransitive* unit (THM-1805, the
  cancelling unit of the Vandermonde) — the tournament face of the both-signs *atom*: a minimal
  balanced (charge-0-reachable) configuration. One atom = a single-character straddle, closed
  here; **two independent atoms** (two coprime-pair straddles with distinct returns) is where the
  collisions of §"residual" begin — the boundary between the functional-agnostic zone and the
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

The single-character (two-charge) both-signs nullcone is now **closed for every charge-graded
functional simultaneously**, cyclotomically — the rigorous form of THM-1830's common parent. The
many-charge case is the residual, pinpointed as the cancellation THM-1770 exhibits and THM-1810
frames. This improves the routing: TNC, GMC, and LRC share a *proved* two-charge base case, and
differ only past it, where the functional's radial weight decides whether the colliding returns
survive.

*Files: `04-computation/functional_agnostic_nullcone_klein_S393.py` (+ `.out`).*
