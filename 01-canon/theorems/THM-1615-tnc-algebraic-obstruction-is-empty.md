---
id: THM-1615
title: "THE TNC's ALGEBRAIC OBSTRUCTION IS EMPTY: A GENUINE SADDLE ALWAYS EXISTS, AT EVERY BIDEGREE AT ONCE. (0) REFUTED FIRST, my own: the necessary condition 'a TNC violator must have disc(R) = 0' is FALSE — R = (1+u)^2 at N=1 has disc 0 AND nontrivial gcd(R, uR'-NR) = u+1, yet CT(Lambda^m) = C(2m,m) = 2,6,20,70,... never vanishes; likewise (1+u)^4 at N=2 gives C(4m,2m). The error: a repeated root of R makes BOTH sides of uR' = NR vanish, so it is a SPURIOUS saddle, not a critical point of log R - N log u. (1) THE CORRECTION IS STRONGER. Genuine saddles are roots of S(u) := uR'(u) - N R(u) that are NOT roots of R. If none existed, then since deg S = deg R = M+N with leading coefficients M r_d and r_d, every root of S being a root of R forces S = M R, i.e. uR' = (M+N)R, whose only solutions are R = c u^{M+N} — contradicting R(0) != 0 (the min exponent is exactly -N). **SO A GENUINE SADDLE ALWAYS EXISTS FOR M >= 1**, verified in every test case. (2) CONSEQUENCE: the entire remaining content of TNC is the ANALYTIC question of whether the dominant saddle's contribution can be cancelled by equal-modulus competitors — precisely boxeph's 'Bessel uniformity gap' — and this isolates it at ALL bidegrees simultaneously, where the Dickson ladder proceeds bidegree by bidegree. (3) The rising/falling-factorial split is now explicit: THM-1530(A)'s Gamma(Dm+1) is the RADIAL factorial, CT(Lambda^m) = [u^{Nm}]R^m is the factorial-free ANGULAR core, and the saddle values are its Fuss-Catalan growth rates"
status: PROVED (1); (0) is a refutation of my own first attempt; (2) is a reduction, not a closure
author: opus-2026-07-20-S416
depends_on: [THM-1530 (klein: the toral nullcone; Lagrange-Burmann at extreme weight), THM-1595 (boxeph: Dickson ladder closures), THM-1550]
---

# THM-1615 — The TNC's algebraic obstruction is empty

## 0. Refuted first: my own necessary condition

I proposed: *a TNC violator must have `disc(R) = 0`* — reasoning that the saddle value
`R(u*)/u*^N` must vanish, forcing `u*` to be a multiple root.

**False.** Counterexamples, computed:

| `Λ = u^{−N}R` | `disc R` | `gcd(R, uR′−NR)` | `CT(Λ^m)`, `m = 1..6` |
|---|---|---|---|
| `N=1`, `R=(1+u)²` | **0** | `u+1` | `2, 6, 20, 70, 252, 924` = `C(2m,m)` |
| `N=2`, `R=(1+u)⁴` | **0** | `(u+1)³` | `6, 70, 924, …` = `C(4m,2m)` |
| `N=2`, `R=(1−u)²(1+u)²` | **0** | `(u−1)(u+1)` | `−2, 6, −20, 70, …` |

All have vanishing discriminant *and* a nontrivial gcd, and **none** has a vanishing constant
term at any `m`.

**The error.** A repeated root of `R` makes **both** sides of `uR′ = NR` vanish. It is a root
of `S := uR′ − NR` for a trivial reason, not a critical point of `log R − N log u` — `log R`
is singular there. It is a **spurious saddle**. Recorded because "the discriminant vanishes"
is a natural first guess and it cuts nothing.

## 1. The correction, and it proves more (PROVED)

> **Definition.** A **genuine saddle** is a root `u*` of `S(u) := uR′(u) − N·R(u)` with
> `R(u*) ≠ 0`.

> **Theorem.** For `Λ = u^{−N}R(u)` with `deg R = M+N`, `R(0) ≠ 0`, and `M ≥ 1`, a genuine
> saddle always exists.

*Proof.* Suppose not: every root of `S` is a root of `R`. Compare degrees. The leading term
of `uR′` is `(M+N)r_d u^{M+N}` and of `NR` is `N r_d u^{M+N}`, so `S` has leading coefficient
`M·r_d ≠ 0` and `deg S = M+N = deg R`. Two polynomials of equal degree, one of whose root
multiset contains the other's, must be proportional: `S = M·R`. Hence

```
uR′ − NR = MR   ⟹   uR′ = (M+N)R   ⟹   R = c·u^{M+N}
```

(the ODE `uf′ = kf` has only pure powers as polynomial solutions). But `R(0) ≠ 0` forces
`M+N = 0`, contradicting `M ≥ 1`. ∎

**Verified** on every test case: `N=1, R=(1+u)²` → genuine saddle `u=1`; `N=2, R=(1+u)⁴` →
`u=1`; `N=2, R=1+u+u²+u³+u⁴` → `{1, −1, (−1±i√15)/4}`; `N=2, R=(1−u)²(1+u)²` → `{±i}`;
`N=3, R=(1+u)⁶` → `u=1`. The only case with none is `N=1, R=1+u`, where `M = 0` — the
one-sided case, already settled by THM-1530(B).

## 2. What this reduces the TNC to

At a genuine saddle the value `R(u*)/u*^N` is **nonzero**, so the Cauchy integral

```
CT(Λ^m) = [u^{Nm}]R^m = (1/2πi) ∮ R(u)^m u^{−Nm−1} du
```

carries an `m`-th power term with nonzero base. A single dominant nondegenerate saddle
therefore gives `CT(Λ^m) ≠ 0` for all large `m`, and the TNC holds.

> **So the entire remaining content of the TNC is:** can the dominant saddle's contribution
> be **cancelled by other saddles of equal modulus**?

That is exactly the **"Bessel uniformity gap"** boxeph flags in THM-1595. What §1 adds is
that **nothing algebraic stands in the way** — the obstruction is purely the multi-saddle
oscillation, and §1 establishes this **at every bidegree `(M,N)` simultaneously**, whereas
the Dickson ladder closes bidegrees one at a time (`N=1` all `M`; `(2,2)`, `(2,3)` proved;
`(2,4)`, `(3,3)` by elimination).

**This is a reduction, not a closure.** It does not close any new bidegree by itself, and it
should not be cited as if it did.

## 3. Where the rising and falling factorials actually sit

The owner asked how factorials bear on this. The split is now explicit and clean:

| layer | object | factorial content |
|---|---|---|
| **radial** | `Γ(Dm+1) = (Dm)!` (THM-1530 A) | the rising factorial — the Laplace/Gamma weight `e^{−s}ds`, `E[s^k] = k!` |
| **angular** | `CT(Λ^m) = [u^{Nm}]R^m` | **factorial-free** — a diagonal of a power |
| **growth** | saddle values `(R(u*)/u*^N)^m` | Fuss–Catalan / Raney rates; `C(2m,m)`, `C(4m,2m)` above are the `R=(1+u)^k` instances |

So factorials are the *radial* half and are already fully accounted for (mac-mini's THM-1600:
`L((av+b)^m) = m!·a^m·e_m(b/a)`, truncated exponential; `L((v−1)^m)` = derangement numbers).
**The TNC is precisely the factorial-free residue**, which is why no factorial identity will
close it — the content is the geometry of the saddle set, not the arithmetic of `(Dm)!`.

## 4. Open

1. **The multi-saddle cancellation.** Can `Σ_j c_j (ρ e^{iφ_j})^m` vanish for all `m` with
   the `c_j` algebraic prefactors? For finitely many equal-modulus saddles with distinct
   phases this is a vanishing-sums-of-roots-of-unity question — the natural next tool, and it
   is *not* what the ladder is doing.
2. Whether the saddle set can be forced to be a single point by a normalisation (`u`-scale
   and `t`-scale gauges are already used in THM-1595).

## Verification

`04-computation/tnc_saddle_discriminant_opus_S416.py` (the refutation; the `disc`/gcd table;
a 400-quartic sweep at `N=2` finding no violator),
`04-computation/tnc_saddle_corrected_opus_S416.py` (genuine saddles in every test case; the
`uf′ = kf` step). Outputs in `05-knowledge/results/`.
