# JC(2): the resonance dictionary, and the manufactured-valuation orbit-product as the natural tool for the multi-root residual

**death-star-2026-07-22-S107** (HYP-8940). Owner: work the 2-D Jacobian Conjecture; try creative new routes; explore
unrelated threads for ideas. **JC(2) is open and famous ("the graveyard of false proofs"); this is exploration,
not a proof.** The fleet has already assembled the main routes (boxeph S225: JC(2)⟺VC(4), descent-termination,
coprime-interval/Lamé; boxeph S146: the H-template / Keller degree-monoid / Euler ledger / ODD-DEGREE conjecture;
codex: Poisson suspension / DC(2), THM-2045/2063; klein/mac-mini: Euler–Zariski, golden corner). This note adds
two things built on those: (1) a **verified "resonance dictionary"** unifying the tame/hard split of *three* repo
problems (DvdK, the Hessian nullcone, JC(2)) as one *single-root vs multi-root* phenomenon; and (2) a concrete
**new attack route** on the multi-root residual — transferring last session's brand-new S106 tool (the
manufactured-valuation Galois orbit-product of THM-2067) to the places at infinity of a JC(2) counterexample,
identifying the valid reduction steps and the one open crux.

## 1. Verified anchors (script `jc2_resonance_dictionary_leadingform_deathstar_S107.py`)

Write a polynomial map `F = I + H`, `H = (H₁,H₂)`. Then `det(DF) = 1 + div(H) + jac(H₁,H₂)` (exact identity,
verified). For a Keller map `det(DF)=1`:

- **Homogeneous 2-D Keller maps are exactly shears.** If `H` is homogeneous of degree `d ≥ 2`, then `div(H)`
  (degree `d-1`) and `jac(H₁,H₂)` (degree `2d-2`) are different degrees, so `det(DF)=1 ⟺ div(H)=0` **and**
  `jac(H₁,H₂)=0`. For same-degree forms `jac(H₁,H₂)=0 ⟺ H₁,H₂ proportional` (verified), and `div(H)=0` then forces
  `H = c·(bx−ay)^d·(a,b)` — a rank-1 shear along `(a,b)`, with invariant `bx−ay` and inverse `I−H` (verified).
  **So homogeneous 2-D JC is trivial; the entire difficulty is non-homogeneous.** (This is why the
  Bass–Connell–Wright cubic-homogeneous reduction *raises dimension* — it has to.)
- **Leading forms are powers of a common binary form.** For a Jacobian pair `(f,g)` with `deg f = n`, `deg g = m`,
  the top part of `jac(f,g)=1` gives `jac(f_n,g_m)=0`, hence `f_n = α h^p`, `g_m = β h^q` for a common form `h`
  (`n=pk, m=qk, k=deg h`), verified via the functional-dependence relation `f_n^{q} = const·g_m^{p}`.
- **The descent dichotomy.** Killing the common leading form reduces degree **unless neither degree divides the
  other**: `n=m` ⟹ proportional ⟹ `g−(β/α)f` drops `deg g`; `m∣n` ⟹ `f−c·g^{n/m}` drops `deg f` (verified,
  4→3); `n∤m ∧ m∤n` is the stuck stratum (Abhyankar–Moh). This is the arithmetic (coprime/gcd) obstruction boxeph
  S225 identified — the Lamé/Fibonacci-worst-case descent.
- **Places at infinity = distinct roots of `h`.** A Keller pair's two components share the same places at infinity
  (roots of the common `h`). **Single-root `h = ℓ^k` ⟹ one place at infinity ⟹ Abhyankar–Moh ⟹ tame** (the affine
  curve embeds as a coordinate line). **A JC(2) counterexample requires a multi-root `h` (≥2 distinct roots).**

## 2. The resonance dictionary (the cross-thread synthesis)

The three hardest "nullcone/degeneracy" problems in the repo have the **same** tame-vs-hard divider — a *single
root* (one direction) versus *multiple coincident roots* (resonance):

| | TAME pole (single-root) | RESONANT pole (multi-root) |
|---|---|---|
| **DvdK / GMC2** (S101, S106) | unique minimal cycle ⟹ single-term `CT ≠ 0` | ≥2 coincident cycles ⟹ cancellation; needs THM-2067 |
| **Hessian nullcone** (S103) | rank-1 nilpotent ⟹ `P` one-sided `(x+iy)^d` | rank-≥2 nilpotent ⟹ dim-≥3 counterexample (THM-1300) |
| **JC(2)** (this note) | single-root leading `h = ℓ^k` ⟹ 1 place at ∞ ⟹ tame | multi-root `h` ⟹ ≥2 places at ∞; the open residual |

In all three, the tame case is *elementary/classical* and the hard case is where the problem lives. The shared
mechanism: a **single dominant direction** (one cycle / one isotropic line / one place at infinity) leaves a lone
uncancellable term, while **coincident directions** open the door to cancellation/counterexamples. This is not a
formal reduction between the problems (MISTAKE-229: no proved `NC2 ⇒ GMC2 ⇒ JC2` chain; de Bondt doubling lands
VC in dim 4, off the rank-1 wall) — it is a **structural dictionary** that says the same *tool* should attack all
three resonant residuals.

## 3. The new route: manufactured-valuation orbit-product at the places at infinity

Last session (S106) mapped the DvdK resonant residual to codex's THM-2067, whose engine is a **manufactured
valuation**: introduce a spectral parameter `t`, form `Φ(X)=X^M − tR(X)`, and derive a contradiction from the
Galois **orbit-product** — the small-root product `Π(t)=ct` has `t`-adic valuation 1, while the full Vieta norm
has valuation 0, so transitivity forces `r·1 = η·0`, impossible. The key move is *creating* a valuation where the
coefficients (units, no natural valuation) offered none.

The JC(2) multi-root residual has the **dual situation**: the places at infinity already carry *natural*
divisorial valuations, but over `ℂ` the roots of `h` have *no Galois structure* to product over — exactly the
obstruction that made DvdK's unit-coefficient case hard. The transfer:

1. **(Lefschetz reduction — valid.)** JC(2) is first-order over an algebraically closed char-0 field; a `ℂ`
   counterexample descends to `ℚ̄`. So WLOG `f,g ∈ ℚ̄[x,y]` and `h ∈ ℚ̄[x,y]`.
2. **(Galois orbit at infinity — valid.)** The ≥2 distinct roots of `h` (the places at infinity of `f=0`) lie in a
   number field; `Gal(ℚ̄/ℚ)` permutes them in orbits. The multi-root residual thus acquires exactly the
   Galois-orbit structure THM-2067 needs — supplied by the arithmetic of the counterexample, not manufactured.
3. **(The crux — open.)** The Jacobian condition `jac(f,g)=1` constrains the *local data* at each place
   (the pole orders and leading coefficients of `g` along the branches of `f`, i.e. the ramification/different at
   infinity). The route seeks a **product of local invariants over a Galois orbit of places** that is rational
   (Galois-invariant) yet carries a valuation incompatible with the global `jac=1` constraint — the direct
   analogue of `Π(t)=ct` vs. the Vieta norm. **Identifying that invariant is the open problem** this route poses.

**Why this is the right shape.** THM-2067's orbit-product is a *reciprocity/norm law*. boxeph S146's Euler ledger
`Σᵢ(d−kᵢ)χ(Cᵢ)=d−1` is a **topological reciprocity** (χ-additivity over the covering). The proposed orbit-product
is its **arithmetic refinement** (a valuation/norm reciprocity over the places at infinity), with THM-2067 as the
exact template. Weil reciprocity (product of local symbols over all places `= 1`) is the classical law of this
type; the JC(2) content would be that a counterexample's places at infinity violate the reciprocity forced by
`jac=1`. This reframes "descent termination" (boxeph S225's coprime-interval *feasibility* engine) as a
*cancellation/valuation* obstruction — the strictly finer S106 tool — localized at the places at infinity.

## 3½. A parallel tool: the monomial/Nullstellensatz certificate (boxeph S231)

Same-day, boxeph S231 (HYP-8932) discharged the DvdK residual for every small support by a **monomial
certificate** — DvdK-free ⟺ the ideal `⟨CT(f^m)⟩ₘ` contains a torus-unit monomial `∏x^e = Σ g_m·CT(f^m)`, an
exact `ℚ`-linear-algebra check per support (no Galois), formalized kernel-pure. This is the *effective/certificate*
face of the same coin: THM-2067 gives the non-effective general theorem (Galois orbit-product), while the monomial
certificate discharges any *given* support finitely, the degree bound being the open effective-DvdK content.

The JC(2) analogue is immediate and worth stating: `F` invertible ⟺ `ℂ[f,g]=ℂ[x,y]` ⟺ `x,y ∈ ℂ[f,g]` ⟺ an
**ideal-membership certificate** `x = P(f,g)` exists — again a Nullstellensatz statement whose *existence* is the
conjecture and whose *degree bound* is the effective (Lamé/descent-termination) crux. So the dictionary gains a
third, unified column — **certificate exists ⟺ tame; effective degree bound = the open part** — across DvdK
(monomial certificate, S231), and JC(2) (coordinate certificate `x=P(f,g)`). The orbit-product (§3) and the
certificate are the *non-effective* and *effective* attacks on the same resonant residual; both are now on the
table for JC(2), the certificate via `ℂ[f,g]`-membership and the orbit-product via the Galois places at infinity.

## 4. Honest scope

- **Not a proof, and not claimed to close any stratum.** §1 anchors are classical (Abhyankar–Moh, van der Kulk),
  here verified in-repo and assembled; §2 is a structural dictionary, not a formal implication chain; §3 is a route
  *sketch* whose crux (the orbit invariant) is explicitly open.
- **What is genuinely new:** the unification of the three residuals as one single-root/multi-root resonance, and
  the observation that S106's manufactured-valuation orbit-product is the natural tool for the JC(2) multi-root
  residual, with the Lefschetz→Galois-places reduction making its two structural hypotheses (a Galois orbit; a
  norm to oppose) available — the same two quantities THM-2067 §5 notes the LRC packet orbit is *missing*.
- **Consistent with the walls:** this does not contradict GMC(2)⇏JC(2) (S205) — it does not route JC(2) through
  GMC's Frobenius wall; it transplants only the *orbit-product tool*, to the places at infinity, where the
  arithmetic (not the rank) supplies the orbit.

Cross-links: S103 (planar-JC = NC2 one-sidedness, the Hessian row of the dictionary), S106 (the manufactured-
valuation orbit-product this transfers), S101 (unique-cycle, the DvdK row), boxeph S225 (descent-termination /
coprime-interval — this sharpens its feasibility engine to a valuation one), boxeph S146 (Euler ledger =
topological reciprocity, the ODD-DEGREE conjecture — the proposed orbit-product is its arithmetic refinement),
boxeph S205 (GMC⇏JC rank wall — respected here), codex THM-2045/2063 (Newton-edge non-cancellation obstruction /
tame one-fiber-linear class — the single-root pole), THM-2067 (the orbit-product template), THM-1300 (dim-3
counterexample, the Hessian rank-≥2 pole), memories `nc2-gmc2-lean-formalization-state`,
`gmc-lrc-same-positivity-manoeuvre`. Script `04-computation/jc2_resonance_dictionary_leadingform_deathstar_S107.py`
(+ `.out`). HYP-8940.
