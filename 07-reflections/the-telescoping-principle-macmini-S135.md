# The telescoping principle: forbidding one side of a grading collapses the problem

*mac-mini-2026-07-20-S135. The owner, working GMC(2), said: "think forbidding one variable
telescopes, see how that pattern can abstractly apply to other problems the repo has touched."
The hint was exactly right, it closed a branch of GMC(2) (THM-1520), and on looking around
the repo the same shape turns out to be doing quiet work in at least four other places.*

---

## The principle

> **Suppose an invariant is computed by a functional that is GRADED — it annihilates every
> nonzero grade. If the object's support in that grading is ONE-SIDED, then the only
> surviving contribution takes the grade-0 piece from every factor, so the invariant
> COLLAPSES to a function of the grade-0 part alone, in strictly fewer variables.**
>
> **Corollary, and this is the useful half: failure REQUIRES two-sided support.**

The mechanism is trivial once stated — grades add under multiplication, so a one-sided
support can only reach `0` by being `0` everywhere. What makes it worth naming is how often
the hard case of a problem turns out to be exactly the two-sided one, and how often a
"parity obstruction" in this repo is this principle wearing a different coat.

## Where it just paid: GMC(2)

For one complex Gaussian, `E[Z^aW^b] = a!δ_{ab}`, so with **charge** `c = deg_Z − deg_W` the
expectation kills every nonzero charge. (Cleanly: `Z = √V e^{iθ}` with `V ~ Exp(1)` and `θ`
uniform *independent*, and charge is the `θ`-Fourier index — the grading is literally the
angular decomposition.)

One-sided charge support therefore gives `E[Pᵐ] = L(p₀ᵐ)` where `L(vᵏ) = k!` — **two variables
collapse to one**. And the one-variable question has a clean answer: `L(pᵐ)/(a_Dᵐ(Dm)!)`
converges to `exp(a_{D−1}/(D a_D))`, never zero, so `L(pᵐ) = 0` for all `m` forces `p = 0`.
Branch closed. What remains of GMC(2) is *exactly* the two-sided case (THM-1520).

## Four places the repo already had it

**Ordinal sums and arborescences (THM-1460).** `T₁ ⊕ T₂` forbids arcs in one direction between
the blocks. That makes the in-Laplacian **block upper triangular**, so
`spec = spec(L₁) ∪ (|T₁| + spec(L₂))` and the arborescence count telescopes into the two
factors: `Σa(T₁⊕T₂) = Σa(T₁)·det(|T₁|I + L₂)`. Same principle, grading = which block a vertex
lies in, one-sidedness = the ordinal sum itself. The *reason* `log Σa` picks up a
size-dependent shift while `log H` does not is that the two invariants read the collapse
differently.

**Skew-Seidel parity (THM-1440).** `Sᵀ = −S` gives `p(−x) = (−1)ⁿp(x)` — a `ℤ/2` grading. At
odd `n` the characteristic polynomial is an *odd* function, so its grade-0 coefficient must
vanish: a **forced zero eigenvalue**. That is this principle at its smallest, and it is the
same fact as "`sin` vanishes at 0 and `cos` does not."

**The master cycle-packing polynomial (THM-506 / HYP-2514).** `Φ(T;{y_k})` has the spectrum as
its *signed vertex-graded* face and `H` as its *unsigned odd-only* face. **Faces are graded
pieces**, and each named invariant is a one-sided restriction of the full sum. The reason the
spectrum is computable and `H` is not is which face collapses to a determinant.

**Cut ⊕ cycle (THM-1405 / THM-1420).** Star flips are the grade-0 "gauge" directions; the
holonomies are what survives the quotient. And the bicycle space vanishing exactly at odd `n`
(THM-1440 §C) is *why* the splitting is canonical there — one-sidedness again, over `F₂`.

## What the pattern is good for

Three concrete uses, in decreasing order of how sure I am:

1. **As a triage step.** Given a vanishing claim, find the grading first. If the support is
   one-sided the claim is usually trivial; if two-sided, that is where the content is. In
   GMC(2) this cut a genuinely open problem down to a precisely delimited remainder in about
   an hour.
2. **As a source of necessary conditions.** "Failure requires two-sided support" is the kind
   of constraint that narrows a search enormously and costs nothing to prove.
3. **As a warning about counterexample hunting.** The `n ≥ 3` GMC counterexamples all have
   two-sided charge; a search restricted (even accidentally) to one-sided objects will find
   nothing and *look* like evidence. That is a specific way to fool yourself, and I nearly
   did it in S133 with a search whose positive control failed.

## The honest limit of the analogy

These are the same *shape*, not the same theorem. The gradings are different in kind — `ℤ` by
charge, `ℤ/2` by parity, a block decomposition, an `F₂` subspace — and nothing transfers
automatically between them. What transfers is the **question to ask**: *what is the grading,
and which side of zero does the support live on?* Every time the repo has recorded something
as vanishing "for parity reasons" or "because the charges do not balance," this is what was
happening, and it is worth naming so the next instance is recognized in minutes rather than
rediscovered.

One caution against over-reading: the principle explains why certain things are *easy*, and it
delimits where the hard case is. It has not, in any of these five instances, solved the hard
case. GMC(2)'s two-sided branch is still open; the arborescence extremals are still
unproved; `H` is still not a determinant.

---

*Cross-links: THM-1520 (GMC(2), the one-sided branch), THM-1500 (the GMC master theorem),
THM-1460 (ordinal sums), THM-1440 (skew-Seidel parity, and the bicycle space at odd n),
THM-506 / HYP-2514 (the master cycle-packing polynomial), THM-1405 / THM-1420 (cut ⊕ cycle).
Artifacts: `04-computation/gmc2_charge_telescoping_macmini_S135.py`,
`gmc2_limit_and_twosided_macmini_S135.py` (+outs).*
