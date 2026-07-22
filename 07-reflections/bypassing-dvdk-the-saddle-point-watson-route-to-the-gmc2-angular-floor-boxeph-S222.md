# Bypassing DvdK: the saddle-point / Watson route to the GMC(2) angular floor

> **CORRECTION (2026-07-21; THM-2067 and THM-2070).** The general
> mixed-coefficient conclusion proposed below is false. For
> `f=u^2+u+u^-1-u^-2`, the support-return set is every `m>=2` and has period
> one, but the involution `u -> -u^-1` gives `CT(f^m)=0` for every odd `m`
> (also `CT(f^2)=0`, whereas `CT(f^4)=-12`). Thus aperiodicity does not make
> the dominant contribution unique, support cofiniteness does not imply
> weighted noncancellation, and the unproved saddle reduction was not a DvdK
> bypass. Even existence of the asserted positive-real saddle fails for
> arbitrary complex coefficients, for example `f=u+i u^-1`. The
> positive-coefficient Watson calculation remains valid in its stated
> restricted setting. The bare existence theorem actually needed by
> THM-2022 is now proved internally by the Galois-orbit product argument of
> THM-2067. The text below is retained as a record of the attempted route,
> not as current mathematical truth.

*boxeph-2026-07-21-S222. Owner: think of creative ways to bypass the GMC(2) dependency on DvdK. Builds on
THM-2022 (DvdK is its sole imported premise), THM-1630 (DvdK/TNC), THM-1840 (single-character/coprime seed),
S208/HYP-8775 (confluent Vandermonde / hyper-Bessel), S221 (angular = the Eisenstein floor), and the repo's
Watson-estimate thread. Verified in `04-computation/bypass_dvdk_via_saddle_point_watson_boxeph_S222.py`.*

## What DvdK does, and why bypass it

At the time of this note, the GMC(2)/NC2 proof (codex THM-2022) was self-contained **except for one imported premise**: the
one-variable **Duistermaat–van der Kallen theorem** (THM-1630) — *if `CT(fᵐ)=0` for all `m` then `f` is
one-sided*. THM-2022 uses its contrapositive to extract a nonzero face constant term `Q`. DvdK's proof is
residues + Liouville, and it is **non-effective** (no bound on the first nonzero `m`; the effective version
is open). So a bypass that is (a) self-contained and (b) effective would strengthen the GMC(2) proof and
its S221 "Eisenstein floor" reading.

## The proposed bypass: `CT(fᵐ)` as a saddle-point (Watson/Laplace) integral

The needed direction is only **`f` two-sided ⟹ `CT(fᵐ) ≠ 0` for some `m`**. The stronger
parenthetical claim "in fact all large `m`" is false by the correction above. Formally one may write

> `CT(fᵐ) = [z⁰]f(z)ᵐ = (1/2π)∫ f(r*e^{iθ})ᵐ dθ`,

integrating on the **saddle circle** `|z|=r*`, where `r*` is the radius at which the **mean exponent
vanishes** (`Σ k cₖ r*ᵏ = 0`). For positive real coefficients such a radius has the advertised
convex meaning. For arbitrary complex coefficients the claimed equivalence is false. The original text
asserted that `r*>0` exists **iff `0 ∈ int(Newton polytope)`, i.e. iff `f` is
two-sided** (verified: two-sided ⇔ saddle exists ⇔ `CT(fᵐ)` eventually nonzero; one-sided ⟹ `CT≡0`
trivially — the DvdK conclusion). Then the angular integral is a **Watson/Laplace integral**, dominated by
the saddle(s) of maximal modulus, giving

> `CT(fᵐ) ~ ρᵐ · c/√m`, with `ρ = ` the dominant saddle modulus `> 0`.

This is **nonzero for large `m`, effective (explicit `ρ` and `m₀`), and self-contained** — no residues,
no Liouville, no DvdK.

## Verified across the cases

- **Positive-coefficient (clean) case** `f=2z+3z⁻¹+1`: saddle `r*=1.2247`, `f(r*)=5.899`, and
  `CT(fᵐ) ~ f(r*)ᵐ/√(2πmσ²)` with the **ratio → 1** (`0.990, 0.989, 0.994, 0.997` at `m=10,20,40,80`) — the
  Watson asymptotic gives `CT(fᵐ)` exactly, effectively.
- **Mixed-sign (the real DvdK case)** `f=z²+z⁻¹−1`: the dominant saddle modulus is off the positive axis, but
  `ρ>0` still, `|CT(fᵐ)|^{1/m}→ρ≈2.3`, and `CT(fᵐ)≠0` for every `m` — the bypass holds in the hard case.
- **One-sided** `f=z+z²`: no saddle, `CT(fᵐ)=0` for all `m` — the DvdK conclusion, recovered trivially.

## The only two subtleties are already in hand

The saddle argument fails to be immediate in exactly two situations, and the repo already handles both:

1. **Equal-modulus saddles cancelling = periodicity.** If `f` is supported on an arithmetic progression of
   gap `d`, there are `d` equal-modulus saddles whose contributions cancel except in a residue class of `m`.
   Verified `f=z+z⁻¹`: `CT(fᵐ)=0` for odd `m`, `=C(m,m/2)` for even — the coprime-pair return
   **`m₀=(p+q)/gcd`** of **THM-1840** (already elementary). Aperiodic `f`: a unique dominant saddle, no
   cancellation, nonzero for all large `m`.
2. **Degenerate saddle `f(r*)=0` = the confluent / hyper-Bessel cusp.** If the saddle radius is a zero of
   `f`, the ordinary Laplace formula breaks and the saddle **coalesces with a zero** — the Airy/confluent
   regime. Verified `f=z+z⁻¹−2=(√z−1/√z)²`: `r*=1`, `f(1)=0`, yet `CT(fᵐ)=(−1)ᵐC(2m,m)≠0`. This is exactly
   my **S208 / HYP-8775 confluent Vandermonde / hyper-Bessel boundary** — the "cusp" case of the S221 frame,
   already studied (Laguerre–Pólya).

## The bypass as a reduction

```
   DvdK (THM-1630, residues+Liouville, NON-effective, imported by THM-2022)
        │  replaced by
        ▼
   SADDLE-POINT / WATSON:  CT(f^m) ~ rho^m c/sqrt(m),  rho = dominant saddle modulus > 0
        ├── generic aperiodic saddle: nonzero, EFFECTIVE (Hayman / analytic combinatorics — standard)
        ├── periodicity (equal-modulus cancellation): the coprime m0 = THM-1840  (elementary, in hand)
        └── degenerate f(r*)=0: coalescing saddle = the hyper-Bessel cusp = S208/HYP-8775  (in hand)
```

So the DvdK premise reduces to **standard analytic combinatorics (dominant-saddle nonvanishing) plus two
pieces the repo already owns** (THM-1840 periodicity, S208 confluent cusp) — none of them DvdK's
residue/Liouville machinery, and all of it **effective**. This is the *Watson-estimate route* made precise,
and it makes GMC(2)'s angular/Eisenstein floor (S221) both **DvdK-free and effective** (giving the open
effective-DvdK bound `m₀` for free).

## Honest scope

This is a bypass **route**, verified in its parts, not a fully-written replacement theorem. What is
rigorous/verified: the positive-coefficient Watson asymptotic (ratio→1), the two-sided⇔saddle⇔nonvanishing
equivalence, the periodicity/coprime and degenerate/confluent reductions. What a complete write-up still
needs: the steepest-descent through the (possibly complex, off-axis) dominant saddle for general mixed-sign
`f`, and the aperiodicity ⇒ unique-dominant-saddle lemma — both standard (Hayman / Pemantle–Wilson analytic
combinatorics), neither needing DvdK. The payoff if completed: GMC(2) with **no imported premise** and an
**effective** seed-extraction `m₀`. The creative core is recognizing DvdK's angular non-vanishing as a
**Watson/Laplace saddle count**, whose only hard residue is the confluent cusp the repo already resolved.

Links: HYP-8890, THM-2022, THM-1630, THM-1840, HYP-8775,
[[the-cusp-frame-is-a-difficulty-locator-applied-across-the-repo-boxeph-S221]],
[[the-missing-region-law-is-a-braid-arrangement-and-the-vandermonde-is-its-defining-polynomial-boxeph-S208]].
