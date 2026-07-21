# The NC2 proof transfers to LRC only after a Frobenius-safe certificate is built

**codex-2026-07-21-NC2-transfer.** This is the repo-wide search requested after
THM-2022. It cross-read the current LRC inverse-theorem frontier (S101--S107),
the relation-lattice/TTNC proposals (THM-1820/1825/1830), exact-period and
Ramanujan packet work (HYP-2978/2979/3036), the polynomial-method composite-14
audit (HYP-2758 and S31ag), p-adic/Hensel zippers (HYP-3024), moment-LP work,
the technique indices, and the broader problem zoo. The result is one proved
transfer theorem, three no-go statements, and a smaller LRC target.

## 1. The reusable mechanism is not "use p-adics"

The load-bearing unit in THM-2022 is:

```text
produce a nonzero base sum Q
  + expose every term that belongs to Q as one tied face
  + make every other term pay a finite-place cost
  + let Frobenius preserve Q as Q^p.
```

The creative move was refusing atomwise isolation. A tied layer is safe when
the **whole layer** is the object transported by Frobenius. This principle is
broader than Gaussian moments, but each application must separately supply a
base nonvanishing theorem and a filtration that removes other layers.

## 2. What transfers rigorously to LRC

THM-2041 proves the exact-period version. In characteristic `p`, with `p`
coprime to the phase period `q`, exponent multiplication permutes the cyclic
group. Primitive-support and Ramanujan kernels are invariant because
multiplication by `p` permutes `(Z/qZ)^*`. Therefore

```text
ct_q(Theta_q Lambda^(p*m)) = ct_q(Theta_q Lambda^m)^p.
```

This is the right formal content behind the old Paley/Frobenius and Ramanujan
carrier language. A whole exact-period packet can be amplified without phase
cancellation reappearing. The result remains valid at composite period 14 for
all coefficient primes coprime to 14.

That last sentence challenges a persistent conflation in the repo. The
Sungkawichai--Trakulthongchai polynomial method fails at `k+1=14` because it
tries to interpolate over `Z/14Z`, which is not a field and has nonunits and
null polynomials. THM-2041 instead uses a separate characteristic `p` while
leaving the phase period `q=14`. It preserves packets at composite period but
does not repair interpolation over `Z/14Z`.

## 3. Why this does not prove LRC(14)

The missing LRC input is now sharper than "signed cancellation."

### No DvdK base theorem in the finite phase algebra

The finite cyclic function algebra is a product of fields and contains many
zero divisors. For a covering configuration, the product of safe indicators is
pointwise zero for the exact reason one is trying to contradict. Frobenius
turns zero into zero. By contrast, DvdK guarantees that a Laurent polynomial
whose support surrounds zero has a nonzero constant term in some power. There
is no analogous automatic base `Q` for a cover.

### No unbounded factorial tower

The NC2 moment order can be multiplied by an arbitrarily large good prime.
LRC's danger count takes values only in `{0,...,13}`. Its factorial/Bonferroni
moments terminate; multiplying the intersection order by `p>13` makes it
identically zero. Thus Kummer must act on a denominator/period or packet lift,
not on the bounded danger-count order.

### A nonzero signed trace is not a lonely phase

Ramanujan projectors preserve exact period but forget the safe inequality and
can retain a nonzero signed trace on a covering packet. The old quotient audits
already found this empirically: raw Ramanujan traces mix routes until the
primitive safe deck, endpoint owner, and boundary/open label are restored.
THM-2041 proves the projector is stable; it also explains why that stability
alone is too weak.

## 4. The revised LRC attack

Use the existing compact/inverse-theorem reduction to select a labelled
medium-modulus packet, not a runner or an individual phase. Its proof carrier
should retain:

```text
exact period
primitive phase deck
safe-phase inequality
endpoint/open-boundary owner
p-adic carry or Hensel status
base integer safe count
terminal witness route.
```

Then aim for a base theorem saying the packet's integer safe count `D_q` is
nonzero modulo some `p>phi(q)`. This immediately implies `D_q>0`, and THM-2041
closes every p-power descendant whose twist stays exact-period invariant. In
the current S101--S104 language, this would turn a global medium-modulus
comparison into an exact phase-count certificate before invoking Freiman. It
does not assume the missing Diophantine-to-energy bridge.

The target is deliberately packet-local. A global assertion that every
covering residual has such a packet would be essentially LRC(14) again.

## 5. Transfer map across other repo problems

| front | usable part of THM-2022 | obstruction / next move |
|---|---|---|
| LRC exact-period / TTNC | THM-2041 preserves a whole primitive-period twist | prove a nonzero base safe count and retain pointwise labels |
| LRC p-adic/Hensel zippers | keep a complete unit-root or lift orbit, not one chosen root | singular roots and CRT coupling still need a filtration theorem |
| S--T finite lifts | coefficient-prime amplification can compress descendants | it does not fix interpolation over composite `Z/14Z` |
| LRC moment-LP / Bonferroni | exposed intersection packets remain useful | Kummer amplification of moment order is impossible because the alphabet stops at 13 |
| Paley/Frobenius tournament carriers | primitive residue orbits and Gauss/Ramanujan packets are stable | no implication from orbit stability to tournament-score or LRC noncancellation; MISTAKE-214 remains active |
| toral/single-character nullcones | whole character packets are exactly the natural Frobenius object | need a return-weight filtration if the functional is not Gaussian |
| polynomial/resultant and Casas--Alvero finite-place threads | use the whole initial form or resultant-unit layer, avoiding unique-term demands | base nonvanishing/resultant-unit hypotheses remain problem-specific |
| convolution, tiling, and finite cyclic covering counts | group-algebra projectors can collapse prime-power descendants | only integer counts with a size bound lift modular nonzero to positivity |
| higher Gaussian dimension | algebraic descent and whole-face philosophy remain valid | one scalar face height no longer controls a product of coordinate factorials |

The general lesson is not that Frobenius solves each problem. It says what a
legal transfer object looks like: a **whole invariant initial layer**, coupled
to a problem-specific theorem that makes the layer nonzero and a filtration
that removes everything else.

## 6. Assumption challenge and Tournament Analysis

Candidate vertices considered were runners, gaps, phases, residues, primitive
periods, Frobenius orbits, relation-lattice channels, projector kernels, and
proof obligations. Runners and individual phases fail because Frobenius
permutes them; the invariant object is the whole exact-period projector.

No tournament was imposed on THM-2041's verification. Its pairwise observable
is equality under a group automorphism, which is symmetric, not a binary
orientation. An arbitrary tie Hamiltonian path would destroy precisely the
whole-packet sum the theorem retains. For LRC route scheduling, proof carriers
may still be tournament vertices, ordered by preservation of the six labels
above; but that is a proof-engineering tournament, not part of the arithmetic
theorem.

## 7. External frontier check

The 2026 Sungkawichai--Trakulthongchai preprint proves `LRC(k)` for `k<=12`
(at most thirteen runners) and explicitly uses a polynomial shortcut only when
`k+1` is prime. Malikiosis--Santos--Schymura give the finite velocity bound;
Giri--Kravitz and Jain--Kravitz describe spectral/relative-spectral strata,
while Bedert's Riesz products improve the general asymptotic bound. These are
compatible with the map above: finite checking supplies compactness, spectral
work supplies strata, and harmonic products supply Archimedean positivity;
none supplies the missing exact-period base safe count for the composite-14
residual.

Primary sources checked in this sweep:

- Sungkawichai--Trakulthongchai, *The lonely runner conjecture is true for at
  most 13 runners*: <https://arxiv.org/abs/2604.23906>.
- Malikiosis--Santos--Schymura, *The Lonely Runner Conjecture turns 60* (finite
  checking and explicit finite reductions): <https://arxiv.org/abs/2411.06903>.
- Giri--Kravitz, *The lonely runner spectra*: <https://arxiv.org/abs/2304.01462>.
- Jain--Kravitz, *Relative lonely runner spectra*: <https://arxiv.org/abs/2411.12684>.
- Bedert, *The lonely runner conjecture via Riesz products*:
  <https://arxiv.org/abs/2511.16636>.
