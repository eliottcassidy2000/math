---
id: THM-2142
title: "THE HALF-ANGLE BRIDGE: the a/b monoid ⟨a=x+1, b=x/2⟩ is the trigonometric skeleton shared by the tournament skew spectrum AND the GMC constant-term projection CT_u. The un-pursued gem the owner's frame set up: b(a(cos θ)) = (1+cos θ)/2 = cos²(θ/2) EXACTLY — the tournament's own coordinate-change b∘a (the {−1,0,1}↔{0,½,1} half-dictionary, THM-1555/1560, = A=(J−I+S)/2) evaluated on a cosine IS the trigonometric half-angle. Forward: a^n folded by b gives the even/odd companions E_n=b(a^n+ā^n)=char_S(transitive), O_n=b(a^n−ā^n), whose spectra are cot of roots of unity and whose eigenvalue equation is the multiple-angle formula (THM-1880/1575); backward: iterating b is dyadic angle-bisection (BS(1,2), THM-1885). The NEW content is the bridge to GMC(2)/directive-2: the charge-0 projection CT_u (constant term = circle average = the trigonometric functional (1/2π)∫f(e^{iθ})dθ) has two-charge value CT_u[(u+u⁻¹)^m]=C(m,m/2) (central binomial = power-reduction charge-0 term) and CT_u[(u^p+u⁻q)^{p+q}]=C(p+q,p) (the UNIQUE balanced composition, THM-1840). These are the SAME binomial family {C(n,2j)} that are the coefficients of E_n and that carry the triangular number C(n,2)=T_{n−1}=#arcs as E_n's subleading coefficient — split by the odd/even (charge-parity) involution x↦−x = the SL₂ Weyl axis (THM-1810). So a shifts, b halves, and folded they produce simultaneously (i) the tournament skew spectrum with triangular-number arc count, (ii) the half/multiple-angle cyclotomic ladder, (iii) via CT_u the two-charge DvdK non-vanishing = the base case of the GMC(2) moment engine. One monoid, read three ways; the owner's 'functional' a/b frame is the trigonometric skeleton of both the tournament and the GMC charge-0 projection."
status: >
  UNIFYING FRAME (verified/classical identities assembled), not an open-problem advance.
  ALL VERIFIED-EXACT in sympy (04-computation/ab_halfangle_ct_bridge_klein_S401.py, n=1..8 and
  (p,q) census): (1) the a/b foundation of THM-1880 — E_n²−O_n²=(x²−1)^n, subleading(E_n)=C(n,2)=
  T_{n−1}, E_n coefficients = even binomials C(n,2j); (2) the half-angle identity b(a(cos θ))−
  cos²(θ/2) = 0; (3) CT_u[(u+u⁻¹)^m] = C(m,m/2) and two-charge CT_u[(u^p+u⁻q)^{p+q}] = C(p+q,p).
  Two honest corrections/traps recorded below: THM-1880(2)'s coupled recursion is stated with
  E_n, O_n SWAPPED on the left (the correct crossed form is O_n=E_{n−1}+x O_{n−1}, E_n=O_{n−1}+
  x E_{n−1}); and CT_u is the FULL circle average and does NOT equal a finite (p+q)-root-of-unity
  average (that sum aliases the higher charges of Λ^m and returns 2^m — the cyclotomy lives in the
  RETURN TIME m₀=(p+q)/gcd, not in a collapse of CT_u). Proves no open problem; it is the
  functional bridge between the a/b archeology (directive 1) and the GMC(2) moment engine
  (directive 2).
source: klein-2026-07-22-S401 (owner: think trigonometric functions, and triangular numbers and thus tournaments as a(x)=x+1 and b(x)=x/2 composed recursively, think functionally; + finish the GMC(2) formalization)
depends_on:
  - THM-1880    # the a/b functional frame: E_n=b(a^n+ā^n), Pell, cot spectra, C(n,2)=T_{n-1} coeff
  - THM-1840    # two-charge DvdK is cyclotomic-clean: CT of (u^p+u^-q)^{p+q} is a single C(p+q,p)
  - THM-1555    # the half-dictionary b∘a=(x+1)/2 = A=(J−I+S)/2 (kind-pasteur)
related:
  - THM-1560    # the halving dictionary and the mod-2 collapse (klein): the ½ is what dies mod 2
  - THM-1575    # the circulant skew spectrum is a tangent: eigenvalue eqn = multiple-angle formula
  - THM-1810    # binary forms / SL₂ Weyl axis x↦−x = odd/even = the charge-parity split
  - THM-1885    # ⟨a,b⟩ ≅ BS(1,2)⁺; backward b = dyadic angle-bisection
  - THM-293     # the succession GF W: tangent numbers A000182 on the odd side (dormant, revival lead)
script: 04-computation/ab_halfangle_ct_bridge_klein_S401.py (+ .out)
---

# THM-2142 — the half-angle bridge: a/b is the CT_u cyclotomic skeleton

The owner's frame — **build tournaments from `a(x)=x+1` and `b(x)=x/2`, composed recursively;
think functionally and trigonometrically** — and the paired directive to **finish the GMC(2)
formalization** meet at one identity the corpus set up across THM-1880/1555/1575/1840 and never
wrote down: the tournament's own coordinate-change `b∘a`, evaluated on a cosine, **is** the
trigonometric half-angle, and the same functional (the constant term = circle average) is the
base case of the GMC(2) moment engine.

## 1. The half-angle identity (the un-pursued gem)

`b∘a = (x+1)/2` is the half-dictionary (THM-1555/1560): the entrywise map `{−1,0,1}↦{0,½,1}` that
turns the skew ±1 Seidel matrix `S` into the 0/1 adjacency matrix `A = (J−I+S)/2`. Evaluate it on a
cosine:

> **`b(a(cos θ)) = (1+cos θ)/2 = cos²(θ/2)`** — verified (the difference simplifies to `0`).

So the coordinate change between the sign world (`S ∈ {−1,+1}`) and the affine world
(`A ∈ {0,1}`) is, on the unit circle, **angle bisection** (`cos²` of the half-angle). This is why
`b` is "the ½ everywhere": the tiling fiber fraction, the Legendre exponent, the `Re=−½` line, and
the trigonometric half-angle are one generator `b`.

- **Forward (`a^n` folded by `b`) = multiple-angle.** `E_n = b(a^n+ā^n) = ((x+1)^n+(x−1)^n)/2 =
  char_S(transitive)` and `O_n = b(a^n−ā^n)`; their eigenvalue equation is the tangent/cotangent
  multiple-angle formula (THM-1575: substituting `x=i·tan θ` into the rotational `char_S` gives
  `i·sin(nθ)/cosⁿθ`), roots at `cot` of roots of unity.
- **Backward (iterating `b`) = dyadic angle-bisection.** `b` is the `BS(1,2)` generator (THM-1885):
  repeated halving is the dyadic solenoid / Stern–Brocot descent on the angle. (This backward
  direction is the one the corpus left untaken; it is named here, not yet developed.)

## 2. The constant-term projection `CT_u` is the trigonometric functional

The GMC(2) moment engine's charge-0 projection is `CT_u[f] = [u^0]f = (1/2π)∫₀^{2π} f(e^{iθ}) dθ`.
On `Λ = u + u⁻¹ = 2cos θ` it is the **power-reduction charge-0 term**:

> `CT_u[(u+u⁻¹)^m] = C(m, m/2)` (central binomial; `0` for odd `m`).

and on a genuine two-charge Laurent binomial it is a single multinomial at the return time
(THM-1840, the cyclotomic-clean case):

> `CT_u[(u^p+u⁻q)^{p+q}] = C(p+q, p)`, nonzero at the return time `m₀=(p+q)/gcd(p,q)`.

Verified across `(p,q) ∈ {(1,1),(1,2),(2,3),(3,5),(2,2),(1,3),(2,4)}`.

## 3. The shared binomial family — the bridge directive-1 ↔ directive-2

| object | binomials | where |
|---|---|---|
| tournament skew char `E_n` | `C(n,2j)` (even index), subleading `C(n,2)=T_{n−1}=#arcs` | THM-1880(5) |
| `CT_u` two-charge value | `C(m,m/2)`, `C(p+q,p)` | §2, THM-1840 |

Both are **the binomial coefficients cut out of the shift `(x+1)^n` / the character sum by the
odd/even (charge-parity) involution `x ↦ −x`** — the `SL₂` Weyl axis of the characteristic binary
form (THM-1810). `E_n` is the even half (secant/cotangent side, `char_S`); the odd companion `O_n`
is the tangent side; `CT_u` reads the charge-0 (even, balanced) part. So:

> `a` **shifts, `b` halves; folded, they produce simultaneously**
> **(i)** the tournament skew spectrum, whose second symmetric function is the **triangular number**
>   `C(n,2)=T_{n−1}=#arcs`;
> **(ii)** the **half/multiple-angle** cyclotomic ladder (`cos²(θ/2)` backward, `sin(nθ)` forward);
> **(iii)** via `CT_u`, the **two-charge DvdK non-vanishing** `C(p+q,p) ≠ 0` — the base case of the
>   GMC(2) moment engine (THM-1840, and the premise-free `GMC2DvdKTwoCharge.dvdk1_pair` in Lean).

The owner's "functional" a/b frame is the **trigonometric skeleton of both the tournament and the
GMC charge-0 projection.** The triangular number `T_{n−1}` (the tournament's arc count) and the
DvdK moment `C(p+q,p)` (the GMC base case) are two readings of the same shift-and-halve binomial
structure, separated only by the charge-parity involution.

## 4. Two honest corrections (measured, then recorded)

- **THM-1880(2) recursion labels.** The coupled recursion is stated there as
  `E_n = E_{n−1}+x·O_{n−1}, O_n = O_{n−1}+x·E_{n−1}`. The correct **crossed** form (verified n=1..8)
  swaps the left sides: `O_n = E_{n−1}+x·O_{n−1}`, `E_n = O_{n−1}+x·E_{n−1}` (each step crosses,
  since `E_{n−1}+x·O_{n−1} = ((x+1)^n−(x−1)^n)/2 = O_n`). The Pell identity, the `T_{n−1}`
  subleading coefficient, and the even-binomial coefficients are all correct as stated; only the
  recursion's E/O labels are transposed.
- **The cyclotomic collapse trap (avoided).** `CT_u` is the **full** circle average; it is **not**
  equal to a finite `(p+q)`-root-of-unity average `(1/N)Σ_{ζ^N=1}Λ(ζ)^m` (which aliases the higher
  charges of `Λ^m` and returns `2^m`-type values, not `C(p+q,p)`). The cyclotomic content is the
  **return time** `m₀=(p+q)/gcd` — *which* `m` first realigns the `+p` and `−q` characters — not a
  finite collapse of the functional. This is the same distinction that makes the GMC detection
  depth grow (THM-1770/1790): the moment is a full radial/angular integral, not a finite sum.

## Scope

A frame, not a bound: every equation is verified-exact or classical. Its value is the **routing**
— it gives the owner's a/b/triangular/trig archeology a single home and connects it, through the
constant-term = circle-average functional, to the two-charge base case that the GMC(2)
formalization proves premise-free in Lean (`GMC2DvdKTwoCharge`).

**Currency note (2026-07-22).** GMC(2) itself was closed the same day this file was written:
`GMC2Main.GMC2.gmc2` is unconditional and kernel-pure (kind-pasteur's `singlePolyCrux_holds`
Ω-wiring, HYP-9020; corroborated by mac-mini-S167, boxeph-S243). Nothing here depends on GMC(2)
being open — the two-charge case cited above is and remains an elementary premise-free lemma, and it
is the *base case* of the moment engine, not the hard content. The forward multiple-angle
direction is THM-1575/1880; the backward dyadic angle-bisection (the `b`-iteration on angles) and
the `O_n`-as-tournament-invariant / tangent-number `W` revival (THM-293) remain the named,
un-pursued next steps.

*Files: `04-computation/ab_halfangle_ct_bridge_klein_S401.py` (+ `.out`).*
