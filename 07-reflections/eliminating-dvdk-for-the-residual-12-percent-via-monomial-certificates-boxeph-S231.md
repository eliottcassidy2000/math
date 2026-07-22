# Eliminating DvdK for the residual 12%: the monomial certificate

*boxeph-2026-07-22-S231. Owner: get rid of DvdK for the remaining ~12% of straddling supports.
Builds on S229/S230 (the unique-channel bypass and `dvdk1_of_uniqueChannel`), death-star-S101/HYP-8878
(unique minimal channel, the 84%), codex THM-2067 (the general DvdK, Galois orbit-product) and
THM-2070 (the dihedral/involution obstruction). New: `04-computation/dvdk_monomial_certificate_
residual_boxeph_S231.py` and Lean `GMC2DvdKMonomialCertificate.lean` + `GMC2DvdKResidualExample.lean`,
both kernel-pure.*

## Where the "12%" comes from, exactly

death-star-S101 measured 84.5% (98/116) of straddling supports (size 3–4 in `[-4,4]`) DvdK-free by a
unique *minimal* balanced channel. My S230 `dvdk1_of_uniqueChannel` needs a unique channel at *some*
mass, not the minimal one — re-scanning at any mass ≤ 40 reclassifies **4** more supports as free
(e.g. `(-4,-3,1,2)`), giving **87.9% (102/116)**. The residual is then exactly **14/116 = 12.1%** —
the "12%". These are the coincident-channel / symmetric supports where the involution `u ↦ -1/u`
(`f(-1/u) = -f(u)`, THM-2070) pairs balanced compositions, forcing `card ≥ 2` at *every* mass, so no
single power can certify non-vanishing coefficient-independently.

## The tool: a monomial certificate — DvdK without Galois

`CT(f^m) = constantTermRelation(q,m)` is homogeneous of degree `m` in the coefficient variables, with
coefficient `multinomial(m;r)` on `x^r` for each balanced composition `r`. The key reframing:

> a support is DvdK-free **elementarily** iff the ideal `I = ⟨CT(f^m) : m ≥ 1⟩` contains a **monomial**
> `μ = ∏ xᵢ^{eᵢ}`. Then on the coefficient torus `μ(c) = ∏ cᵢ^{eᵢ} ≠ 0`, so from `μ = ∑ₘ gₘ·CT(f^m)`
> some `CT(f^m)(c) ≠ 0` — for the *same finite* mass set, for every torus point.

"`I` contains a monomial" is `V(I) ∩ torus = ∅`, which *is* DvdK1 for the support (true by THM-2067) —
but here produced as an **explicit finite, formalizable certificate** per support, replacing the
Galois proof by exact rational linear algebra. Because each `CT(f^m)` is homogeneous, "`I` contains a
degree-`D` monomial" is a graded linear-algebra question: is some degree-`D` monomial in the ℚ-span of
`{ t·CT(f^m) : m ≤ D, deg t = D-m }`? I answered it by exact Gaussian elimination over `Fraction`
(no numpy/sympy).

**Result (`dvdk_monomial_certificate_residual_boxeph_S231.py`): all 14 residual supports carry a
certificate at degree ≤ 6.** So every straddling support of size 3–4 in `[-4,4]` is now DvdK-free
elementarily — 88% by a unique channel, the remaining 12% by a monomial certificate — with **no
Galois orbit-product used for any of them**.

The paradigm `{-2,-1,1,2}` (charges `-2,-1,1,2`, coefficients `a,b,c,d`):
```
   CT(f²) = 2(ad + bc),   CT(f⁴) = 6a²d² + 24abcd + 6b²c²
   12·b²c² = (3ad + 9bc)·CT(f²) − CT(f⁴).
```
A two-mass certificate `{2,4}`: on the torus (`bc ≠ 0`), `CT(f²)=0 ⟹ CT(f⁴) = -12 b²c² ≠ 0`.

## Formalized, kernel-pure

Two Lean modules (`#print axioms = [propext, Classical.choice, Quot.sound]`):

- **`GMC2DvdKMonomialCertificate.dvdk1_of_monomialCertificate`** — the reusable engine. Given the
  polynomial identity `∏ Xᵢ^{eᵢ} = ∑_{m∈M} gₘ · constantTermRelation(q,m)` and `cᵢ ≠ 0`, it returns
  `∃ m ≥ 1, CT(f^m) ≠ 0` (over any field ℚ-algebra). No positivity, no uniqueness, no Galois. The
  unique-channel case is the special case `M = {m₀}`, `μ` the single balanced monomial; this strictly
  generalizes it to the coincident-channel stratum.
- **`GMC2DvdKResidualExample.dvdk1_neg2_neg1_1_2`** — DvdK1 fully discharged for the paradigm residual
  `{-2,-1,1,2}`. I evaluate `CT(f²)` and `CT(f⁴)` at a generic `c` (isolating the balanced subsets by
  `decide` on the `piAntidiag` filter), then the certificate `12b²c² = (3ad+9bc)CT(f²) − CT(f⁴)`
  (a `linear_combination` step) gives the two-mass disjunction. Kernel-pure, no external axiom.

## Scope and honesty

Concrete: the residual 12% is eliminated. Computationally, **all 14** residual straddling supports of
size 3–4 in `[-4,4]` carry a finite monomial certificate (degree ≤ 6), so DvdK holds for each with no
Galois. In Lean, the general certificate engine is proved kernel-pure, and one representative residual
support (`{-2,-1,1,2}`, the symmetric paradigm) is fully discharged kernel-pure. The other 13
certificates are script-verified; mechanizing each is the same `decide`-the-`piAntidiag` +
`linear_combination` recipe (mechanical, not yet all in Lean). What is *not* proved is a **uniform**
bound — that every straddling support (unbounded size/range) has a bounded-degree certificate; that is
exactly effective DvdK (Sturmfels/ESV open, THM-2067 §), and remains the residue. But for any explicit
support the certificate is a finite computation, so DvdK is dischargeable support-by-support with no
Galois orbit-product — which is the concrete sense in which the 12% (and, by the same method, any
explicit case) is done.

Links: HYP-8932, HYP-8931, HYP-8878, THM-2067, THM-2070,
[[bypassing-the-gmc2-dvdk-dependency-for-the-unique-channel-class-boxeph-S230]],
[[the-unique-channel-dvdk-in-lean-and-the-cancellation-inclusion-exclusion-dictionary-boxeph-S229]].
