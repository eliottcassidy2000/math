# The Vitali wall + Brouwer's equioscillation: measure handles the bulk, the CYCLOTOMIC construction handles the core — and my recent cyclotomic work IS the construction side

*mac-mini-2026-06-28-S75f. The owner: see things as functions/information; let Brouwer's fixed-point theorem and
Vitali inspire; explore the unknown for proof gems. The Vitali wall (oracle S551o) + Brouwer's max-min
equioscillation + my cyclotomic work (S75d/e) lock together: the proof splits BULK(measure) ⊕ CORE(construction),
the core IS the cyclotomic primitive-root witnesses, and Brouwer locates them at the equioscillation. Builds on
[[the-vitali-wall-measure-vs-set-and-why-the-core-needs-construction-not-measure-s551o]],
[[the-cap-is-a-totally-real-cyclotomic-quantity-and-n14s-two-cyclotomic-heads]], HYP-2840.*

## The function / information lens (the frame the owner asked for)
The LRC is the **safety function** `S(t) = min_i ‖s_i t‖` on the circle; `M(S) = max_t S(t)`, and LRC(14) is
`M(S) ≥ 1/14`. Two ways to read it:
- **as MEASURE (circle-method):** control `meas{t : S(t) > 1/14}` — an *information-rich* tool on the bulk, but it
  carries **zero information at a measure-zero maximum** (it cannot tell empty from measure-zero-nonempty).
- **as a SET / CONSTRUCTION:** exhibit a `t` with `S(t) ≥ 1/14` — an arithmetic witness.
The Vitali set is the sharp diagnosis of when measure has no information: a full transversal of `ℝ/ℚ` is nonempty
yet has no measure. **Measure cannot certify set-nonemptiness.**

## The Vitali wall (verified) = the BULK/CORE split
`lrc_vitali_core_brouwer_macmini_S75.py` (this) + oracle S551o:
- **BULK** (generic/spread sets): `meas{S(t)>1/14} = 0.136 > 0` — the circle-method/equidistribution proves
  nonemptiness. (My S74 analytic resolution lives here.)
- **CORE** (the AP / n-gon, the tight locus): `meas{S(t)>1/14} = 0` — the open lonely set is EMPTY, measure is
  **blind**. Yet the closed witnesses are NONEMPTY and **explicitly cyclotomic**:
```
  AP {1..13}, n=14:  closed witnesses t=a/14 with S>=1/14  ⟺  a ∈ {1,3,5,9,11,13}
                     = the UNITS mod 14 = the PRIMITIVE 14th roots = Φ₁₄  (φ(14)=6 of them)
```
So **the core witnesses are exactly the primitive 14th roots `Φ₁₄`** — the very factor I flagged in S75d as the
"half-tiling 2·7 mixing." The Vitali wall says: *measure ends at the core; the cyclotomic construction begins.*

## Brouwer / EVT: the core is the EQUIOSCILLATION (the max-min saddle)
`S(t)` is continuous on the compact circle, so `max_t S(t)` is **attained** (the 1-D Brouwer / extreme-value
theorem) — existence is free; the *value* (`≥1/14`) is the content. At the core the max is attained on the
cyclotomic `Φ₁₄` points with an **equioscillation**:
```
  at t=1/14:  runners ±1 (i.e. 1 and 13) sit at EXACTLY 1/14  — the max-min SADDLE.
```
This equioscillation is the **Brouwer/Kakutani saddle of the observer-vs-runners game** and is exactly kps's
de Moivre/Chebyshev equioscillation (`Φ₇`, S31ao/an): the max-min value is pinned by the runners that
simultaneously touch `1/14`. **`Φ₁₄` = where the core lives; `Φ₇` (de Moivre) = how it equioscillates.**

## THE GEM: my cyclotomic work IS the Vitali-construction side
The Vitali wall says the proof must be **BULK (measure) ⊕ CORE (construction)**, and the core construction must be
arithmetic (cyclotomic), not measure. My recent sessions built exactly that construction:
- **S75e:** the cap is a totally-real cyclotomic quantity in `Q(cos 2π/7)`; the magic function `F₇ =
  (de Moivre cubic)²` is the totally-positive square. → the **`Φ₇` equioscillation** side.
- **S75d:** the cap char poly `= (x−1)^depth · Φ₂ · Φ₇ · Φ₁₄`; `Φ₁₄` = the primitive-root mixing. → the **`Φ₁₄`
  core-witness** set (now identified as the units mod 14 = the closed witnesses).
- **S75b:** the three-gap kernel = the Diophantine/Stern-Brocot recursion. → the **construction algorithm** for the
  witnesses (the continued-fraction / largest-gap witness).
So the **proof skeleton** the Vitali wall demands is:
> **LRC(14) = [BULK: `meas{S>1/14}>0` by circle-method/equidistribution (S74)] ⊕ [CORE: the cyclotomic `Φ₁₄`
> witnesses `t=a/14` are nonempty & safe, with the `Φ₇` de Moivre equioscillation (S75d/e), Vitali marking the
> boundary, Brouwer attaining the saddle].**
My cyclotomic program is not a side-structure — it is the **measure-blind core half** of the proof.

## Exploring the unknown (gems to chase)
- **Brouwer degree / index:** the equioscillation saddle has a topological index; the number of core witnesses
  is `φ(14)=6` (the `Φ₁₄` units). A degree/index count of the saddle could *force* a witness (a topological
  lower bound on `M`), the construction analogue of the measure bound.
- **The information handoff is sharp:** measure-information `→ 0` exactly as the cyclotomic-arithmetic information
  turns on (the conductor `7²` of S75e appears only in the binding/core rows). Information is *conserved across
  the Vitali wall* — measure on the bulk, arithmetic on the core.
- **Besicovitch/Vitali covering (HYP-2840):** the bulk bound itself is a Vitali-covering of the circle by the
  danger combs; the overlap kernel `K(a,b)` (S75b) is the covering multiplicity. The two Vitali uses (covering
  for the bulk, the wall for the core) are the two halves again.

## Honest status
- **VERIFIED:** the AP core open-lonely measure `=0` (Vitali wall); the closed witnesses `= Φ₁₄` units mod 14;
  the equioscillation `±1` at `1/14` (Brouwer saddle); the bulk has positive measure.
- **SYNTHESIS/GEM:** the Vitali wall splits LRC(14) into measure-bulk ⊕ cyclotomic-construction-core; **my S75d/e
  cyclotomic work is the construction side** (Φ₁₄ witnesses, Φ₇ equioscillation); Brouwer/EVT attains the saddle;
  the function/information lens shows information is conserved across the wall (measure→arithmetic).
- **NOT a proof.** Both halves retain open pieces (the bulk's effective equidistribution; the core's full
  witness construction across all covering sets, not just the AP). But the proof's *architecture* is now pinned to
  the Vitali wall, and my cyclotomic work is correctly placed as the core-construction half. LRC(14) open.

Related: HYP-3237 (this), S551o Vitali wall, HYP-2840 (Vitali covering), HYP-3235 (cap in Q(cos2π/7)), HYP-3233
(cyclotomic factors, Φ₁₄), HYP-3221 (the one obstruction), OPEN-Q-108.
