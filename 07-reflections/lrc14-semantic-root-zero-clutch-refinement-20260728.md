# The semantic root-zero clutch survives, but its coefficients shear

> **STATUS: FINITE-EXACT + VERIFIED; AWAITING INDEPENDENT HOSTILE AUDIT.**
> On the canonical rail-8 adjacent-chart overlap, insert the complete `E3`
> source mask, the delayed `Q_(3,{1,2})` word, and the common lawful target
> sheet `U_(0,3)`.  Independently recomputed source-root-12 and target-root-1
> seven-clock coefficients are both primitive units.  They are not equal.
> This strengthens the physical-survival side of the clutch while sharply
> locating the remaining debt in coefficient/phase transport.  It is not an
> endpoint current, row exclusion, or proof of LRC(14).

## 1. Inheritance and the precise question

The exact root-zero overlap sidecar proves that the scalar lift from semantic
residue `7` to residue `8` lands in the common physical strip

```text
H^R_12 intersect T^(-1)H^L_1=(169,181)/182
  <-> (1,13)/182=H^R_0 intersect H^L_1.                 (1)
```

Before adding the semantic and target masks, the two recomputed overlap
vectors agree and are units.  Pointwise, the canonical source and target both
have semantic record `E3 -> D^6 -> Q_(3,{1,2})`, and `(s,t)=(0,3)` belongs to
both lawful target-label banks.  Pointwise coexistence does not imply that the
globally integrated, further-cut coefficients still agree or even survive.

The probe therefore inserts all three requirements *inside* the integral.
On the physical coordinate it intersects with the full `E3` mask and the four
shifted safe factors defining `U_(0,3)`.  On the delayed coordinate it rebuilds
the sector-zero prefix with both lower blockers dangerous and the source
blocker safe, exactly `Q_(3,{1,2})`.  Source and translated-target rail weights
are retained separately, and both delayed-carry vectors are recomputed.

## 2. Exact outcome

The seven source and target weighted masses are constant in the shallow-clock
index but differ across the chart transport:

```text
source mass = 929934280541992017600,
target mass = 929934304688494607040.                    (2)
```

The exact coefficient vectors are

```text
A_12=(0,a,a,a,a,a,a),   a=1812281403506324508838080,
B_1 =(0,b,b,b,b,b,b),   b=1826551436254490256030720.    (3)
```

Both have canonical content `26`.  After the THM-2640 root normalization and
reduction modulo `Phi_7`, their profiles are

```text
A_12 -> (5,0,0,0,0,0),
B_1  -> (9,0,0,0,0,0).                                 (4)
```

Thus both multiplication determinants are nonzero: the complete semantic
source, delayed word, target sheet, and partial physical rechart coexist with
private units at both nonzero roots.  However `a!=b`; the full mask is not
translation invariant and destroys the uncut coefficient intertwiner.

This is a useful positive/hostile pair.  The clutch no longer lacks a lawful
coefficient on either side.  What it lacks is a canonical map identifying or
phase-pairing the two unequal coefficients while retaining the endpoint
current.

## 3. Boundary and next test

The result uses one common sheet on one canonical semantic rail.  It does not
show that all common sheets survive, that the two coefficients have a lawful
ratio in a common residue field, or that either coefficient is the marked
endpoint current of the global scalar cover.  Equality of the uncut overlap
vector cannot be transported through an arbitrary semantic mask.

The cheapest next test is a target-character bank over all common `(s,t)`
labels, retaining the two amplitudes separately.  The relevant object is not
another support intersection but the cross coefficient

```text
sum_(s,t) chi(s,t) A_12(s,t) conjugate(B_1(s,t)),        (5)
```

with the inherited delayed carry and endpoint typing kept fixed.  A nonzero
value would supply the missing phase-pairing sidecar; vanishing would expose a
new finite character obstruction.

## 4. Exact reproduction

Run

```bash
python 04-computation/lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py
python -O 04-computation/lrc14_semantic_root_zero_clutch_refinement_probe_20260728.py
```

Both modes byte-match
`05-knowledge/results/lrc14_semantic_root_zero_clutch_refinement_probe_20260728.out`.
The script uses exact rational/integer interval arithmetic and explicit
exceptions, with no truth-bearing Python `assert`.

LF-normalized SHA-256:

```text
script  cf8f0757f5097405c2ec86688624366fb4aea9ae010cd37fe8f80bb358070056
output  ff76e371820d3b7721804f4abcfaf9d16114a5d3c404d8b2ab02a1f407e776e4
```
