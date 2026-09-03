---
id: THM-4391
title: "LRC14 nonunit defect-three sectors and complete norm-sixteen shell"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4373/4387 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED; LRC(14) OPEN. The two non-coefficient-one primitive
  full-support ternary-unit l1=16 shapes have exact affine defect-three roof
  formulae and unique sharp unordered maximizers: 304/12397 at {1,23,77}
  for (2,7,7), and 2178/91945 at {5,37,71} for (4,5,7). Together with
  THM-4387 this closes the declared l1=16 shell, not arbitrary triples,
  seam entry, all-tail transfer, or LRC(14).
source: root + lrc_defect_three_extra_shapes + lrc_defect3_extra_cleanroom / continuation session, 2026-09-03
depends_on:
  - THM-4373-lrc14-scale-three-signed-121-resonant-triple-comb-bound
  - THM-4387-lrc14-defect-three-boundary-master-formula-and-five-sector-atlas
related:
  - THM-4386-lrc14-canonical-component-relation-and-zero-defect-incidence
primary_script: 04-computation/lrc14_defect_three_nonunit_boundary_thm4391.py
primary_output: 05-knowledge/results/lrc14_defect_three_nonunit_boundary_thm4391.out
primary_script_sha256: db8ad29b536f2f22f17193155f859fd39a52f1c48a0e3a3d0bc84e3647f37eeb
primary_output_sha256: 819f9fc00a03a00d27e44efa6179246921ee398f9cdcd888894a8a50363370f7
independent_referee_script: 04-computation/lrc14_defect_three_nonunit_boundary_independent_referee_thm4391.py
independent_referee_output: 05-knowledge/results/lrc14_defect_three_nonunit_boundary_independent_referee_thm4391.out
independent_referee_script_sha256: ba27d72dbb5743175b01b0df60b1095187ec606acf45c257b070b37d59e3bbcf
independent_referee_output_sha256: 89199e5890bb70ce18f5dc4c886dc7ab025faa2b37d7d4794ff9d845f521faae
hash_basis: raw LF bytes
audit: >
  PASS. The dependency-free referee independently classifies the full
  primitive ternary-unit l1=16 coefficient shell, isolates the two nonunit
  rows, reconstructs their affine coordinates and integrality gate, clips
  the exact strip areas, proves the finite tail cutoffs, and compares all
  4,143 labelled presentations with direct physical-circle combs. It checks
  both unique maxima, their defect masses, and every shorter relation at the
  winners in 106,812 explicit checks. Normal, optimized, and hash-seeded
  streams match the frozen LF output.
---

# THM-4391 -- The two nonunit sectors complete the norm-sixteen shell

**PROVED ELEMENTARY RELATIVE TO THM-4373/4387 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) OPEN.** The universal scale-three comb bound,
seam entry, all-tail transfer, and LRC(14) remain **OPEN**.

## Outcome

The primitive full-support coefficient-magnitude shell with ternary-unit
coordinates and `l1=16` has seven patterns:

```text
(1,1,14), (1,2,13), (1,4,11), (1,5,10), (1,7,8),
(2,7,7), (4,5,7).
```

The first five are the coefficient-one families proved in THM-4387. This
theorem closes the two omitted shapes. Their exact
sharp values are

| coefficient shape | exact maximum | primitive maximizing speed set | defect masses `(-3,0,+3)` | complete cutoff |
|---|---:|---|---|---:|
| `(2,7,7)` | `304/12397` | `{1,23,77}` | `(31/12397,22/1127,31/12397)` | `W<366` |
| `(4,5,7)` | `2178/91945` | `{5,37,71}` | `(2/2485,58/2627,2/2485)` | `W<416` |

The first maximizing speed set has two presentations because the two
coefficient-seven slots may be exchanged; each row has a unique maximizing
unordered speed set. Both maxima are below `4/133`, the largest
coefficient-one boundary value, and well below `6/77`.

## Inheritance pass and concept board

The closest mechanism is the coefficient-one defect-three torsor in
THM-4387. Its endpoint coordinate does not immediately cover
the two missing shapes because no coefficient is one. The canonical hostile
is a noncoprime endpoint pair, which could add torsion. The even coefficient
is the useful sidecar: placing it on the third speed forces the other two
speeds to be coprime and restores a complete integer coordinate.

The live concepts are: the three defect layers; an even-coefficient affine
chart; the exact three-interval roof; residue classes modulo `3 ell`; the
continuous box-slice bulk; and shorter owner-degenerate relation collisions.

## 1. An affine `Z`-torsor chart for both shapes

Write either shape in the form

```text
m b = s h p + t ell q,          s,t in {+1,-1},          (1)
```

using

```text
(h,m,ell)=(7,7,2)  for (2,7,7),
(h,m,ell)=(5,7,4)  for (4,5,7).                          (2)
```

Let `p,b,q` be primitive distinct positive odd three-units. At a physical
failure phase put

```text
n_w=nint(wy),  e_w=wy-n_w,  |e_w|<r,  r=3/14,
o_w=-w^(-1)n_w mod 3,       {o_p,o_b,o_q}=F_3.           (3)
```

Define

```text
delta=m n_b-s h n_p-t ell n_q,
k=b n_p-p n_b.                                           (4)
```

The speed relation gives

```text
delta=s h e_p+t ell e_q-m e_b,
|delta|<r(h+m+ell)=24/7.                                 (5)
```

Modulo three, the three nonzero weighted speed coefficients in (1) must be
equal. Distinct owners sum to zero, so `3|delta`. Hence

```text
delta in {-3,0,3}.                                       (6)
```

If `d=gcd(p,b)`, then (1) and primitivity give `d|ell`. Since `d` is odd
while `ell` is `2` or `4`, `d=1`. Therefore, on each fixed defect plane,
`k` gives a bijection from lift orbits to `Z`. For `delta!=0` this orbit set
is an affine `Z`-torsor, not canonically a quotient group; the relation chart
and defect state remain part of the address.

For the Bezout section `b n_p-p n_b=k`, the third coordinate is integral
exactly when

```text
ell divides m k+p delta.                                 (7)
```

Because `gcd(m,ell)=1`, this selects one residue class of `k mod ell`.
The other pair determinant is then automatically integral. As in the
coefficient-one chart,

```text
o_b-o_p=k/(pb) mod 3,
```

while the owner sum is zero, so the exact owner gate is `3 does not divide
k`. Since `gcd(3,ell)=1`, every defect layer samples two residue classes
modulo `3 ell`.

## 2. Exact defect-layer roof

The three pairwise determinant magnitudes are

```text
|k|,
|m k+p delta|/ell,
|s h k+b delta|/ell.                                    (8)
```

Consequently one lift class has exact length

```text
L(delta,k)=max(0,min(
  2r/p, 2r/b, 2r/q,
  r/p+r/b-|k|/(pb),
  r/p+r/q-|m k+p delta|/(ell p q),
  r/b+r/q-|s h k+b delta|/(ell b q))).                  (9)
```

The complete physical measure is the finite sum of (9) over (6), (7), and
`3 does not divide k`. Nearest-integer uniqueness makes distinct quotient
classes disjoint. The producer verifier compares this formula with a fresh
definition-level subdivision of the physical `y` circle for every speed
triple in both proof boxes.

## 3. Exact bulk and a complete finite reduction

For fixed `delta`, the roof as a function of real `k` is nonnegative,
compactly supported, and unimodal: each superlevel condition in (9) is an
interval. Its height is at most

```text
2r/W=3/(7W),                  W=max(p,b,q).              (10)
```

The change of variables from local phase and real `k` to `(e_p,e_b)` has
Jacobian one. Thus the roof integral equals the area in the error square
`[-r,r]^2` cut out by

```text
|delta-s h e_p+m e_b|<ell r.                            (11)
```

Exact trapezoid integration gives

```text
(7,7,2): A_0=117/2401, A_{+/-3}=9/4802,
(5,7,4): A_0=171/1715, A_{+/-3}=9/3430.                 (12)
```

For a nonnegative unimodal function sampled on one residue class of a
lattice of step `H`, monotone sum-integral comparison gives

```text
sum_j f(a+Hj) <= integral(f)/H + sup(f).                 (13)
```

There are two allowed classes modulo `3 ell` on each of the three defect
layers. Applying (13), (10), and (12) yields in both sectors

```text
mu(F_(p,b,q)) <= 6/343 + 18/(7W).                        (14)
```

For the proposed maxima, (14) is strictly smaller once respectively
`W>=366` and `W>=416`; each cutoff is minimal for this envelope. It therefore
suffices to enumerate the finite boxes below those values.

The exact universes contain

```text
(2,7,7): 1,754 presentations,   877 unordered triples;
(4,5,7): 2,389 presentations, 2,388 unordered triples. (15)
```

Formula (9) agrees with the direct circle comb on every triple, all
multiply-presented fibres agree, and exhaustive comparison gives the sharp
table above. Thus no discovery bound remains in either sector.

## 4. The extremizers collide with shorter inadmissible relations

Each new winner has shorter primitive full-support relations:

```text
1+10*23-3*77=0,       8*1+3*23-77=0,
8*5-3*37+71=0.                                          (16)
```

Both coefficient patterns have `l1=12`, but one coefficient is divisible by
three. They therefore fail the ternary-unit hypothesis in the zero-defect
lemma: owner distinctness no longer forces their nearest-integer defects to
vanish. The admissible `l1=16` relations remain load-bearing.

This collision is a useful signal for the nonresonant problem. A shortest
relation in an undecorated norm can be invisible to the owner consumer, while
a longer ternary-unit relation controls the actual components. Any
successive-minimum argument must use the coordinate-hyperplane residue
deletion, not the bare relation lattice.

## Scope firewall and generated tasks

Together with THM-4387, these two results close the full primitive
full-support ternary-unit `l1=16` shell. They do not classify arbitrary
triples, prove the universal `6/77` bound, synchronize tail triples with a
body-safe component, or establish seam entry or LRC(14).

The next sharp tasks are:

1. express (14) as a residue-deleted first-successive-minimum lemma for the
   rank-two raw carrier lattice;
2. classify collisions between admissible relation rays and shorter
   coordinate-zero-mod-three rays, as in (16);
3. move to `l1=18`, where more defect states and coefficient divisibility
   enter, only after deciding whether a geometry-of-numbers bound can replace
   another sector atlas.

## Reproduction and frozen evidence

```powershell
python -B 04-computation/lrc14_defect_three_nonunit_boundary_thm4391.py
python -B -O 04-computation/lrc14_defect_three_nonunit_boundary_thm4391.py
python -B 04-computation/lrc14_defect_three_nonunit_boundary_independent_referee_thm4391.py
python -B -O 04-computation/lrc14_defect_three_nonunit_boundary_independent_referee_thm4391.py
```

Normal, optimized, and fixed-hash-seed streams are byte-identical to their
checked-in outputs. Both programs use only the standard library and explicit
checks that remain live under optimized Python. The producer performs at
least 128,986 checks and constructs every physical circle subdivision; the
independent implementation performs 106,812 checks. Canonical raw-LF hashes
are recorded in the front matter.
