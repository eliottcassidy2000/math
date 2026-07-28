---
id: THM-2647
title: "Endpoint-anchored two-point deconvolution and the thirteen-halves signed tax"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.  On an odd
  cycle, the two-point kernel I+T_d has the explicit inverse
  (1/2)sum_j(-1)^j T_(jd).  For C_13 this is (I-H_d)/2, the inverse already
  latent in THM-2368 and THM-2526; it is unique, necessarily signed, and has
  l1 norm 13/2.  Hence the full THM-2645 multiplicity profile together with
  one absolute labelled endpoint two-set A recovers the other endpoint B
  exactly.  No endpoint leaves the thirteen-fold origin gauge, while one
  marked member leaves two representatives per gauge orbit.  The theorem
  supplies a sharp minimal algebraic sidecar, not the missing physical
  endpoint reference or an LRC row exclusion.
source: wild-holotopy-mining-2026-07-28-endpoint-deconvolution
depends_on:
  - THM-2645-eleven-sheet-multiplicity-full-character-spectrum-and-energy-split
related:
  - THM-2366-retained-probe-target-covariance-and-subthirteen-budget-bridge
  - THM-2368-owner-pivot-root-fibre-radon-invertibility
  - THM-2521-k13-drift-k14-potential-module-bridge
  - THM-2526-affine-skew-orientation-gauge-boundary
  - THM-2634-endpoint-pair-two-carry-cospan-and-single-carry-no-go
script: 04-computation/lrc14_endpoint_anchored_two_point_deconvolution_thm2647.py
output: 05-knowledge/results/lrc14_endpoint_anchored_two_point_deconvolution_thm2647.out
script_sha256: b2a8b9346141685e023629c3a651a1f1cc8e5daa24798b925eb7e572a1b50a3b
output_sha256: f3a6e66efd2b539e0b53788409424aa83e023a90224ce6314a3185e3c3f9aab7
hash_basis: LF-normalized bytes
---

# THM-2647 -- one endpoint deconvolves the saturated relation pair

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT AUDIT.**

THM-2645 proves that exact multiplicity retains every carry character but
forgets how the two endpoint relations share a common origin.  An old local
inverse from THM-2368 and the skew Cayley transform of THM-2526 together give
the exact missing algebraic consumer: one absolute endpoint two-set is enough
to reconstruct the other.  The inverse is necessarily signed, and its total
variation is exactly `13/2`.

## 1. Odd-cycle two-point inversion

Let `d` be an element of odd order `q` in a finite group, and work in its
rational group algebra.  Put

```text
K_d=delta_e+delta_d,
V_d=(1/2)sum_(j=0)^(q-1)(-1)^j delta_(d^j).              (1)
```

The finite geometric series gives

```text
K_d V_d
 =(1/2)(delta_e+delta_d)
          (delta_e-delta_d+...+delta_(d^(q-1)))
 =(1/2)(delta_e+delta_(d^q))
 =delta_e.                                                (2)
```

Thus `K_d` is a unit and `V_d` is its unique inverse.  Since the `q` powers
of `d` are distinct,

```text
||V_d||_1=q/2,          sum_g V_d(g)=1/2.                 (3)
```

For `q>1`, exactly `(q-1)/2` coefficients are negative.  No nonnegative
translation-equivariant inverse exists because any such inverse would equal
the unique `V_d`.

Odd order is sharp.  If `d` has even order, the alternating vector on a
`d`-orbit is killed by `K_d`; equivalently `-1` is an eigenvalue of the
translation `T_d`.

## 2. The THM-2368/2526 identity is this inverse

On `C_13`, write translations additively and define

```text
H_d=sum_(j=1)^12 (-1)^(j+1) T_(jd).                      (4)
```

THM-2526 proves the telescoping Cayley identity

```text
(I+T_d)H_d=T_d-I.                                        (5)
```

Rearranging (5) gives

```text
(I+T_d)(I-H_d)=2I,
(I+T_d)^(-1)=(I-H_d)/2=V_d.                              (6)
```

This is exactly the one/two-hole inverse `V` printed earlier in THM-2368.
The present theorem does not invent a new local inverse; it identifies its
previously unused endpoint-allocation consequence.

The circulant determinant is also exact:

```text
det(I+T_d)=prod_(k=0)^12(1+zeta_13^(kd))=2.               (7)
```

Hence the half-integral coefficients in (6) are not a normalization
artifact.  They record the index-two integral image of the two-point
convolution operator.

## 3. One anchored endpoint recovers the other

Let `p` be an odd prime, let `A,B` be two-point subsets of `C_p`, and put

```text
S=C_p\A,                  T=C_p\B,
r=1_S*1_T,                m=1_A*1_B.                     (8)
```

Inclusion--exclusion, as in THM-2645, gives

```text
m=r-(p-4)1.                                               (9)
```

Suppose the **left** endpoint `A` is supplied in one absolute common-origin
gauge.  Choose either ordering

```text
A={a_0,a_1},             d=a_1-a_0.                     (10)
```

Since

```text
1_A=T_(a_0)(delta_0+delta_d),                             (11)
```

equations (1)--(2) recover the other endpoint by

```text
1_B=V_d*T_(-a_0)m
   =((I-H_d)/2)T_(-a_0)[r-(p-4)1]             (p=13).     (12)
```

The answer is independent of the ordering chosen in (10): both formulas are
inverses of convolution by the same mask `1_A`.  More strongly, among all
rational functions on `C_p`, there is exactly one `B` solving (8) once `A`
and `m` are fixed.  Thus one labelled endpoint two-set, not a complete
endpoint pair, kills the allocation ambiguity.

For normalized Fourier transforms the same statement is

```text
Bhat(k)=mhat(k)/(p Ahat(k)),                               (13)
```

and every denominator is nonzero because an odd-order two-term root sum
cannot vanish.  Formula (12) is stronger for the present purpose: it exposes
the exact rational signs and the sharp norm invoice.

## 4. The common-origin sidecar ladder is sharp

The THM-2645 gauge action is

```text
(A,B) -> (A+u,B-u),                u in C_p.              (14)
```

It preserves `m` and `r`.  For a proper two-point `A` on `C_p`, the action is
free, so no absolute endpoint information leaves `p` representatives in each
gauge orbit.

Fixing only one absolute point `x` and requiring `x in A+u` leaves exactly
two representatives **within that orbit**:

```text
u=x-a_0                 or                 u=x-a_1.       (15)
```

They are distinct because `p` is odd and `a_0!=a_1`.  Supplying the full
absolute set `A` leaves one representative, and then (12) determines `B`.
This proves the natural ladder

```text
no endpoint anchor: p;
one marked member: 2 per common-origin orbit;
one full labelled endpoint two-set: 1.                   (16)
```

Equation (16) does not claim that the common-origin orbit is the only
possible factorization ambiguity when no endpoint is given.  It says that
the full `A` sidecar removes both that orbit and every remaining ambiguity in
`B`, by invertibility.

If the left/right endpoint role is forgotten, convolution is commutative and
the chronology swap remains.  The word **labelled** in (16) is load-bearing.

## 5. The `13/2` signed tax

For `p=13`, (3) says

```text
V_d has seven +1/2 coefficients and six -1/2 coefficients,
||V_d||_1=13/2,       ||(V_d)_+||_1=7/2,
||(V_d)_-||_1=3.                                         (17)
```

This gives a structural explanation for a recurring constant.  THM-2366
shows that a hypothetical **positive** role-to-target intertwiner must pay an
amplification of at least `13/2`; THM-2521/2526 find the same half-cycle scale
in an aligned Radon/Cayley operator.  Here `13/2` is the exact total variation
of the unique signed inverse that solves endpoint deconvolution.

These are comparisons, not identifications.  The three theorems act on
different packets and no canonical map equates their gauges.  What THM-2647
does prove is the local obstruction underneath the numerical coincidence:
odd-cycle two-point inversion sits exactly at a signed `13/2` boundary, so a
mass-preserving positive inverse cannot exist.

## 6. Hostile controls and physical boundary

Every retained datum matters.

1. **Even cycle:** `I+T` kills the alternating vector, so no inverse exists.
2. **No endpoint:** all thirteen transforms (14) have the same `r`.
3. **One point only:** the two choices in (15) remain distinct.
4. **Energy only:** with `A={0,1}`, the endpoints `B={0,2}` and `B={0,3}`
   have different multiplicity profiles but the same THM-2645 charged energy
   `36/169`.
5. **Positivity:** the unique inverse (17) has six negative coefficients.
6. **Unlabelled chronology:** swapping `A,B` preserves `m`.

Most importantly, current LRC canon supplies neither input of (12) on the
needed physical object.  THM-2645 is conditional on a same-base positive
two-edge transition table; the current eleven-sheet rows are static
coefficient fibres.  THM-2368 supplies local rooted inverses but explicitly
loses their phase under integration.  No theorem supplies an absolute
left-endpoint two-hole mask in the same carry gauge as a physical THM-2645
multiplicity profile.

Thus (12) is an exact minimal-sidecar consumer, not a construction of its
sidecar.  It does not produce a canonical THM-2334 current, endpoint phase,
semantic owner, terminal word, or scalar row exclusion.  The LRC ledger
remains `165`.

## 7. Exact reproduction

Run

```bash
python 04-computation/lrc14_endpoint_anchored_two_point_deconvolution_thm2647.py
python -O 04-computation/lrc14_endpoint_anchored_two_point_deconvolution_thm2647.py
```

The dependency-free exact companion uses `Fraction` arithmetic and explicit
optimization-safe guards.  It checks every nonzero translation on odd cycles
of sizes `3,5,7,9,11,13`, the even-cycle alternating hostile through size
`12`, all twelve `C_13` Cayley identities and determinant/index controls, and
all `6,084` ordered endpoint pairs.  For each pair it reconstructs `B` in
both orientations of `A`, proves uniqueness against all `78` candidate
two-sets, and checks the `13 -> 2 -> 1` gauge ladder.  Normal and optimized
executions must byte-match the stored transcript and end in `PASS`.

QED (candidate; independent audit pending).
