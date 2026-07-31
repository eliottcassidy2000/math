---
id: THM-2767
title: "Split degrees ten, fourteen, and eighteen componentwise monicization closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every planar
  Keller pair in the complete chosen-sheet split polynomial exact-square-
  prefix reduced-degree 10, 14, or 18 terminal family is an automorphism.
  The proof exhausts all coefficient strata and boundary branches,
  including the exceptional degree-18 root-zero end, constant U, and
  reducible/nonreduced ambient response schemes.  It does not derive this
  chart for an arbitrary Keller pair, treat degree at least 26, or prove
  JC(2) or DC(2).
source: root/split-degrees-10-14-18-closure-2026-07-28
audit: >
  split-101418-hostile-audit-2026-07-28 independently reconstructed the
  finite boundary atlas and every Faber row through degree 18; audited all
  DVR slopes, zero-coordinate/no-active/cancellation branches, the unique
  degree-18 q^6 section pole, projective and vertical arguments, constant U,
  and component scope; and replayed 44 explicit gates under normal and
  optimized Python with byte-identical stored output and no assert nodes:
  PASS.
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2181-exact-square-prefix-compression-and-monic-depressed-quartic-closure
  - THM-2202-uniform-all-degree-quartic-pole-closure
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
related:
  - THM-2752-all-even-zero-first-flux-response-regularization-closure
  - THM-2759-split-degree-six-componentwise-monicization-closure
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2764-all-degree-chosen-sheet-even-zero-flux-componentwise-closure-away-from-six-mod-twelve
script: 04-computation/jc2_split_degrees10_14_18_componentwise_monicization_thm2767.py
output: 05-knowledge/results/jc2_split_degrees10_14_18_componentwise_monicization_thm2767.out
script_sha256: a44dbf13edc0c6f1b9e3015a6bb9ef970d174f4f81105055da5d24d16024339b
output_sha256: 46e991410602e1381dd01841f177170f8159bc5b8c18f84090934c83cd3eb089
hash_basis: LF-normalized bytes
---

# THM-2767 -- componentwise closure in split degrees 10, 14, and 18

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The argument below closes, rather than merely reduces, the complete
chosen-sheet split polynomial exact-square-prefix reduced Faber families in
degrees `10`, `14`, and `18`.  It is the degree-six mechanism of THM-2759
with two additions: finite exact augmented-ideal exclusion of every
nondegenerate top end, and a uniform sub-slope-four response-pole lemma for
the full odd/even bank.  The finite atlas needed to verify prefix exclusion,
response activity, and transversality is checked by the exact companion; the
DVR and projective steps are proved here.  The now independently audited
THM-2760 corroborates the prefix exclusion in all degrees, but remains
related rather than load-bearing because this proof reconstructs its finite
instances directly.

The former THM-2764 reservation advertised only the all-even, zero-first-
flux family away from `6 mod 12`; it is now superseded by the stronger
all-degree THM-2778.  The finite argument here remains an independently
audited proof in degrees `10,14,18`.

This does not derive the quartic terminal chart for an arbitrary Keller pair,
does not treat a nonsplit or nonpolynomial prefix, and does not yet claim any
degree at least `26`.

## 1. Exact family

Assume the chosen split sheet of the polynomial exact-square-prefix chart:

```text
P=(w^2+d)^2+(qw-s),
d=C-B^2/(4U^2),  q=A/U,  s=AB/(2U^2)-E.
```

In full Faber gauge a reduced mate of degree `M in {10,14,18}` is

```text
Q=E_M + sum_(1<=j<M, 4 does not divide j) a_j E_j,
```

with constant `a_j`.  The target shear in `C[P]` has zero response and is
irrelevant.  THM-2129 and THM-2723 give constants `lambda,W,kappa`, with
`kappa!=0`, such that

```text
Phi_Q=lambda,       Psi_Q=W,       U R_Q'=kappa.       (1)
```

If `U` is constant, `(x,z)->(x,w)` is a polynomial source automorphism and
THM-2181 makes `(P,Q)` an automorphism.  For contradiction, assume the only
other THM-2723 case:

```text
U=u_0(x-a)^m,       m>=2,
R_Q=r_0+r_1(x-a)^(1-m),       r_1!=0.                 (2)
```

Thus every response pole is the single finite point `x=a`, where all five
original prefix coefficients are regular.

Homogenize with `wt(h,d,q,s)=(1,2,3,4)`:

```text
F_(M+1)=Phi_M+sum a_j h^(M-j)Phi_j-lambda h^(M+1),
G_(M+2)=Psi_M+sum a_j h^(M-j)Psi_j-W h^(M+2),
R_(M+3)=R_M+sum a_j h^(M-j)R_j.                       (3)
```

The top forms `Phi_M,Psi_M` are coprime in each of the three degrees.  Hence
no member of `(3)` has a common surface component or a curve trapped in
`h=0`.  The nonconstant response in `(1)` makes the generic image a curve;
its reduced integral projective closure and normalization are the only
target objects used below.

## 2. Complete finite atlas of the top boundary

At `q=0`, `Phi_M=0` and `Psi_M` is a nonzero multiple of
`s^k`, where

```text
k=(M+2)/4 = 3,4,5.
```

Thus the only `q=0` top point is

```text
P_infty=[d:q:s]=[1:0:0].                              (4)
```

On the finite `q=1` index cover, exact lexicographic bases for
`(Phi_M,Psi_M)` are:

```text
M=10:
  d-3s^2,
  s^3-1/64.

M=14:
  d-160s^5+15s^2,
  s^6-s^3/10-1/6400.

M=18:
  d-(20480/7)s^8+896s^5-(65/7)s^2,
  s^10-(5/16)s^7+(5/1024)s^4-s/32768.                (5)
```

In all three quotient rings, both the top response and the Jacobian
`det d(Phi_M,Psi_M)/d(d,s)` are units.  Therefore every listed top point is
transverse, `h` is a local parameter on every response curve germ through
it, and `R_(M+3)/h^(M+3)` has a pole there.  The index-cover counts are
`3,6,10`; these are not asserted to be coarse weighted-orbit counts.

At a source point above one of these ends, the response pole and `(2)` place
the source point at `a`.  Put

```text
beta=hB/(2U).
```

The two polynomial exact-prefix identities are

```text
beta^2+d=h^2 C,             q beta-s=h^4 E.           (6)
```

On the `q=1` cover, `beta` is DVR-regular.  If `omega` is its residue, then

```text
d=-omega^2,                 s=omega.                  (7)
```

The load-bearing finite augmented-ideal calculation proves that no listed end
with `q omega!=0` can satisfy both top fluxes.  THM-2760 proves the same
exclusion uniformly and is a corroborating conceptual explanation.  Direct
substitution `d=-omega^2,s=q omega` gives the hard-checked common factors and
reduced pairs:

```text
M=10: gcd=q^3,
  (5 omega^2,  -5(32omega^3+q)/32).

M=14: gcd=q^4,
  (-7(80omega^3+q)/64, 35omega(16omega^3+q)/64).

M=18: gcd=omega q^5,
  (63(32omega^3+q)/128,
   -63omega(32omega^3+3q)/128).                       (8)
```

For `M=10,14`, the augmented ideals from `(5),(7)` are the unit ideal.  For
`M=18`, their common scheme is exactly

```text
d=s=omega=0,
```

i.e. the single point

```text
P_q=[d:q:s]=[0:1:0].                                  (9)
```

All other nonzero-`q` ends are therefore impossible.

## 3. The exceptional degree-18 point is killed by the original section

At a source point over `P_q`, `(6)` and regularity of `C` give
`v(beta)>0` on the homogeneous `q=1` cover.  If `alpha=v(h)>0`, then

```text
v(q_aff)=v(q/h^3)=-3alpha,
v(beta_aff)=v(beta/h)>-alpha.                          (10)
```

The exact polynomial value of the top Faber seed at the original section
`z=0` has the form

```text
E_18(0)=( ... -21 q_aff^6)/1024,                      (11)
```

where every omitted monomial is weighted-homogeneous of weight `18` in

```text
wt(beta_aff,C,E,q_aff)=(1,2,4,3)
```

and is not the pure `q_aff^6` monomial.  Hence `(11)` has the unique valuation
`-18alpha`.  Every monomial in every lower `E_j(0)` has valuation strictly
greater than `-j alpha`, or at worst `-j alpha`, and `j<18`.  Thus no lower
Faber row can cancel the pole.  This contradicts polynomiality of the
target-sheared `Q(x,0)`.  Point `(9)` is impossible.

The exact companion checks the coefficient `-21/1024`, weight `18` for every
top-section monomial, the unique pure-`q^6` support, and weighted degree `j`
for every monomial of every lower `E_j(0)` in the complete degree-18 bank.

This is exactly the degree-six `q_aff^2` mechanism of THM-2759, with
`q_aff^6` replacing `q_aff^2`.

## 4. Uniform full-bank slope-four lemma at `P_infty`

Work on the finite `d=1` index cover over `(4)`.  Let

```text
alpha=v(h)>0,       u=v(q),       v_0=v(s),
b=min(u,v_0),       k=(M+2)/4.                         (12)
```

The following statement includes every odd coefficient, every lower even
coefficient, `lambda`, `W`, unequal slopes, zero coordinates, and coefficient
specializations:

```text
b<4alpha  ==>  R_(M+3)/h^(M+3) has a pole.            (13)
```

The proof uses four exact row facts.

1. For every even row `j=4ell+2`, put `k_j=(j+2)/4`.  Both fluxes and the
   response have `(q,s)`-order `k_j` at `d=1`, while

   ```text
   ord_(q,s)(R_j+Phi_j/2)>=k_j+1
   ```

   (with the zero difference allowed), and

   ```text
   (M-j)+4k_j=M+2.
   ```

   Thus every lower even flux **and response** row is strictly later than the
   top order `kb` when `b<4alpha`.

2. The top leading faces are the odd and even parts of `(-s+iq)^k`; they are
   coprime.  `Psi_M` has a pure `s^k` term.  If `q` has smaller valuation,
   the pure `q^k` term is in `Psi_M` for even `k` and in `Phi_M` for odd `k`.

3. Every odd row has

   ```text
   Phi_j(1,0,0)=c_j!=0,       Psi_j=O(q,s),
   R_j(1,0,0)/c_j=(j+1)/(2(j+3)).                     (14)
   ```

   The gaps `M-j` are distinct.  Treat `-lambda h^(M+1)` as one additional
   Phi-only row of gap `M+1` and response ratio zero.

4. The all-degree even syzygy from THM-2752 gives

   ```text
   ord_(q,s)(R_M+Phi_M/2)>=k+1,                       (15)
   ```

   so the top response/first-flux ratio on the leading face is `-1/2`.

For these three degrees, every fact in this list is also reconstructed
directly from the finite Faber recurrence by the exact companion: the proof
does not need to inherit a finite identity from THM-2760 or from an unaudited
extrapolation of `(15)`.

Here is the exhaustive valuation argument.  If `s` has smaller valuation,
the pure `s^k` term is the unique top `Psi` term at order `kb`.  Any odd
`Psi` row early enough to cancel it has its constant `Phi` term strictly
before `kb`, where that term is unique; hence no branch exists.  The same
argument applies when `q` has smaller valuation and `k` is even.

Suppose `q` has smaller valuation and `k` is odd.  The pure `q^k` term lies
in `Phi`.  Let `g` be the least active odd/`lambda` gap.  If `g alpha` is
less or greater than `kb`, respectively the odd or top `Phi` term is unique.
At equality they can cancel, but `(14)--(15)` leave response multiplier

```text
-(1/2+(j+1)/(2(j+3)))=-(j+2)/(j+3) !=0               (16)
```

relative to the top first-flux coefficient; for the `lambda` row the
multiplier is `-1/2`.  Every lower even response row is strictly later by
item 1, every other odd response row has a larger gap, and the response
numerator therefore has exact order `kb<(M+3)alpha`.

Finally suppose `u=v_0=b`.  If the top `Phi` face vanishes, coprimality makes
the top `Psi` face nonzero.  An odd `Psi` row capable of reaching order `kb`
would put its constant `Phi` term strictly earlier, impossible.  If the top
`Phi` face is nonzero, a viable first equation forces `g alpha=kb`; then the
top `Psi` face must vanish because odd `Psi` starts one order later.  The
response again survives by `(16)`.  This proves `(13)`.

If **no odd coefficient and no `lambda` row is active**, there is no `g`:
for unequal slopes the appropriate pure top `Psi` term is unique (except the
odd-`k`, smaller-`q` case, where the pure top `Phi` term is unique), and for
equal slopes coprimality makes at least one top face nonzero with nothing at
order `kb` that can cancel it.  Thus the no-active case has no sub-slope
branch at all; it is not silently omitted from the response-ratio split.

By `(2)`, every pole in `(13)` is the finite point `a`, so the exact-prefix
identities `(6)` are available.  Their residues at `P_infty` give

```text
beta_0^2=-1.                                         (17)
```

If `b<4alpha`, the second identity forces

```text
u=v_0=b,                s_0=beta_0 q_0.              (18)
```

At `(18)`, both top faces are nonzero: when `beta_0=+i` one of
`-s+iq,-s-iq` vanishes and the other does not, and the odd/even parts of its
`k`-th power are both nonzero; the case `beta_0=-i` is symmetric.  The first
flux can cancel only against a least odd/`lambda` row with `g alpha=kb`.
Its `Psi` contribution has order `g alpha+b>kb`, all lower even rows are
later, and `W h^(M+2)` is later.  The nonzero top `Psi` term is unique, a
contradiction.  With no active odd/`lambda` row, the nonzero top `Phi` is
already unique, so this last step also covers that coefficient boundary.
Therefore every physical branch over `P_infty` satisfies

```text
min(v(q),v(s))>=4v(h),
v(q/h^3)>=v(h)>0.                                    (19)
```

Finite index covers multiply all valuations by one positive ramification
index, so `(13),(19)` descend to the coarse normalization.

## 4A. Normalization, boundary preimages, and ambient scheme structure

The coefficient functions define a nonconstant rational map from
`P1_x` to the projective response scheme `(3)`: if the generic coefficient
image were constant, then `R_Q` would be constant, contrary to `(1)`.  Let
`X` be the normalization of the reduced integral closure of that generic
image.  Properness extends the rational map uniquely to a morphism

```text
gamma:P1_x -> X.
```

It is a nonconstant morphism of complete integral curves, hence finite and
surjective.  Therefore every normalization point of `X` over `h=0`--that is,
every boundary branch used in Sections 2--4--has a nonempty source fibre.
The local response pole and exact-prefix arguments may consequently be
pulled back at every such branch; no abstract boundary point is being used
without a physical source preimage.

No integrality or reducedness of the whole ambient complete intersection is
assumed.  A map from the reduced source kills nilpotents and factors through
the target reduction, while its nonconstant generic image lies in one reduced
irreducible component.  Embedded components, nilpotents, and unrelated
irreducible components cannot change `(1)`, the exact identities `(6)`, or a
valuation on `X`.  Coprimality of the top forms excludes a generic-image
curve contained in `h=0` and excludes a common surface component.

## 5. Projective regularity forces the vertical branch

Sections 2 and 3 exclude every boundary point except `P_infty`; Section 4
makes `q_aff=q/h^3` regular and zero at every point above `P_infty`.  Every
projective image component meets `h=0`, since the top gcd excludes a boundary
component and a positive-dimensional projective curve cannot lie entirely in
the affine chart.  Thus `q_aff` is a global regular function on the complete
integral normalization and vanishes at a boundary point.  Hence

```text
q_aff=0.                                              (20)
```

## 6. The full vertical bank is triangular and impossible

Equation `(20)` gives `A=0` and `s=-E in C[x]`.  At `q=0`, every odd
`Psi_j` vanishes, while every even row `j=4ell+2` is a nonzero scalar times
`s^(ell+1)`.  The top coefficient is normalized to one, so the second flux
is a nonzero polynomial in `s` with constant coefficients.  Therefore
`s in C`.

For odd `j`, `Phi_j(d,0,s)` has degree exactly `(j+1)/2` in `d`, with
nonzero leading coefficient

```text
4 binom(j/2,(j+1)/2).                                 (21)
```

As `j=1,3,...,M-1`, these degrees are `1,2,...,M/2`.  Consequently the first
flux either gives a nonzero polynomial equation over `C` and forces `d in C`,
or is identically zero only when every odd `a_j` and `lambda` vanish.  In the
first case every response input is constant, contradicting `(1)`.  In the
second case only even response rows remain, and all are divisible by `q`, so
`R_Q=0`, again contradicting `(1)`.

The nonconstant-`U` branch is empty.  The constant-`U` branch is an
automorphism.  This closes the complete chosen-sheet split polynomial
exact-square-prefix reduced-degree `10`, `14`, and `18` families.

## 7. Hostile boundaries and cheapest rechecks

1. **THM-2760's scope is only the top exact-prefix exclusion.**  It is proved,
   cited, and verified-exact, with independent hostile audit pending.  It
   does not supply the top response-unit/transversality atlas, the root-zero
   section obstruction, or the lower-bank slope argument.  The finite
   augmented ideals here are therefore load-bearing; THM-2760 is
   corroborative rather than a dependency.
2. **Degree 18 really has an exception.**  The exact-prefix augmented ideal
   is not the unit ideal; it is exactly `(d,s,omega)`.  Omitting the
   `Q(x,0)` sidecar loses `P_q`.
3. **Pole capacity alone is insufficient.**  The proof needs exact-prefix
   regularity at the unique finite response pole and projective regularity of
   `q/h^3`.
4. **Odd rows are not discarded.**  They are essential in the response-pole
   lemma and the final vertical triangle.
5. **No total-integrality assumption.**  Work only on the reduced integral
   closure of the actual generic image; nilpotents die on pullback.
6. **Do not extrapolate to `M>=26` without auditing the top scheme.**  The
   all-degree prefix-resultant pattern looks favorable through `M=58`, but
   transversality and response nonvanishing at every top end have only been
   checked here for `10,14,18` (and are separately canonical at `22`).

Reproduction:

```bash
python3 04-computation/jc2_split_degrees10_14_18_componentwise_monicization_thm2767.py
python3 -O 04-computation/jc2_split_degrees10_14_18_componentwise_monicization_thm2767.py
```

The two outputs are byte-identical and match the transcript in
`05-knowledge/results`.  Every
finite identity used above is reconstructed directly from the recurrence,
including all lower odd/even rows in each complete bank, rather than merely
sampled from an all-degree claim.  The LF-normalized script and output hashes
are pinned in the front matter.
