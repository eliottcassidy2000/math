---
id: THM-2356
title: "Finite-field chirp Gram tomography and the Bockstein pairing debt"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. For equal
  finite abelian groups G,A and a planar map phi:G->A, the complete
  character-intensity table of the masks eta(phi(x)) reconstructs every
  off-diagonal entry z(x)conjugate(z(y)) by an explicit double Fourier
  inversion. Labelled singleton energies therefore reconstruct the full
  Gram matrix and z up to global phase. Over every odd finite field,
  phi(x)=x^2/2 is planar. Linear masks merely relabel the ordinary Fourier
  intensity and retain only autocorrelation. Chirp intensities alone
  cannot locate a singleton and have an exact unequal-magnitude
  two-support swap ambiguity; singleton energies are a sharp uniform
  sidecar. For the LRC target and first Bockstein spaces, both of order
  169, the theorem identifies a coefficient-respecting planar
  target--jet coupling as a sufficient global alternative to pair probes.
  Current canon supplies neither that coupling nor the aligned singleton
  ledger, so no target landing, scalar-row exclusion, or LRC(14) closure
  is proved.
source: codex-2026-07-25-finite-field-chirp-tomography
depends_on: []
related:
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
script: 04-computation/finite_field_chirp_gram_tomography_thm2356.py
output: 05-knowledge/results/finite_field_chirp_gram_tomography_thm2356.out
script_sha256: 9af9deda870754d520a9b040568bc4903c31562bbe351622087a150cc384f8b1
output_sha256: af261b877df50dae6f9fd22a2d16da0dffaed646eb37e4c6aa86b54329790145
hash_basis: working-tree bytes (LF)
---

# THM-2356 -- finite-field chirp Gram tomography

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2355 gives a local quadratic repair for the terminal-component phase
debt: singleton energies and cyclically twisted pair energies recover a
phase tree. There is a complementary global repair. A planar finite-group
mask turns every translated pair into a distinct character. Double Fourier
inversion then recovers all off-diagonal Gram entries simultaneously.

The distinction between a linear mask and a planar mask is exact:

```text
linear mask:
  translates the Fourier label and retains the same autocorrelation;

planar mask:
  gives every nonzero displacement a bijective derivative and separates
  every ordered base point on that displacement.                     (1)
```

This is the first nonlinear response level that can destroy the
perfect-autocorrelation hostile in THM-2344 without selecting individual
component pairs.

## 1. Planar masks reconstruct every off-diagonal Gram entry

Let `G,A` be finite abelian groups of the same order

```text
|G|=|A|=q.                                             (2)
```

Write their operations additively. A map

```text
phi:G->A
```

is **planar** if, for every nonzero `h in G`, the derivative

```text
D_h phi(y)=phi(y+h)-phi(y)                            (3)
```

is a bijection from `G` to `A`.

Let `z:G->C`. For characters

```text
eta in A^,                 chi in G^,
```

define the masked Fourier amplitude and its real intensity by

```text
F_z(eta,chi)
 =sum_(x in G) z(x) eta(phi(x)) chi(x),

E_z(eta,chi)=|F_z(eta,chi)|^2.                       (4)
```

Then for every `h!=0` and every `y in G`,

```text
z(y+h)conjugate(z(y))

 =1/q^2 sum_(eta in A^) sum_(chi in G^)
    E_z(eta,chi)
    conjugate(chi(h))
    conjugate(eta(D_h phi(y))).                      (5)
```

### Proof

Expand the intensity:

```text
E_z(eta,chi)
 =sum_(x,y)
   z(x)conjugate(z(y))
   eta(phi(x)-phi(y))
   chi(x-y).                                        (6)
```

Multiply by `conjugate(chi(h))` and average over `chi`. Character
orthogonality kills every term except `x-y=h`, leaving

```text
sum_y z(y+h)conjugate(z(y))
       eta(D_h phi(y)).                              (7)
```

Now multiply by `conjugate(eta(D_h phi(y_0)))` and average over `eta`.
Because `D_h phi` is bijective, character orthogonality kills every
summand in (7) except `y=y_0`. This is exactly (5). QED.

Thus the `q^2` real intensities recover the complete off-diagonal Hermitian
Gram data. If the labelled singleton energies

```text
s(y)=|z(y)|^2                                       (8)
```

are also retained, they supply the diagonal, so the full matrix

```text
Gamma_(x,y)=z(x)conjugate(z(y))                     (9)
```

is known. A rank-one positive semidefinite Gram matrix determines `z` up
to multiplication by one common element of `U(1)`. In particular it
determines every grouped current, every relative phase, every Fourier
coefficient, and every vanishing verdict.

The singleton ledger is redundant when `z` has support at least three.
Indeed, the nonzero off-diagonal entries first recover the support. For
three distinct supported points `x,y,w`,

```text
|z(x)|^2
 =Gamma_(x,y) Gamma_(w,x) / Gamma_(w,y).            (10)
```

The right side is defined and equals the positive real number on the left.

## 2. Odd finite fields supply the canonical planar chirp

Let `K=F_q` be a finite field of odd characteristic, let

```text
Tr:K->F_p
```

be its absolute trace, and fix the standard nontrivial additive character

```text
psi(t)=exp(2*pi*i Tr(t)/p).                          (11)
```

Identify `G=A` with the additive group of `K` and put

```text
phi(x)=x^2/2.                                       (12)
```

For every `h!=0`,

```text
D_h phi(y)=hy+h^2/2.                                (13)
```

Multiplication by `h` is a bijection of `K`, so (13) is planar. Writing
the two character parameters as `a,b in K`, the amplitudes are the
quadratic chirps

```text
F_z(a,b)
 =sum_(x in K)z(x)
   psi(bx+a x^2/2).                                 (14)
```

Formula (5) becomes the explicit field inversion

```text
z(y+h)conjugate(z(y))

 =1/q^2 sum_(a,b in K)|F_z(a,b)|^2
   psi(-b h-a(hy+h^2/2)),             h!=0.         (15)
```

No genericity, full-support, or nonvanishing assumption enters (15).
The construction is a finite-field coded-diffraction identity, but the
proof is only the two orthogonality steps above.

## 3. Linear masks and missing singleton magnitudes are sharp boundaries

If `phi` is a homomorphism, then `eta composed with phi` is a character
of `G`. Hence

```text
F_z(eta,chi)=z_hat(chi (eta composed with phi)),    (16)
```

up to the chosen Fourier sign. Varying `eta` only relabels the ordinary
Fourier intensities. Their inverse transform is the periodic
autorrelation

```text
C(h)=sum_y z(y+h)conjugate(z(y)),                   (17)
```

not the labelled pair products. Thus arbitrary banks of independent
linear target and jet characters do not implement (5). The load-bearing
coordinate is the bijective nonlinear derivative (3).

The diagonal sidecar in (8) is also sharp for a uniform reconstruction
theorem. For every `x in G`,

```text
z=delta_x
```

has

```text
E_z(eta,chi)=1                                     (18)
```

for every mask and character. Chirp intensities alone cannot locate a
singleton.

There is a second exact ambiguity on a fixed two-point support. At distinct
points `x_0,x_1`, compare the real signals

```text
z(x_0)=2, z(x_1)=3,

z'(x_0)=3, z'(x_1)=2.                               (19)
```

Both have total diagonal mass `13`, and every ordered off-diagonal product
equals `6`. Equation (6) therefore gives

```text
E_z(eta,chi)=E_(z')(eta,chi)                        (20)
```

for every planar mask, although `z,z'` are not common-phase multiples.
The labelled singleton energies distinguish them. These hostiles explain
exactly why (5) recovers all off-diagonal data but not the separate
diagonal entries.

For support at least three, (10) removes this ambiguity. For support one
or two, the intensity table still decides many phase-insensitive verdicts,
but it does not reconstruct the labelled current without (8).

## 4. The exact LRC target--jet sidecar

THM-2334's target quotient is

```text
G isomorphic to F_13^2,
```

and THM-2337's ordered first Bockstein target-jet space is

```text
B isomorphic to F_13^2.                             (21)
```

Each has `169` elements. After choosing compatible bases, both additive
groups may be identified with

```text
K=F_169.                                            (22)
```

Let `Z(q)` denote the current grouped by target residue. Existing target
characters provide the ordinary Fourier amplitudes

```text
sum_q Z(q)psi(bq).                                  (23)
```

Their energies recover only the autocorrelation in THM-2355. Independent
linear characters of the Bockstein jet do not change this conclusion:
Section 3 shows that a linear response bank can only relabel (23).

The sufficient nonlinear service is now completely explicit:

```text
for every a,b in F_169, retain the lawful semantic amplitude

sum_q Z(q)psi(bq+a q^2/2),                          (24)

together with the labelled singleton target energies |Z(q)|^2.         (25)
```

Equations (15) and (25) reconstruct the whole target Gram matrix and hence
the target current up to common phase. In particular they eliminate the
real full-support perfect-autocorrelation hostile of THM-2355.

There is an equivalent target--jet formulation. A coefficient-respecting
planar section

```text
q -> beta(q)=q^2/2 in B                             (26)
```

would let the linear jet character `a` produce precisely the chirp in
(24). More generally any weight-preserving coupling whose derivative in
`q` is a permutation supplies the planar hypothesis of Section 1.

THM-2337 proves only that every target and every first-jet value is
termwise occupied. Its exact address gauge changes the jet while fixing
the address and is not weight preserving. Surjective occupancy therefore
does not supply the coefficient identity (26), a semantic restriction to
its graph, or the amplitudes (24). Nor does current canon align the
component singleton energies of THM-2303 with the labelled target energies
in (25).

The new loss ledger is

```text
ordinary target twists
  -> target autocorrelation only;

independent linear target x jet twists
  -> relabelled joint autocorrelations only;

coefficient-respecting planar target--jet coupling
  -> every off-diagonal target Gram entry;

planar coupling + labelled target singleton energies
  -> complete target current up to common phase.                    (27)
```

Thus THM-2356 gives a precise global alternative to THM-2355's lawful
pair-twist tree. The next canonical LRC task is no longer “find more
characters.” It is to construct the nonlinear planar coupling (26), or to
prove that the actual interval weights confine the target current to a
phase cone. No such construction is claimed here; the ledger remains
`165` and LRC(14) remains open.

## 5. Exact companion

The dependency-free companion works in

```text
F_169=F_13[theta]/(theta^2-2)
```

and represents all cyclotomic values exactly in `Q(i,zeta_13)`. It:

- exhausts all `168` nonzero derivatives of `x^2/2` on `F_169`;
- constructs the complete `28,561`-entry chirped intensity table for a
  five-site Gaussian-integer signal;
- performs the literal double inversion on every ordered supported pair
  and five zero controls;
- checks all `28,561` linear-mask relabellings; and
- verifies the `169` invisible singleton locations and the two-site
  unequal-magnitude swap.

Reproduce with

```bash
python3 04-computation/finite_field_chirp_gram_tomography_thm2356.py
python3 -O 04-computation/finite_field_chirp_gram_tomography_thm2356.py
```

Both transcripts must match

```text
05-knowledge/results/finite_field_chirp_gram_tomography_thm2356.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent audit is pending. QED.
