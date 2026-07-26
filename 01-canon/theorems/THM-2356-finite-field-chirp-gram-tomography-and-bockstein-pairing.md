---
id: THM-2356
title: "Finite-field chirp Gram tomography and Bockstein planar graph slices"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For equal
  finite abelian groups G,A and a planar map phi:G->A, the complete
  character-intensity table of the masks eta(phi(x)) reconstructs every
  off-diagonal entry z(x)conjugate(z(y)) by an explicit double Fourier
  inversion. Labelled singleton energies therefore reconstruct the full
  Gram matrix and z up to global phase. Over every odd finite field,
  phi(x)=x^2/2 is planar. Linear masks merely relabel the ordinary Fourier
  intensity and retain only autocorrelation. Chirp intensities alone
  cannot locate a singleton and have an exact unequal-magnitude
  two-support swap ambiguity; singleton energies are a sharp uniform
  sidecar. The LRC target and first Bockstein spaces both have order 169.
  After choosing identifications G,B~=F_169, their joint Abel array
  partitions lawfully but noncanonically into the 169 planar graphs
  beta=q^2/2+c; on each graph the linear target--jet characters become
  quadratic chirps up to a scalar phase. Their intensities are derived
  quadratic functionals, not pre-existing physical probes. At least one graph signal
  survives. Any graph with target support at least two therefore forces a
  nonzero target through its recovered off-diagonal Gram entry. The exact
  no-landing locus inside the one-sparse boundary is the vertical tensor
  A(q,z)=delta_0(q)B(z). A factorized full-support hostile can occupy all
  28,561 joint fibres termwise and make all 169 planar graphs survive
  while every graph is a singleton at q=0; it simultaneously saturates
  finite-group uncertainty and the separate-degree footprint bound. The
  missing observable is a graph-channel phase ratio or equivalent
  target--jet coefficient coupling, not construction of the chirp service
  itself. No scalar-row exclusion or LRC(14) closure is proved.
source: codex-2026-07-25-finite-field-chirp-tomography
depends_on:
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
related:
  - THM-2285-centered-grid-footprint-and-generic-keller-lines
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2333-abel-target-fibre-sum-landing-and-zero-fibre-boundary
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
script: 04-computation/finite_field_chirp_gram_tomography_thm2356.py
output: 05-knowledge/results/finite_field_chirp_gram_tomography_thm2356.out
script_sha256: f83cfbbafd56ecd8ae82554a80b30e588c3fdf005f3d5af94ff47e1ed0770949
output_sha256: 22f6cc95cde44f7dc7d65b31dfe9eb7979cefc6467cece28eb2af16ae39f146e
hash_basis: working-tree bytes (LF)
---

# THM-2356 -- finite-field chirp Gram tomography and planar graph slices

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

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

Let

```text
A(q,z),                    (q,z) in G x B,           (23)
```

be THM-2337's full-semantic joint target/first-jet Abel array. Every joint
limit exists and the array is not identically zero.

For each `c in K`, define its planar graph signal

```text
Z_c(q)=A(q,q^2/2+c).                                (24)
```

The `169` graphs in (24) are pairwise disjoint and partition `K x K`:
every `(q,z)` lies on the unique graph

```text
c=z-q^2/2.                                         (25)
```

Consequently

```text
A is not zero
 => Z_c is not zero for at least one c.             (26)
```

This is a coefficient-respecting partition of the already existing joint
array, not a gauge choice or a termwise occupancy assertion. It is lawful
after the two identifications in (22), but those identifications and hence
the planar foliation are not canonical LRC data.

Restrict a joint linear target--jet character to one graph. Its amplitude
is

```text
sum_q A(q,q^2/2+c) psi(bq+a(q^2/2+c))

 =psi(ac) sum_q Z_c(q)psi(bq+a q^2/2).              (27)
```

The leading factor has modulus one. Thus finite coefficient grouping by
the graph label followed by the existing linear characters gives the full
chirp amplitudes of Section 2 for `Z_c`. Squaring their moduli produces a
derived quadratic functional of THM-2337's coefficient array. No theorem
here asserts that this intensity table is already a physical pair-probe
service in the canonical LRC construction.

This corrects a tempting false dichotomy. Without graph restriction,
independent linear target and jet characters recover only the
autocorrelation of the joint array. On a retained graph the same linear
characters become nonlinear in `q`, because the graph itself is
quadratic.

Formula (15) now has a direct LRC consequence. If some surviving graph
signal has two distinct supported targets `q_0,q_1`, then

```text
Z_c(q_0)conjugate(Z_c(q_1))!=0                     (28)
```

is one of the reconstructed off-diagonal Gram entries. Two distinct
targets cannot both be zero, so some nonzero full target survives.
Equivalently,

```text
no nonzero target survives
 => every Z_c is supported in {0}.                  (29)
```

More generally, if the chirp table has any nonzero off-diagonal Gram
entry, target landing follows without a singleton-energy sidecar.

The residual case is sharp. A graph signal supported at one point has
constant chirp intensities, independent of that point, by (18). Thus the
chirp table cannot distinguish

```text
support(Z_c)={0}
```

from a singleton at a nonzero target. Labelled singleton graph energies

```text
|Z_c(q)|^2                                          (30)
```

would locate it and complete the reconstruction, but current canon does
not align THM-2303's component magnitudes with (30).

The new loss ledger is

```text
unrestricted linear target x jet twists
  -> joint autocorrelation only;

retain the planar graph label c=z-q^2/2
  -> lawful quadratic chirps on every graph;

some graph has support at least two
  -> reconstructed off-diagonal Gram entry
  -> nonzero target;

every graph is one-sparse
  -> chirp intensities lose its singleton location;

planar graph chirps + labelled graph singleton energies
  -> complete graph currents up to phase.                          (31)
```

Thus THM-2356 supplies a finite-field algebraic alternative to THM-2355's
lawful pair-twist tree, not merely a proposed identity. It reduces the LRC
phase/target problem to the **planar graph singleton boundary**:

```text
force one graph to contain two surviving target coefficients,
or locate a surviving singleton away from q=0.                      (32)
```

THM-2337's termwise occupancy of every `(q,z)` does not settle (32),
because the Abel sums may cancel separately in every joint fibre. No
scalar profile is excluded; the ledger remains `165` and LRC(14) remains
open.

## 5. The singleton boundary is exactly a vertical tensor

The one-sparse boundary has a complete classification. If every graph
signal is at most one-sparse, write uniquely

```text
Z_c=0,

or

Z_c=a_c delta_(r_c),       a_c!=0, r_c in K.        (33)
```

If some nonzero graph has `r_c!=0`, then

```text
A(r_c,r_c^2/2+c)=a_c!=0                            (34)
```

is already a surviving nonzero target. If no such graph exists, every
nonzero graph has `r_c=0`, and necessarily

```text
A(q,z)=delta_0(q)B(z),       B(c):=a_c,             (35)
```

where `a_c=0` on a zero graph. Conversely, every array of the form (35)
has

```text
Z_c(q)=A(q,q^2/2+c)=B(c)delta_0(q).                 (36)
```

Thus the exact no-landing locus inside the one-sparse boundary is the
vertical tensor ideal

```text
delta_0 tensor B.                                  (37)
```

It is not an unclassified family of missing singleton positions.

### 5.1 A factorized, termwise-full hostile

The vertical tensor boundary is maximally robust at the abstract
group-algebra level. Let `G` be any finite abelian group of order `N`,
fix `p in G`, and use THM-2333's rational endpoint weights

```text
U=delta_0+1_G,

V_0=delta_0-(1/(N+1))1_G,

V(x)=V_0(x-p).                                      (38)
```

Thus

```text
U(0)=2,             U(x)=1                  (x!=0),

V(p)=N/(N+1),       V(x)=-1/(N+1)           (x!=p). (39)
```

Both weights are nowhere zero. Put

```text
H(q)=sum_(u in G)U(u)V(u+p-q).                      (40)
```

Every sum contains exactly `N` nonzero products. At `q=0`, its numerator
over `N+1` is

```text
2N-(N-1)=N+1,
```

whereas at `q!=0` it is

```text
-2+N-(N-2)=0.
```

Therefore

```text
H(q)=delta_0(q).                                    (41)
```

Equivalently, `V_0` is the correlation inverse of `U`.

Now let `B:K->C` be arbitrary and define the atomic joint terms

```text
T_B(u;q,z)=U(u)V(u+p-q)B(z)
```

and their aggregate

```text
A_B(q,z)
 =sum_u T_B(u;q,z)
 =H(q)B(z)
 =delta_0(q)B(z).                                  (42)
```

If `B(z)!=0`, the coefficient at `(q,z)` has exactly `N` nonzero atomic
endpoint pairs. If `B` has full support and `N=169`, all

```text
N^2=28,561
```

joint fibres are termwise occupied, for

```text
N^3=4,826,809
```

nonzero atomic incidences. Nevertheless only the `N` coefficients
`(0,z)` survive aggregation. Every planar restriction is

```text
Z_c(q)=B(c)delta_0(q).                              (43)
```

A full-support `B` therefore makes all `169` graphs nonzero and
one-sparse at target zero.

This hostile also defeats two apparently stronger support arguments.
Take `B=1_K`. Its joint character transform is

```text
A_B^(xi,eta)
 =sum_(q,z)delta_0(q)xi(q)eta(z)
 =sum_z eta(z).                                    (44)
```

Its spatial and Fourier supports both have size `N`, so

```text
|support A_B| |support A_B^|
 =N^2
 =|K x K|.                                         (45)
```

It saturates the finite-group support uncertainty inequality. In target
coordinates `q=(x,y) in F_13^2`,

```text
delta_0(q)=(1-x^12)(1-y^12).                       (46)
```

As an interpolation polynomial on the two target and two jet coordinates,
(46) has separate degrees `(12,12,0,0)` and exactly

```text
1*1*13*13=169
```

nonzero grid values. It is the sharp THM-2285 interpolation-footprint
extremizer with surplus vector `(1,1,13,13)`. This is not a claim about
the number of monomials in its expanded polynomial.

Finally, the graph-chirp amplitude of (42) is

```text
M_c(a,b)
 =sum_q A_B(q,q^2/2+c)psi(bq+a(q^2/2+c))

 =B(c)psi(ac),                                     (47)
```

so

```text
|M_c(a,b)|^2=|B(c)|^2.                             (48)
```

For `B=1`, all `169^3` labelled graph-chirp intensities are one. Complete
termwise occupancy, full jet support, uncertainty equality, a sharp
polynomial footprint, and the complete intensity bank still do not
locate a singleton away from zero.

The scope is exact. Equations (38)--(48) form an abstract rational/complex
group-algebra hostile. They do not assert that THM-2337's canonical
interval coefficients realize `U,V,B`, or that a scalar LRC row realizes
this array. The signed endpoint weights, actual word factor, Abel-limit
compatibility, and terminal-component phase remain extra canonical
structure. The hostile rules out any proof using only the listed abstract
inputs, including any sidecar depending solely on `z`: arbitrary jet-only
data can be absorbed into `B`.

### 5.2 The exact missing graph-channel observable

There is a positive phase-sensitive coefficient functional which detects
target landing without a sparsity assumption. For a graph signal put

```text
F_c(b)=sum_q Z_c(q)psi(bq),

S_c=F_c(0).
```

Normalized Parseval gives

```text
D_c
 :=1/N sum_b |F_c(b)-S_c|^2

 =sum_(q!=0)|Z_c(q)|^2
  +|sum_(q!=0)Z_c(q)|^2.                           (49)
```

Hence

```text
sum_c D_c>0
 iff some nonzero target survives.                 (50)
```

Under (33), `D_c=0` for `r_c=0` and

```text
D_c=2|a_c|^2
```

for `r_c!=0`. Two dual-basis phase ratios

```text
F_c(b)/F_c(0)=psi(b r_c)                            (51)
```

locate a one-sparse target exactly. Intensities erase these ratios. The
minimal missing quadratic observable is

```text
F_c(b)conjugate(F_c(0)),                            (52)
```

which THM-2355 would recover from a lawful graph-channel pair twist.
Separate intensities `|F_c(b)|^2` do not determine (49), because they lose
the cross-frequency phase with `F_c(0)`. Current canon does not supply
that canonical pair probe. At coefficient level `D_c` is only a derived
quadratic functional of THM-2337. Thus (49)--(52) are a precise sidecar
specification, not a physical target-landing theorem.

## 6. Exact companion

The dependency-free companion works in

```text
F_169=F_13[theta]/(theta^2-2)
```

and represents all cyclotomic values exactly in `Q(i,zeta_13)`. It:

- exhausts all `168` nonzero derivatives of `x^2/2` on `F_169`;
- exhausts the `169` disjoint planar graph translates partitioning all
  `28,561` target--jet pairs;
- constructs the complete `28,561`-entry chirped intensity table for a
  five-site Gaussian-integer signal;
- performs the literal double inversion on every ordered supported pair
  and five zero controls;
- checks all `28,561` linear-mask relabellings; and
- verifies the `169` invisible singleton locations, the two-site
  unequal-magnitude swap, the factorized vertical tensor with all
  `28,561` joint fibres termwise occupied, uncertainty equality, and the
  sharp finite-field footprint.

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

## 7. Independent audit

The independent audit rederived the two character-orthogonality steps,
the signs and normalization in (5) and (15), all `169` disjoint graph
slices after chosen identifications `G,B~=F_169`, and the scalar phase in
(27). It explicitly separated lawful coefficient grouping from a
pre-existing physical measurement and checked that the graph foliation is
basis-dependent.

The same audit proved the at-most-one-sparse dichotomy (33)--(37),
recomputed the factorized vertical tensor and every term count in
(38)--(43), and checked the uncertainty and interpolation-footprint
equalities without treating (46) as a monomial-support assertion. It
independently verified `35` Gaussian-integer instances of (49) and `373`
vertical cells, confirming that `D_c>0` is equivalent to nonzero-target
support while separate chirp intensities do not determine `D_c`.

Normal, optimized, and stored transcripts, LF hashes, and the exact
companion all agree. The accepted scope is: chosen finite-field
identifications; lawful but noncanonical graph grouping; derived
quadratic chirp and `D_c` functionals; and a formal group-algebra hostile.
No scalar profile is decremented until canonical structure excludes the
vertical tensor or realizes the graph-channel phase ratio. QED.
