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
  After noncanonical identifications G,B~=F_169, the joint Abel array
  partitions into 169 planar graphs beta=q^2/2+c. Conditional on retaining
  one graph, the restricted linear
  target--jet characters are quadratic chirps up to a scalar phase. At
  least one graph signal survives, and a graph with support at least two
  forces a nonzero refined joint coefficient A(q,z) with q nonzero. It
  does not force the coarser THM-2334 target coefficient
  C(q)=sum_z A(q,z): different graph/jet fibres at the same q can cancel.
  Graphwise chirp intensities and labelled singleton energies retain one
  independent global phase per graph. An exact F_169 hostile has a
  two-supported graph and nonzero refined q while C=delta_0; a same-data,
  same-total-current pair has zero versus positive coarse target
  Dirichlet energy. The exact average over first-jet characters is
  (1/169)sum_eta E_eta=sum_(q,z)W(q)|A(q,z)|^2, so
  any refined q!=0 survivor does force a
  nonzero edge defect for some exact factor-coloured jet-character probe.
  A nontrivial jet probe is not the original THM-2334 physical/coarse
  current and need not force the trivial jet polarizer. Cross-graph phase
  transport, jet depolarization, a target-preserving physical realization,
  a target row supported on one graph, or another coefficient-sensitive sidecar is
  still required. No scalar-row exclusion or LRC(14) closure is proved.
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
script_sha256: 56750547e1040d5898c72dd7a4f32f418bed779c224fe43a135f6ff8b3886238
output_sha256: b600f32c8caba8f3ea28058cc4bec274e679938572c2787f2d1c52f38ad42389
hash_basis: working-tree bytes (LF)
---

# THM-2356 -- finite-field chirp Gram tomography and planar graph slices

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2355 gives a local quadratic repair for the terminal-component phase
debt: singleton energies and cyclically twisted pair energies recover a
phase tree. There is a complementary fibrewise nonlinear repair. A planar
finite-group mask turns every translated pair inside one retained signal
into a distinct character. Double Fourier inversion then recovers all
off-diagonal Gram entries of that signal simultaneously.

The distinction between a linear mask and a planar mask is exact:

```text
linear mask:
  translates the Fourier label and retains the same autocorrelation;

planar mask:
  gives every nonzero displacement a bijective derivative and separates
  every ordered base point on that displacement.                     (1)
```

This is the first nonlinear response level that can destroy the
perfect-autocorrelation hostile inside one retained signal without
selecting individual component pairs. It does not compare the independent
global phases of different retained signals.

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

The leading factor has modulus one. Thus, **conditional on retaining the
graph-restricted current**, the linear target--jet characters give exactly
the full chirped intensity table of Section 2 for `Z_c`. This condition is
substantive: the unrestricted joint character amplitude is the sum of the
right side of (27) over every `c`, and does not isolate the graph label.
Squaring these amplitudes produces a derived quadratic functional of
THM-2337's coefficient array, not a pre-existing physical pair probe.

This corrects a tempting false dichotomy. Without graph restriction,
independent linear target and jet characters recover only the
autocorrelation of the joint array. On a retained graph the same linear
characters become nonlinear in `q`, because the graph itself is
quadratic.

Formula (15) has the following exact **refined** consequence. If some
surviving graph signal has two distinct supported targets `q_0,q_1`, then

```text
Z_c(q_0)conjugate(Z_c(q_1))!=0                     (28)
```

is one of the reconstructed off-diagonal Gram entries. Two distinct
targets cannot both be zero, so for some `q!=0` and some `z`,

```text
A(q,z)!=0.                                         (29)
```

This is a nonzero **joint target/jet fibre**. It is not yet a nonzero
coarse target fibre in the sense of THM-2334. Indeed, summing the joint
decomposition over the jet coordinate gives

```text
C(q):=sum_z A(q,z)=sum_c Z_c(q).                   (30)
```

The implication from (29) to `C(q)!=0` is false.

Choose `a,c in K` nonzero and define a joint array by its only nonzero
entries

```text
A(0,0)=1,
A(a,a^2/2)=1,
A(a,a^2/2+c)=-1.                                  (31)
```

Its graph signals are

```text
Z_0=delta_0+delta_a,       Z_c=-delta_a.            (32)
```

Thus `Z_0` has support two and off-diagonal Gram entry `1`, while

```text
C=delta_0.                                         (33)
```

The total current is nonzero, every coarse target twist is constant, and
no nonzero coarse target survives. This three-atom `F_169` hostile is
real-valued and uses the literal planar graph partition. Each active graph
is separately sign-definite, but the two graphs occupy opposite real
cones at the common target `a`; thus graphwise cone control does not
supply the familywise cone needed to prevent (30) from cancelling.
It is atom-minimal for this failure pattern: support two already spends
two nonzero joint entries on one graph, and cancelling its nonzero-target
entry in the coarse sum requires a third entry on a different graph.

The singleton boundary from Section 3 remains relevant to reconstructing
one graph. A graph signal supported at one point has constant chirp
intensities, independent of that point, by (18). Thus the chirp table
cannot distinguish a singleton at zero from a singleton at a nonzero
target. Labelled singleton graph energies

```text
|Z_c(q)|^2                                          (34)
```

locate it and, together with the chirp table, reconstruct each graph
signal up to a global phase. They do not glue those graph phases.

## 5. The exact cross-graph phase debt

For arbitrary phases `u_c in U(1)`, the transformation

```text
Z_c -> u_c Z_c                                     (35)
```

preserves every graphwise chirp intensity and every labelled singleton
energy. It changes the coarse sum (30). Hence complete tomography on each
graph has a `U(1)` gauge for every nonzero graph, not just one common
global gauge.

This loss can be written as an exact Dirichlet invoice. Identify
`K` additively with `F_13^2`, choose its two coordinate characters
`gamma_o,gamma_d`, and put

```text
W(q)=|gamma_o(q)-1|^2+|gamma_d(q)-1|^2.             (36)
```

The two characters separate points, so

```text
W(q)=0 iff q=0.
```

For the target transform

```text
T(ell)=sum_q chi_ell(q)C(q),
```

finite Parseval gives the two-generator target Dirichlet energy

```text
E_coarse
 =1/169 sum_ell (
    |T(ell+e_o)-T(ell)|^2
   +|T(ell+e_d)-T(ell)|^2)
 =sum_q W(q)|C(q)|^2.                              (37)
```

Consequently `E_coarse>0` is equivalent to a nonzero coarse target. The
sum of the separate graph energies is instead

```text
E_sep=sum_c sum_q W(q)|Z_c(q)|^2.
```

Expanding (30) gives

```text
E_coarse
 =E_sep
  +2 Re sum_(c<d) sum_q
      W(q)Z_c(q)conjugate(Z_d(q)).                 (38)
```

Planar chirp tomography supplies the within-graph Gram entries but none
of the cross-graph terms in (38). In (31), `E_sep>0` while
`E_coarse=0`.

The ambiguity remains even if the complete scalar current is retained.
Choose distinct nonzero graph labels `c,d` and a nonzero target `a`.
Compare the two graph families

```text
cancel:
  Z_0=delta_0-delta_a,
  Z_c=-delta_0+delta_a,
  Z_d=delta_0;

escape:
  Z_0=delta_0-delta_a,
  Z_c= delta_0-delta_a,
  Z_d=delta_0.                                     (39)
```

The two families have identical graphwise chirp intensities and labelled
singleton energies: only the global sign of `Z_c` changed. Their total
currents are both `1`. Their coarse target arrays are respectively

```text
C_cancel=delta_0,
C_escape=3delta_0-2delta_a,                         (40)
```

so (37) is zero for the first and positive for the second.

There is nevertheless an exact positive bridge at the **refined**
target/jet level. For a jet character `eta`, define

```text
C_eta(q)=sum_z eta(z)A(q,z),
H_eta(ell)=sum_q chi_ell(q)C_eta(q),                (41)
```

and let `E_eta` be the two-generator Dirichlet energy (37) of
`H_eta`. These are precisely the fixed-jet-character slices of
THM-2337's joint target--jet character bank. In graph coordinates,

```text
C_eta(q)
 =eta(q^2/2)sum_c eta(c)Z_c(q).                    (41a)
```

Thus the trivial character is exactly the coarse row sum, while the
nontrivial characters diagonalize its mean-zero kernel. Parseval first in
the target and then in the jet variable gives

```text
E_eta=sum_q W(q)|C_eta(q)|^2,

1/169 sum_eta E_eta
 =sum_(q,z) W(q)|A(q,z)|^2.                        (42)
```

Equivalently,

```text
E_1=E_coarse,             1/169 sum_eta E_eta=E_sep. (42a)
```

Because `W(q)>0` exactly for `q!=0`, the following are equivalent:

```text
A(q,z)!=0 for some q!=0;

E_eta>0 for some jet character eta;

some jet-polarized target current H_eta has a nonzero
  coordinate-edge defect.                          (43)
```

Thus the first Bockstein character is a lawful polarizer only in the
narrow, factor-coloured sense proved by THM-2337: it converts any nonzero
refined target coefficient into at least one exact jet-character
Dirichlet defect. For `eta!=1`, this probe reweights first-jet fibres and
is realized by independently shifted word factors. It is not the original
THM-2334 physical/coarse current and does not by itself certify the
canonical terminal-word current. In particular it does not force the
**trivial** polarizer `eta=1`. In the minimal hostile (31), the trivial
polarizer cancels at `a`, while any `eta` with `eta(c)!=1` leaves

```text
C_eta(a)=eta(a^2/2)(1-eta(c))!=0.                  (44)
```

The remaining passage is therefore jet depolarization `eta->1`,
graph-phase gluing, or a proof that some nontrivial probe has a
target-preserving physical realization—not construction of another
within-graph chirp.

The strongest valid loss ledger is therefore

```text
unrestricted linear target x jet twists
  -> joint autocorrelation only;

retain the planar graph label c=z-q^2/2
  -> a conditional quadratic-chirp table on every graph;

one graph has support at least two
  -> a nonzero refined joint coefficient A(q,z), q!=0;

jet-character averaging
  -> some jet-polarized magnetic edge defect;

graph chirps + labelled graph singleton energies
  -> each graph current up to its own global phase;

sum over graph/jet labels
  -> cross-graph phases can cancel every nonzero coarse target;

cross-graph phase transport, or a nonzero target occurring in exactly
  one active graph
  -> the missing gluing sidecar for eta=1.                         (45)
```

At a fixed target `q`, THM-2355 applies literally to the scalar components
`Z_c(q)`: labelled singleton magnitudes and cyclic pair-twist energies on
a spanning tree of the nonzero graph support recover their relative
phases and decide whether `C(q)` vanishes. Current canon supplies no
lawful graph-isolating pair probe. A common open half-plane for the active
`Z_c(q)` would also prevent cancellation, but physical factor positivity
does not by itself provide such a familywise cone.

Thus the remaining boundary is not merely one-sparse graph location. It
is the **refined-to-coarse graph-phase gluing problem**:

```text
prove C(q)!=0 for some q!=0,
or retain enough cross-graph phase to evaluate one such row sum.    (46)
```

THM-2337's termwise occupancy of every `(q,z)` does not settle (46),
because cancellation occurs both inside joint fibres and across the jet
fibres in (30). No scalar profile is excluded; the ledger remains `165`
and LRC(14) remains open.

## 6. The coarse row-sum kernel and the vertical refined-zero subspace

Write the graph coordinates rowwise as

```text
Z_q(c):=Z_c(q).
```

The coarse projection is the row sum

```text
C(q)=sum_c Z_q(c).                                  (47)
```

For `N=169`, its zero-only target space is therefore

```text
{Z:sum_c Z_q(c)=0 for every q!=0}.                 (48)
```

It decomposes as an arbitrary row at `q=0` plus one copy of
`1^perp subset C^N` for each of the other `N-1` rows. Its dimension is

```text
N+(N-1)^2=28,393.                                  (49)
```

The vertical tensors

```text
A(q,z)=delta_0(q)B(z)                              (50)
```

form only the first `N=169` dimensions. Even the one-sparse-per-graph
condition does not reduce the coarse kernel to (50): for distinct graph
labels `0,c,d` and `a!=0`,

```text
Z_0=delta_a,       Z_c=-delta_a,       Z_d=delta_0 (51)
```

has every active graph one-sparse and `C=delta_0`, but has refined
support at `q=a` and is not vertical.

What the vertical tensors classify exactly is the stronger
**refined-zero** condition

```text
A(q,z)=0 for every q!=0.                           (52)
```

Indeed (52) is equivalent to (50), and then every graph restriction is
`Z_c=B(c)delta_0`. The following hostile shows that this smaller boundary
is itself maximally robust.

### 6.1 A factorized, termwise-full refined-zero hostile

The vertical refined-zero boundary is maximally robust at the abstract
group-algebra level. Let `G` be any finite abelian group of order `N`,
fix `p in G`, and use THM-2333's rational endpoint weights

```text
U=delta_0+1_G,

V_0=delta_0-(1/(N+1))1_G,

V(x)=V_0(x-p).                                      (53)
```

Thus

```text
U(0)=2,             U(x)=1                  (x!=0),

V(p)=N/(N+1),       V(x)=-1/(N+1)           (x!=p). (54)
```

Both weights are nowhere zero. Put

```text
H(q)=sum_(u in G)U(u)V(u+p-q).                      (55)
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
H(q)=delta_0(q).                                    (56)
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
 =delta_0(q)B(z).                                  (57)
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
Z_c(q)=B(c)delta_0(q).                              (58)
```

A full-support `B` therefore makes all `169` graphs nonzero and
one-sparse at target zero.

This hostile also defeats two apparently stronger support arguments.
Take `B=1_K`. Its joint character transform is

```text
A_B^(xi,eta)
 =sum_(q,z)delta_0(q)xi(q)eta(z)
 =sum_z eta(z).                                    (59)
```

Its spatial and Fourier supports both have size `N`, so

```text
|support A_B| |support A_B^|
 =N^2
 =|K x K|.                                         (60)
```

It saturates the finite-group support uncertainty inequality. In target
coordinates `q=(x,y) in F_13^2`,

```text
delta_0(q)=(1-x^12)(1-y^12).                       (61)
```

As an interpolation polynomial on the two target and two jet coordinates,
(61) has separate degrees `(12,12,0,0)` and exactly

```text
1*1*13*13=169
```

nonzero grid values. It is the sharp THM-2285 interpolation-footprint
extremizer with surplus vector `(1,1,13,13)`. This is not a claim about
the number of monomials in its expanded polynomial.

Finally, the graph-chirp amplitude of (57) is

```text
M_c(a,b)
 =sum_q A_B(q,q^2/2+c)psi(bq+a(q^2/2+c))

 =B(c)psi(ac),                                     (62)
```

so

```text
|M_c(a,b)|^2=|B(c)|^2.                             (63)
```

For `B=1`, all `169^3` labelled graph-chirp intensities are one. Complete
termwise occupancy, full jet support, uncertainty equality, a sharp
polynomial footprint, and the complete intensity bank still do not
locate a singleton away from zero.

The scope is exact. Equations (53)--(63) form an abstract rational/complex
group-algebra hostile. They do not assert that THM-2337's canonical
interval coefficients realize `U,V,B`, or that a scalar LRC row realizes
this array. The signed endpoint weights, actual word factor, Abel-limit
compatibility, and terminal-component phase remain extra canonical
structure. The hostile rules out any proof using only the listed abstract
inputs, including any sidecar depending solely on `z`: arbitrary jet-only
data can be absorbed into `B`.

### 6.2 A graph-channel detector for refined support

There is a positive phase-sensitive coefficient functional which detects
nonzero refined target support inside one graph. For a graph signal put

```text
F_c(b)=sum_q Z_c(q)psi(bq),

S_c=F_c(0).
```

Normalized Parseval gives

```text
D_c
 :=1/N sum_b |F_c(b)-S_c|^2

 =sum_(q!=0)|Z_c(q)|^2
  +|sum_(q!=0)Z_c(q)|^2.                           (64)
```

Hence

```text
sum_c D_c>0
 iff A(q,z)!=0 for some q!=0.                      (65)
```

For a one-sparse graph `Z_c=a_c delta_(r_c)`,
`D_c=0` for `r_c=0` and

```text
D_c=2|a_c|^2
```

for `r_c!=0`. Two dual-basis phase ratios

```text
F_c(b)/F_c(0)=psi(b r_c)                            (66)
```

locate its refined target coordinate exactly. Intensities erase these ratios. The
minimal missing quadratic observable is

```text
F_c(b)conjugate(F_c(0)),                            (67)
```

which THM-2355 would recover from a lawful graph-channel pair twist.
Separate intensities `|F_c(b)|^2` do not determine (64), because they lose
the cross-frequency phase with `F_c(0)`. Current canon does not supply
that canonical pair probe. At coefficient level `D_c` is only a derived
quadratic functional of THM-2337. Thus (64)--(67) specify a refined-support
sidecar, not a coarse target-landing theorem. In the one-sparse hostile
(51), `sum_c D_c>0` while `E_coarse=0`.

## 7. Exact companion

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
  sharp finite-field footprint;
- verifies the three-atom refined/coarse hostile (31)--(33); and
- verifies that its trivial jet polarizer has zero target edge while an
  explicit separating jet character has a nonzero edge; and
- verifies the same-graph-data, same-total-current pair (39)--(40), with
  zero versus nonzero coarse target response.

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

## 8. Independent audits and correction lineage

The first independent audit rederived the two character-orthogonality steps,
the signs and normalization in (5) and (15), all `169` disjoint graph
slices after chosen identifications `G,B~=F_169`, and the scalar phase in
(27). It explicitly separated lawful coefficient grouping from a
pre-existing physical measurement and checked that the graph foliation is
basis-dependent. It also checked the vertical refined-zero hostile,
uncertainty equality, interpolation footprint, and `D_c` formula.

That audit nevertheless accepted the false implication from a nonzero
refined coefficient to the coarse row sum in commits `e0839fe1e` and
`59c933aae`. MISTAKE-261 retracts precisely that interpretation, not the
planar tomography theorem or those refined-support calculations.

A second hostile audit verified the three-atom counterexample, the
independent graph-phase gauge, both Dirichlet normalizations, the
`28,393`-dimensional coarse kernel, and the exact formula

```text
C_eta(q)=eta(q^2/2)sum_c eta(c)Z_q(c).
```

It confirmed that `eta=1` is the coarse row sum while the nontrivial
characters diagonalize its mean-zero kernel. Normal, optimized, and
stored transcripts and the refreshed LF hashes agree. The accepted scope
is: audited planar tomography; chosen noncanonical graph grouping;
derived refined-support functionals; exact refined-to-coarse hostiles;
and a precisely typed, still-open graph-phase/jet-depolarization debt.
No scalar profile is decremented. QED.
