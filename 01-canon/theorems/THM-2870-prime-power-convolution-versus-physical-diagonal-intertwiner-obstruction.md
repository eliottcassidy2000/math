---
id: THM-2870
title: "Prime-power convolution versus physical-diagonal intertwiner obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  If an
  invertible operator C and a projection D act on the same
  finite-dimensional space, every intertwiner CT=TD factors from im(D)
  through ker(C-1), while every reverse intertwiner DT=TC factors from
  coker(C-1) into im(D).  For a nontrivial finite p-group, a Boolean mask
  of p-unit mass is therefore an invertible convolution kernel but a
  proper singular physical diagonal, so no invertible same-mask
  intertwiner exists.  The result explains why the full Cayley spectra
  in THM-2852 do not supply a physical LRC action.
source: root/cayley-physical-intertwiner-2026-07-28
audit: >
  thm2870-two-orientation-intertwiner-audit-2026-07-28 (independent
  factorization and rank-nullity derivation, modular/nonsemisimple scope,
  prime-power Boolean and Fourier typing, Cayley C5/C9 replay, sharp C3
  maps, and mixed-prime C6 hostile: ACCEPT)
depends_on:
  - THM-2852-prime-power-orbit-spectrum-harvest-and-cayley-tournament-nonsingularity
related:
  - THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2602-commutative-vertex-insertion-and-ordered-transition-curvature-no-go
script: 04-computation/prime_power_convolution_physical_intertwiner_thm2870.py
output: 05-knowledge/results/prime_power_convolution_physical_intertwiner_thm2870.out
script_sha256: 8db743b70bb3b7942abcdc77b4561741937435dd65792d0332855cd203fed823
output_sha256: a946a81bda5dde97869b40d3fed9fcca760e0f2f018d68eecd4357eab6014e44
hash_basis: LF-normalized bytes
---

# THM-2870 -- convolutional fullness is not a physical diagonal action

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2839 and THM-2852 prove that a `p`-integral mask of `p`-unit mass on a
finite `p`-group is a unit convolution kernel.  In the LRC application,
however, the same mask is physically a diagonal support predicate.
Convolutional fullness has repeatedly looked close to the missing
translation action.  The following elementary classification shows that
the two uses of the mask are maximally different.

## 1. Projection--unit intertwiner classification

Let `K` be a field, let `V` be a finite-dimensional `K`-vector space, let

```text
D^2=D in End_K(V),        C in GL_K(V),                 (1)
```

and put

```text
V_1=im(D),        V_0=ker(D),
E_1(C)=ker(C-1),
Q_1(C)=V/(C-1)V,
r=dim(V_1),       e=dim(E_1(C)).                        (2)
```

Rank--nullity gives

```text
dim Q_1(C)=e                                             (3)
```

without any semisimplicity hypothesis.

There are two possible intertwining orientations.  They admit exact
factorizations:

```text
{T: C T=T D}  = { i_E A D:
                   A in Hom_K(V_1,E_1(C)) },             (4)

{T: D T=T C}  = { i_1 B pi:
                   B in Hom_K(Q_1(C),V_1) }.             (5)
```

Here `i_E:E_1(C)->V` and `i_1:V_1->V` are the inclusions and
`pi:V->Q_1(C)` is the quotient.

Indeed, if `CT=TD` and `v in V_0`, then `CTv=0`; invertibility of `C`
forces `Tv=0`.  For `v in V_1`, one has `CTv=Tv`, so `Tv in E_1(C)`.
This proves `(4)`, including its converse.  If `DT=TC`, then

```text
(1-D)TC=0.
```

Right composition by `C^(-1)` gives `(1-D)T=0`, so the image lies in
`V_1`.  The original equation then becomes `T(C-1)=0`, which is exactly
the quotient factorization `(5)`.

Consequently

```text
dim {T:CT=TD}=dim {T:DT=TC}=r e,                         (6)
rank(T)<=min(r,e)                                        (7)
```

in either orientation.  In particular, if `D` is a proper projection,
there is no invertible intertwiner: in `(4)` it kills `V_0`, and in `(5)`
its image lies in `V_1`.

## 2. One Boolean mask, two inequivalent actions

Let `G` be a finite group and let `S subset G`.  On `V=K^G`, the Boolean
mask `h=1_S` has two natural actions:

```text
(D_h f)(x)=h(x)f(x),                 pointwise/physical;
R_h f=f*h,                           regular convolution. (8)
```

Thus

```text
rank(D_h)=|S|,       nullity(D_h)=|G|-|S|.               (9)
```

Assume now that `G` is a nontrivial finite `p`-group, that `h` is viewed
in `Z_(p)[G]`, and that

```text
p does not divide |S|.                                  (10)
```

THM-2839, in the exact form used by THM-2852, says that `R_h` is a unit
over `Z_(p)`, hence over `Q` and `C`.  Since `p` divides `|G|`, condition
`(10)` also forces

```text
empty != S != G.                                        (11)
```

Apply `(4)--(7)` with `C=R_h` and `D=D_h`.  There is no invertible
same-mask intertwiner in either direction:

```text
R_h T=T D_h       or       D_h T=T R_h.                 (12)
```

More precisely, every forward intertwiner kills all physical coordinates
outside `S` and lands in the convolutional `1`-eigenspace; every reverse
intertwiner is supported physically inside `S` and kills `(R_h-1)V`.

For abelian `G` over a splitting field of characteristic zero, using the
unnormalized Fourier convention

```text
hhat(chi)=sum_(g in G) h(g)chi(g),                       (13)
```

one has

```text
e=#{chi in G^: hhat(chi)=1}.                             (14)
```

Equation `(14)` is only a convenient abelian description.  The
kernel/cokernel classification `(4)--(7)` remains valid in modular
characteristic and for nonabelian groups.

## 3. Why Fourier inversion is not the missing physical sidecar

For abelian `G`, under the chosen Fourier/inversion convention a Fourier
transform gives

```text
F R_h F^(-1)=D_(hhat),                                  (15)
```

not `D_h`.  The opposite convention composes `hhat` with character
inversion, which changes none of the conclusions.  Under `(10)`, every
value of `hhat` is nonzero because `R_h`
is a unit, whereas a proper Boolean `h` has zeros by `(11)`.  Hence no
permutation, diagonal gauge, or invertible linear change of basis can
turn `(15)` into a same-mask identification.

The information destroyed in passing from the physical diagonal to
convolution is therefore not rank.  It is the coefficient/character
reference:

```text
coefficient label g  <---- missing reference ---->  character label chi.
                                                               (16)
```

A lawful sidecar would have to construct that reference from the physical
source/target action and prove that the resulting diagonal is
`D_(hhat)` with the required phase, not merely observe that `R_h` is
nonsingular.

## 4. Cayley tournaments

Let `p` be odd, `|G|=p^d`, and let `S` choose exactly one element from
each inverse pair in `G\{1}`.  Then

```text
|S|=(p^d-1)/2=-1/2 mod p.                               (17)
```

THM-2852 proves that the Cayley-tournament adjacency matrix is, up to
transpose convention, the regular matrix of `R_h` and is nonsingular.
The physical diagonal has

```text
rank(D_h)=(p^d-1)/2,
nullity(D_h)=(p^d+1)/2.                                 (18)
```

Thus tournament nonsingularity supplies no same-mask physical
orientation.  Every linear bridge satisfying either equation in `(12)`
has rank at most

```text
min((p^d-1)/2, dim ker(R_h-1)).                          (19)
```

This is a precise stopping rule for importing tournament adjacency into
the LRC diagonal packet.

## 5. Sharp boundaries and controls

1. The bound can be attained.  On `G=C_p`, let `S={g}` for a generator
   `g`.  Then `R_h` is the cyclic shift and `e=1`.  The maps

   ```text
   T_+(f)=f(g) * 1_G,
   T_-(f)=(sum_x f(x)) delta_g                              (20)
   ```

   have rank one and realize the two orientations in `(12)`.

2. Invertibility of `C` is load-bearing.  For an abstract proper
   projection, taking `C=D` and `T=1` gives an invertible intertwiner;
   the cancellation by `C^(-1)` in the proof of `(5)` is unavailable.

3. The `p`-group and `p`-integrality hypotheses are load-bearing for
   obtaining convolutional invertibility.  THM-2852's `C_6` mask
   `1+Z^2+Z^4` and its rational normalized subgroup mask are the sharp
   counterexamples.

4. A nonnegative inverse would be stronger than an intertwiner.
   THM-2852 proves that a nonnegative mask has a nonnegative convolution
   inverse only in the monomial case.

5. Nothing here constructs a physical translation, positive endpoint
   current, target/carry intertwiner, row exclusion, or LRC decrement.

## 6. Exact evidence and status

The exact companion enumerates every `5`-unit Boolean mask on `C_5`, all
Cayley tournaments on `C_5` and `C_9`, and the two sharp rank-one maps on
`C_3`.  Over the relevant prime field it verifies convolutional
invertibility, both Sylvester-intertwiner nullities in `(6)`, the
rank/nullity split `(9)`, and the `C_6` mixed-prime hostile.  It must use
explicit exceptions rather than truth-bearing Python `assert` statements.
It also exhausts all `144` pairs consisting of an invertible operator on
`F_3^2` and a standard coordinate projection.  Run

```text
python 04-computation/prime_power_convolution_physical_intertwiner_thm2870.py
python -O 04-computation/prime_power_convolution_physical_intertwiner_thm2870.py
```

Both modes LF-normalized byte-match

```text
05-knowledge/results/prime_power_convolution_physical_intertwiner_thm2870.out.
```

```text
PROVED IN THE CANDIDATE: exact two-orientation factorization;
                         Hom-space dimensions and rank bound;
                         prime-power same-mask no-go;
                         abelian Fourier-label boundary;
                         Cayley-tournament specialization;
                         sharp rank-one intertwiners.

NOT PROVED:              a coefficient/character physical reference;
                         target-active ancestry or carry transport;
                         a positive inverse or endpoint current;
                         any LRC row exclusion or decrement;
                         LRC(14), JC(2), DC(2), or G1.             (21)
```

The independent hostile audit rederived both factorization orientations
and the common `re` dimension without semisimplicity.  It checked the
prime-power Boolean specialization, the Fourier-convention boundary,
every Cayley control, both sharp `C_3` maps, and the `C_6` hostile.
Normal, optimized, and stored transcripts agree exactly; both declared
LF-normalized hashes match; the companion compiles without an optional
CAS, truth-bearing `assert`, or floating-point decision; and the
documentation gate passes.

**QED.**
