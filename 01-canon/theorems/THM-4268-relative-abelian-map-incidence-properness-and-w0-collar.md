---
id: THM-4268
title: "Relative abelian-map incidence properness and W=0 collar"
status: >
  PROVED RELATIVE TO THM-4230/4260 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. On the complete exact-M=12 gate, the fixed-degree 34/42 map
  incidence which collapses the twelve contacts is proper over the top-
  coefficient base; after fixing the common target value it is finite. Its
  W=0 fibre is empty by THM-4260, so its image is closed and misses the whole
  gate-interior W=0 divisor. The complementary Zariski open is an actual,
  though non-effective, W=0 collar with no inherited exact-M=12 Keller
  candidate; every formal trait centred on W=0 is excluded as well. No
  explicit collar equation, boundary-wall crossing, seam entry, JC(2), or
  DC(2) follows.
source: codex-longer-frontiers-20260827
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4260-w0-canonical-node-reciprocal-denominator-attachment-exclusion
related:
  - THM-4255-specialization-kernel-and-transverse-hasse-jet-repair
  - THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport
  - THM-4264-w0-visible-incidence-two-edge-attachment-observer
  - THM-4265-w0-reduced-wall-factor-and-transverse-jacobian-reduction
external: >
  The graph form and local finite-presentation statement for Grothendieck
  relative-Hom are recorded in Stacks Project Tags 0D1B/0D1C.
  Bosch--Lutkebohmert--Raynaud, Neron Models (1990),
  Propositions 1.2.8 and 1.4.4, supply that an abelian scheme over a DVR is
  the Neron model of its generic fibre and that generic maps from smooth
  sources extend uniquely. These inputs prove the valuative step; none
  supplies the exact M=12 face or the W=0 emptiness.
script: 04-computation/jc23_relative_incidence_collar_thm4268.py
output: 05-knowledge/results/jc23_relative_incidence_collar_thm4268.out
script_sha256: 470d056ac6832be6b8e7c084e93b5ff2b680b8dc94715dd3948a98c0de017a15
output_sha256: 7e6f1371b2b3e1162748b84a3747f46017e79c1e3b8da7f45816f7054f105d35
independent_script: 04-computation/jc23_relative_incidence_collar_independent_audit_thm4268.py
independent_output: 05-knowledge/results/jc23_relative_incidence_collar_independent_audit_thm4268.out
independent_script_sha256: a7914d34e899f2033fdde11487869fcd4de6e318e8755da1932b8a4cf2ad237c
independent_output_sha256: ca0c4abb1c61672be1e4c8f8567df26a5c2cdeb3ae479154c869ddb9b6700d1c
hash_basis: raw LF bytes
audit: >
  PASS. The primary symbolic companion reconstructs the face derivatives,
  critical identity, four edge gates, rank-twelve contact divisor,
  homogeneous V4 wall quotient, and nonproper empty-fibre hostile. A
  standard-library path independently executes 38,880 critical/contact and
  360 Kummer checks. The graph-Hilbert/Neron argument below is the
  load-bearing proof of properness; the scripts audit its exact face input.
---

# THM-4268 -- relative abelian-map properness and the `W=0` collar

**PROVED RELATIVE TO THM-4230/4260 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. THE COLLAR EXISTS BUT IS NOT EFFECTIVELY EQUATED. `JC(2)` REMAINS
OPEN.**

## 1. Statement and inheritance

Retain

```text
D=W^2-4UZ,                    Lambda=U+W+Z,
B=Spec C[U,W,Z,(UZD Lambda)^(-1)].                    (1)
```

Let `Cbar -> B` be the simultaneous smooth projective genus-seven face
constructed in Section 2, let `A_att <= Cbar` be its twelve-contact divisor,
and put

```text
E=E_0 x B,                    E_0:Y^2=X^3+1.
```

the marked special target used in THM-4230's definition of `H_0`. The proof
only needs an abelian scheme. Retaining an unnormalized `j=0` Weierstrass
twist gives the same argument after a finite etale Kummer trivialization;
finiteness, properness, and geometric emptiness descend.

For `d in {34,42}`, let `I_d` parametrize maps

```text
phi:Cbar_b -> E_b,             deg(phi)=d,             (2)
```

whose restriction to `A_att,b` is constant. Let `I_d^0` require that common
value to be the zero section.

> **Theorem.** The incidences `I_d -> B` are proper and the based incidences
> `I_d^0 -> B` are finite. Translation gives
>
> ```text
> I_d ~= E x_B I_d^0.                                  (3)
> ```
>
> On `B_0=V(W)<=B`,
>
> ```text
> I_34 x_B B_0=I_42 x_B B_0=empty.                    (4)
> ```
>
> Hence, for the finite scheme-theoretic images
>
> ```text
> Z_d=image(I_d^0 -> B),
> Omega=B \ (Z_34 union Z_42),                         (5)
> ```
>
> `Omega` is a Zariski-open neighbourhood of the whole gate-interior
> divisor `B_0`. No inherited exact-`M=12` Keller candidate has top
> coefficient point in `Omega`.

The conclusion is uniform in every lower coefficient allowed by THM-4230.
It does not give an analytic bound `|W|<epsilon(U,Z)` or a single uniform
epsilon on the noncompact base.

The inheritance board is:

- closest proved mechanism: THM-4260's complete intrinsic `W=0` exclusion;
- hostile: `Spec R[1/W] -> Spec R`, with empty special and nonempty generic
  fibre;
- corrected near miss: `partial_W d_ell(t)` is ill-typed because the reduced
  denominator and Hom basis were built only after setting `W=0`; and
- least-used sidecar: properness of fixed-degree maps to an abelian target,
  which retains the whole graph rather than one specialized coefficient.

## 2. A simultaneous face and unsplit contact divisor

On the torus the face is

```text
F=1-UP^6-WS^2P^5-ZS^4P^4.                            (6)
```

Its derivatives and an exact syzygy are

```text
F_S=-2SP^4(WP+2ZS^2),
F_P=-P^3(6UP^2+5WS^2P+4ZS^4),
2Z pbr-(4ZS^2+3WP)sbr=-3DP^2.                        (7)
```

Thus `D!=0` excludes torus critical points. The four boundary edge schemes
of THM-4230 are

```text
X-1,  1-ZX^4,  (X-1)(UX^2+WX+Z),  U-X^6,            (8)
```

with relevant discriminants/resultant

```text
D,  -256Z^3,  46656U^5,
Res_X(X-1,UX^2+WX+Z)=Lambda.                          (9)
```

Take one fixed regular subdivision of the normal fan of THM-4230's face
polygon and close `(6)` in the resulting smooth projective toric surface
over `B`. Equations `(7)--(9)` say the closure is fibrewise nondegenerate and
meets every toric boundary stratum transversely. The relative Jacobian
criterion therefore makes that closure `Cbar -> B` smooth and projective.
Its fibres have the constant Pick genus seven. This supplies the simultaneous
model rather than assuming that fibrewise normalizations commute with base
change.

The rational component meets it at `P=S^2`. Substitution gives

```text
1-Lambda S^12=0,              derivative=-12Lambda S^11. (10)
```

Hence `A_att -> B` is finite etale of rank twelve. It need not split over
`B`: saying `phi|A_att,T` is constant means it factors through `T`. After an
etale cover this is equality of twelve evaluations, so no radical or node-
orbit choice is silently discarded.

## 3. Relative map scheme and closed incidence

Choose relatively ample `L_C` on `Cbar` and the degree-one zero-section
bundle `L_E` on `E`. A graph of a degree-`d` map has a fixed genus-seven
Hilbert polynomial and polarization degree

```text
deg(L_C|Cbar_b)+d.                                    (11)
```

The graph locus is open in that fixed projective relative Hilbert space.
Thus fixed-degree maps form a separated finite-type scheme over `B`. The
conditions `phi|A_att=0` and “all evaluations agree” are closed: etale-
locally they are intersections of equalizers into the separated scheme `E`,
and closedness descends. This constructs `I_d^0` and `I_d`.

For an unbased incidence, the common value is a unique section `a` of `E`.
Translation by `-a` is based and translating back is inverse, proving `(3)`.
The Hilbert construction supplies finite type only; its graph locus is open
and is **not** being mistaken for a proper Hilbert component.

## 4. The abelian valuative closure

Let `R` be a DVR with fraction field `K` and a map `Spec R -> B`. Then
`Cbar_R/R` is smooth projective and `E_R/R` is an abelian scheme. A generic
based incidence is a map

```text
phi_K:Cbar_K -> E_K.                                  (12)
```

The abelian scheme is the Neron model of its generic fibre. Its mapping
property extends `(12)` uniquely to

```text
phi_R:Cbar_R -> E_R.                                  (13)
```

The relative line bundle `phi_R^*L_E` has constant fibre degree, so the
special map still has degree `d`. Its restrictions to `A_att,R` and the zero
section agree generically; flatness of `A_att,R` and separatedness of `E_R`
make them agree everywhere. Thus `(13)` stays in `I_d^0`. The valuative
criterion proves `I_d^0 -> B` proper, and `(3)` proves `I_d -> B` proper.

Every geometric fibre of `I_d^0` is finite. Over an algebraic closure choose
one attachment as Abel--Jacobi basepoint. A based map determines

```text
J(Cbar_b) -> E_b,                                     (14)
```

and its degree is a positive-definite Rosati form on the finitely generated
free Hom lattice. A fixed value `d` contains finitely many lattice vectors.
Hence `I_d^0` is quasi-finite; proper plus quasi-finite makes it finite.

The abelian target is load-bearing. The generic degree-one map
`[X:Y] |-> [pi X:Y]` on `P^1_R` acquires a base point and a vertical rational
component at `pi=0`; the corresponding graph locus is not proper. Neron
extension forbids exactly that escape here.

## 5. THM-4260 empties the special incidence

On `B_0=V(W)`, the gate is

```text
U*Z*(U+Z)!=0,                  D=-4UZ,                (15)
```

and every geometric source fibre scales to `C_0:x^6+y^4=1`. THM-4241 gives
the full integral Hom lattice; THM-4249's necessary projection and THM-4253's
last profile deletion reduce every degree-`34/42` collapsed-contact map to
THM-4260's exhaustive 280 classes. THM-4260 audits both node orbits and
excludes every class. Translation does not change degree or equality.
Therefore its conclusion `S_34=S_42=empty` identifies with the full geometric
fibres of the present functor, not merely a selected denominator table.

If `I_d^0 x_B B_0` were nonempty, finite type over `C` would give a closed
complex point, contradicting that exhaustive theorem. This proves `(4)`.
Since `I_d^0 -> B` is finite, `Z_d` in `(5)` is closed and disjoint from
`B_0`; hence `Omega` is the claimed collar.

THM-4230's response interface says a hypothetical nonautomorphic exact-
`M=12` Keller pair supplies a degree-`34` or degree-`42` component map which
collapses all twelve contacts. Thus no inherited candidate lies above
`Omega`.

For any trait centred on `B_0`, a generic incidence would extend by Section
4 and contradict `(4)`. Equivalently,

```text
(I_34 union I_42) x_B completion_(B_0)(B)=empty.       (16)
```

This controls all transverse orders without constructing `partial_W d_ell`.

## 6. Homogeneous wall coordinate

The same audit clarifies, but does not cross, the three `W=0` walls.
THM-4260's coordinate

```text
Z/U=((t^2-1)/(2t))^2                                  (17)
```

homogenizes for `t=T/R` to

```text
[T:R] |-> [U:Z]=[4T^2R^2:(T^2-R^2)^2],
U+Z=(T^2+R^2)^2.                                     (18)
```

This degree-four map is the Klein-four quotient by
`{t,-t,1/t,-1/t}`. Its three branch fibres, each consisting of two points of
ramification index two, are

```text
U=0: t=0,infinity;   Z=0: t=+1,-1;   U+Z=0: t=+i,-i. (19)
```

Thus `t`, `t^2-1`, and `t^2+1` have different geometric owners. Moreover,
homogenizing `(10)` as

```text
B_att^12=Lambda A_att^12                              (20)
```

shows that `U=0` or `Z=0` with the other endpoint nonzero still has twelve
reduced contacts. At `U+Z=0`, all twelve form one length-twelve contact at
`S=infinity`. This is only a contact-divisor statement; source normalization
and Keller exclusion on all three walls remain open.

## 7. Firewall, survivor, and reproduction

The nonproper family `WX-1=0` has empty `W=0` fibre and nonempty generic
fibre. Empty specialization alone therefore proves nothing. The legal map is

```text
source:       fixed-degree graph plus all contact evaluations,
target:       top-coefficient base B,
map:          proper incidence projection,
preserved:    existence, degree 34/42, twelve-fold equality,
lost:         individual map and translation,
sidecar:      smooth source, abelian target, etale contacts.          (21)
```

No coefficientwise vanishing, Cartier commutation, or derivative of a
specialized denominator occurs. The method does not repair the external
p-adic-zeta density proofs in THM-4255/4263 because no analogous proper
finite-type Cartier-incidence carrier has been constructed there.

The strongest survivor is

```text
PROVED: a non-effective Zariski/formal W=0 collar in the exact-M=12 gate;
OPEN:   equations for its complement, U=0/Z=0/U+Z=0 walls, M=12 seam entry,
        JC(2), and DC(2).                              (22)
```

Replay with

```bash
python3 -B 04-computation/jc23_relative_incidence_collar_thm4268.py
python3 -B -O 04-computation/jc23_relative_incidence_collar_thm4268.py
python3 -B 04-computation/jc23_relative_incidence_collar_independent_audit_thm4268.py
python3 -B -O 04-computation/jc23_relative_incidence_collar_independent_audit_thm4268.py
```

The scripts audit Sections 2 and 6 and the hostile. Sections 3--5 are the
geometric proof; the scripts do not claim to mechanize the cited standard
representability and Neron theorems. Normal and optimized outputs byte-match
the frozen files. **QED.**
