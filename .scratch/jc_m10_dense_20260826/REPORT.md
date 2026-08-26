# Exact-`M=10` reduced `(2,3)` seam: scratch frontier report

**Status:** `FINITE-EXACT / TRACKED SCRATCH`, 2026-08-26.  This directory is
deliberately outside canon.  The generic exclusion is a theorem candidate
relative to the already-proved THM-4147 transport, not canon until
independently audited and promoted.  `JC(2)`, seam entry, and every
coefficient wall remain open.

## Inheritance

- Closest mechanism: THM-4147, `generic-exact-weight-nine-planar-jacobian-
  monodromy-exclusion.md`, especially its source-critical length, prime
  carrier, and finite/full permutation gates.
- Canonical hostile: THM-4045, `live-two-three-max-seven-hidden-elliptic-
  tail-no-go.md`; a highest face can miss a lower positive-genus tail.
- Corrected near miss: MISTAKE-519; search ambient equations and universal
  quantifiers before calling a specialization new or open.
- Least-used sidecar: the endomorphism algebra of the normalized face
  Jacobian together with labelled attachment differences on the hidden tail.
- Meta-patterns used: `Search the statement before the method`, `Certify
  before projection`, and `Find the hidden second coordinate`.

## Exact object

On the inherited `b=d=0` seam let

```text
u=[p^5]H,        v=[p^2 y^2]H,        zeta=[y^3]H.
```

Assume

```text
u*v*zeta*(u+v) != 0.                                  (1)
```

The complete monomial universe through weight ten is

```text
p,p^2,p^3,y^2,p^2y,p^4,py^2,p^3y,y^3,p^5,p^2y^2.
```

For `F_Q=(s^2-p)(1-QH)-Qs^2/2`, exact lower-hull enumeration gives only

```text
M: h=(i+2j-2)/10,
T: h=(i+j-2)/6.                                       (2)
```

This stays true for all `2^8=256` choices of lower support once the three
terms in `(1)` are retained.  The global polygon is

```text
(0,1),(2,0),(5,3),(4,4),(0,6),
(2Area,boundary,g)=(36,12,13).                         (3)
```

Its complete toric packet is

```text
e=(9,9,6,2,2,2,1),       n=31,       defect=24=2g-2. (4)
```

The main face factors as

```text
(S^2-P)(1-uP^5-vS^2P^4).                              (5)
```

The nonrational factor normalizes via `Y=SP^2` to

```text
vY^2=1-uP^5.                                          (6)
```

It is a smooth genus-two curve.  Its order-five automorphism acts on
holomorphic differentials with characters `zeta_5,zeta_5^2`, so its
Jacobian contains `Q(zeta_5)` in its rational endomorphism algebra.
That Jacobian is simple: if it were isogenous to nonisogenous elliptic
curves, a degree-four field could not inject into either endomorphism
factor; if it were isogenous to `E^2`, then `M_2(Q)` has maximal
commutative degree two, while for CM `E`, adjoining the imaginary-quadratic
center to `Q(zeta_5)` would give a commutative algebra of dimension eight
inside a degree-two central simple algebra, whose commutative maximum has
dimension four.  Hence it has no elliptic quotient and in particular no
nonzero Hom to the target `j=0` curve.

The hidden side face is nevertheless

```text
1-zeta S^3P^3-vS^2P^4=0.
```

With `T=SP,W=SP^2`, this is

```text
vW^2=1-zeta T^3,                                      (7)
```

an elliptic curve isomorphic to `E_0:y^2=x^3+1`.  The length-two internal
edge attaches it at `(0,1),(0,-1)`.  Their difference is nonzero
three-torsion: doubling `(0,1)` gives `(0,-1)`.  Thus the endomorphism
`1-zeta_3` identifies the attachments.  The highest-face Hom test and its
attachment repair both fail sharply, and any compatible specialized map has
degree divisible by three.

## Dense critical-open exclusion candidate

Use the lawful rational control

```text
Delta=1, K=2848/45-(7/6)Delta=5591/90,
Phi=2, Theta=5, eta=7, zeta=11, u=13, v=17.            (8)
```

For the exact source pair

```text
A=(-sp+t^2 H_s)/p,
C_0=s^2+2t^2H_p,
B=(C_0+sA)/t^2,
t=p-s^2,
```

the two independent eliminations agree:

```text
Res_s(A,B)=p^6 R_25(p),
Res_s(A,C_0)=p^8 R_25(p),                              (9)
```

where `R_25` is squarefree and

```text
R_25(0)=-189675421056/5,
[p^25]R_25=2731549392000000000.                       (10)
```

There is no `p=0`, `t=0`, or source-infinity loss.  The disjoint `(X,T)`
projection is

```text
Res_X(G_X/T,G_T)=T^72(6T+1)^2 Q_25(T),                (11)
```

with squarefree `Q_25`, nonzero zero and `-1/6` endpoints.  Restoring the
two universal Morse points at each of `T=0,-1/6` gives

```text
L=25+2+2=29.                                          (12)
```

The edge joining the base term to `y^3` is the same prime separable cubic
carrier as in THM-4147:

```text
q-1/2=K W^2+zeta W^3.                                 (13)
```

Therefore the THM-4147 response arithmetic becomes

```text
finite: (n,beta)=(25,3),
        2n-L-1+beta=23 < n-1=24;
full:   n=31,
        2(n-L)=4 < defect=24.                          (14)
```

Both alternatives are strictly impossible on the critical-open chamber,
provided an independent audit verifies the complete packet, source
critical scheme, carrier transport hypotheses, and the use of THM-4147's
already-proved permutation lemmas.  The rational control proves that this
open chamber is nonempty; it is not itself a Keller pair.

## Connection contract

```text
source:     complete exact-M=10 H and its unsaturated critical ideal
target:     lower Newton component inventory + labelled monodromy response
map:        Q-adic lower model; independently, source ideal -> resultants
preserved:  genus, puncture indices, prime carrier degree, critical length
destroyed:  source coordinates in each resultant; lower coefficients in hull
sidecars:   second critical projection, p/t/infinity checks, universal fibres,
            carrier residue field, attachment labels
cheap test: independent expansion and the two strict inequalities (14)
```

## Ranked next fronts

1. **Audit/promote the critical-open chamber.**  Re-expand `F_Q` without the
   scratch support helper, recompute both resultants from source coordinates,
   audit every outer edge and THM-4147 dependency, then promote a narrowly
   scoped exact-`M=10` theorem.
2. **Close `zeta=0`, `uv!=0` off its edge walls.**  Exhaustion shows that its
   only optional side face is
   `1-(SP)^2(K+Theta P+vP^2)=0`, whose normalization is
   `Y^2=K+Theta P+vP^2`, hence rational.  The cyclotomic genus-two component
   is then the only abelian component.  A complete lower-model audit should
   give a no-Hom exclusion away from `u+v=0` and
   `Theta^2-4Kv=0`.
3. **Stratify critical walls in the `zeta!=0` chamber.**  Full response still
   closes whenever `L>=20`, but the elementary finite-carrier inequality
   needs `L=29`.  Reuse the source-resultant endpoint cascade and carrier-
   orbit refinements rather than treating every discriminant zero alike.
4. **Exploit the hidden-tail arithmetic on surviving walls.**  The complete
   specialization forces an Eisenstein-norm cover degree divisible by three
   and generic ramification `2g-2=24`.  Combine this with critical-value and
   nonproperness allocations; it is a real sidecar, not a contradiction by
   itself.
5. **Then handle support contractions.**  Separate `u=0`, `v=0`, and the
   repeated ten-node wall `u+v=0`.  Do not call singleton highest faces
   closed without a complete lower model: THM-4045 is the hostile.
6. **Keep entry separate.**  Every statement here begins inside the inherited
   reduced `(2,3)`, `b=d=0` seam.  No result forces an arbitrary planar
   Keller counterexample into it.

## Reproduction

```bash
python3 -B .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.py \
  | diff -u \
      .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.out -
python3 -B -O .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.py \
  | diff -u \
      .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.out -
sha256sum \
  .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.py \
  .scratch/jc_m10_dense_20260826/jc_m10_lower_hull_probe.out
```

Current hashes:

```text
449a0f3db9e305225744e7629ad1be0aaef7c7aedba97f27868a8b6f0fe87619  script
82470e40ba60cc6cc7dc165c36ad87d96b0d8413975595f8e1b4b32d6584e72f  output
```
