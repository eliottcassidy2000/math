---
id: THM-2465
title: "The G1 exclusion package: degree-four 2-jet Keller maps are pinned"
status: >
  PROVED (universal layer: no 2-jet Keller map has A != 0 with
  A x B == 0; B != 0 always; universal z-rationality -- for A != 0
  the generic fiber injects into the (x,y)-plane, so field degree =
  degree of a plane 0-cycle and the z-quadratic never doubles it; on
  THM-2633 excludes C4,V4,D4, leaving A4/S4; on both surviving branches
  there are no proper intermediate fields (A3 max in A4, S3 max in S4).
  The historical D4 root-field and matching quadratics remain valid
  conditional algebra but no longer define a live Keller lane.  Every
  quadratic Keller map is injective, so any G1
  witness has total degree >= 3) + PROVED (stratum exclusions,
  QQ-exact Groebner with cross-prime adversarial re-derivation:
  point-cap b == 0 face empty; line-cap doubly-degenerate face
  empty; conic-cap empty of Keller maps of ANY field degree in four
  boxes through (A=2,B<=2,C<=3); box (1,2,3) r2a stratum forces
  field degree exactly 1; on the affine point-cap face deg(f),deg(b)<=1,
  D3=-B1^4 Jac(f/B1^2,B2/B1): rank-two b forces f=0, rank-one
  varying projective direction is exactly removable by a quadratic target
  shear; after normalizing its z-affine top frame to (1,x,a), the remaining
  Keller equations integrate to an explicit polynomial inverse, so the
  whole varying-direction branch has field degree one; on the constant-
  projective face, rank-zero scale with nonconstant f is also invertible and
  rank-zero scale with constant f is exactly a suspension of a planar Keller
  pair, while rank-one scale-varying direction is empty; hence every affine
  point-cap Keller map here is an automorphism or an exact planar suspension;
  the (1,2,3) cusp engine needs deg B >= 3)
  + PROVED (transfer/hardness floor: affine point-cap field degree four
  exists iff planar Keller field degree four exists; field degrees 1 and 3
  are realized on both point- and line-cap strata by explicit elementary
  extensions of the THM-1310 wild map, and a broader unconditional degree-4
  exclusion on those strata implies the degree-4 case of the order-{1,3}
  conjecture) + PROVED (typed tower constraints: purity forces
  a nonempty Jelonek set and a nontrivial matching-quotient layer to ramify
  over it; S4 witnesses need odd Jelonek valuation of Delta_4, and A4
  witnesses have square Delta_4 with cyclic cubic layer.  The now-excluded
  D4 lane conditionally has a 1+2 matching algebra, an auxiliary quadratic,
  and an unramified split owner pair.  A
  THM-1310-as-resolvent witness must satisfy the five conditions
  (N1)-(N5)). SUPERSEDED GLOBAL VERDICT: THM-3438 gives an explicit
  z-quadratic field-degree-four Keller witness, so G1 is positive.  The
  exclusions and detection floors here remain valid on their stated strata;
  the z-affine order-{1,3} conjecture and JC(2) remain open.
source: kind-pasteur-2026-07-26-S134
depends_on:
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
related:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - THM-1310-conic-pair-fibers-and-design-equations
  - THM-1340-engine-trichotomy-zaffine-keller
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
  - MISTAKE-297
script: 04-computation/jacobian_g1_stratum_exclusions_kps_S134.py
output: 05-knowledge/results/jacobian_g1_stratum_exclusions_kps_S134.out
point_cap_script: 04-computation/jacobian_g1_point_cap_affine_dichotomy_thm2465.py
point_cap_output: 05-knowledge/results/jacobian_g1_point_cap_affine_dichotomy_thm2465.out
point_cap_script_sha256: 9e7a7f6d561da17166f79d53db4ddd04808a056c16b20bb40e9f8fd9721c385c
point_cap_output_sha256: c75aa39713dc2ac523c987f8cf5e5a0df0494b7fe5a5d979b21826e2c8d7aad3
companions: >
  jacobian_g1_degree_arithmetic_kps_S134, jacobian_g1_strata_kps_S134,
  jacobian_g1_line_strata_kps_S134,
  jacobian_g1_staircase_replication_kps_S134,
  jacobian_g1_conic_cap_hunt_kps_S134 (+ master out),
  jacobian_g1_conic_kernel_structure_kps_S134,
  jacobian_g1_point_cap_affine_dichotomy_thm2465,
  jacobian_g1_gb_wildcert_kps_S134 (GB wildcert boxes: PENDING,
  no verdict is claimed from them)
hash_basis: working-tree bytes (LF); per-file hashes in INDEX entry
---

# THM-2465 -- what a degree-four 2-jet Keller map must be

**PROVED LOCAL EXCLUSIONS + GLOBAL VERDICT SUPERSEDED/POSITIVE** as itemized
in the status. Produced by a
six-agent workflow (theory strata, degree arithmetic, staircase and
conic-cap hunts, adversarial verification, synthesis); every
load-bearing computational claim was independently re-derived over
QQ and fresh primes by the adversarial pass, with zero refutations.

> **GLOBAL G1 CORRECTION (THM-3438 / MISTAKE-396).**  The weighted-lift map
> `G` is polynomial, quadratic in `z`, has constant determinant `-6`, generic
> field degree four, and the collision `G(1,0,0)=G(-1,0,2)`.  It lies outside
> the boxes eliminated here, in the previously unsearched high-degree
> line-cap tail, and has geometric monodromy `S_4`.  Consequently the old sentence “G1 remains
> open” was false globally; every local exclusion below survives unchanged,
> and `G` does not touch the z-affine order-`{1,3}` conjecture.

> **MONODROMY CORRECTION (THM-2598 / MISTAKE-297).**  The original
> synthesis inherited THM-1375's unconditional self-normalizing rule and
> silently removed `D4`.  THM-1365 proves only the polynomial-deck case.
> That correction restored the provisional list `D4,A4,S4`.  THM-2633 now
> removes `D4` by a different point-stabilizer abelianization argument that
> never extends the deck map.  The live list is `A4,S4`; every
> no-intermediate-field assertion applies to both.  The graded exclusions and
> computations are unchanged, and the D4 paragraphs remain conditional
> field-theory controls.

## 1. Universal layer (PROVED)

For `F = A z^2 + B z + C` Keller: (U0) `B != 0` (else `D0 == 0`);
(U1) `A != 0` with `A x B == 0` is impossible (`B = rho A` makes
`det J` vanish along `z = -rho/2`); (U2) with `A != 0`, the frame
identity `z (A x B) = A x (P - C)` (THM-2455) makes `z` rational
over `C(P)(x,y)` on the generic fiber, which therefore injects into
the `(x,y)`-plane: **field degree = degree of the plane 0-cycle
`{g1 = 0} intersect {lambda = mu^2}`** -- the z-quadratic shape
never doubles the degree; (U3) on the `S4` and `A4` branches there is
no proper intermediate field (`S3` maximal in `S4`, `A3` in `A4`).
On the historical conditional `D4` branch the point stabilizer `C2` lies in an
order-four subgroup, so a quadratic intermediate layer exists and the
matching resolvent factors `1+2`; THM-2598 proves that its auxiliary
matching quadratic is distinct from that root-field intermediate.  A quartic
`D4` root field has only one proper intermediate field.  Since (U2) gives
`K_4=K(x,y)`, if neither coordinate generated `K_4`, then both `K(x)` and
`K(y)` would lie in that same quadratic intermediate, forcing
`K(x,y)` to lie there as well.  Hence at least one of `x,y` is a primitive
quartic coordinate in the `D4` model too; the old no-intermediate-field proof
was false, but its useful conditional coordinate conclusion survives by
uniqueness.  THM-2633 excludes this model from affine Keller realization.

The elementary quartic

```text
T^4-2,                    resolvent W(W^2+8)                 (U4)
```

is the sharp field-theory hostile: `Q(2^(1/4))` contains `Q(sqrt(2))`, while
the nontrivial resolvent factor gives the different quadratic
`Q(sqrt(-2))`.  It is not a Keller example; it proves that the two visible
quadratic sidecars cannot be identified formally.  Every quadratic Keller map is
injective (midpoint argument), so **any G1 witness has total degree
at least 3**; the staircase boxes have even BKK ceilings
`(d+2)^3 - d^3 = 8, 26, 56, 98` with the excess over the realized
degree dropping by exactly 2 per rank-<=1 point of the fiber system.

## 2. Stratum exclusions (PROVED, QQ-exact)

### 2.1 Affine point-cap dichotomy (PROVED)

Put a point-cap into constant target-linear normal form

```text
F=(B1 z+C1, B2 z+C2, f z^2+B3 z+C3),                         (PC)
```

and assume exactly `deg(f),deg(B1),deg(B2)<=1`, with `f!=0`.  No
degree hypothesis on `C1,C2,B3,C3` is needed for the first identity.  Let
`Jac(g,h)=g_x h_y-g_y h_x`; on the chart `B1!=0`, the coefficient `D3` of
`z^3` in `det JF` is

```text
D3=-B1^4 Jac(f/B1^2,B2/B1).                                  (PC1)
```

Indeed, with `b=(B1,B2)`,

```text
L=det(b_y,b),  M=det(b_x,b),  N=det(b_x,b_y),
D3=f_x L-f_y M+2fN,
X=L partial_x-M partial_y,
X(B1)=-N B1.
```

Expanding the rational Jacobian gives
`B1^4 Jac(f/B1^2,B2/B1)=-(X(f)+2fN)=-D3`.  Swapping `B1,B2`
covers the other nonzero chart.  Thus (PC1) is a polynomial identity after
clearing denominators, not a generic-point heuristic.

Let `rho=rank(dB1,dB2)`.

1. **Rank two is empty for a genuine point-cap.**  Affine source
   coordinates make `(B1,B2)=(x,y)`.  Writing `f=ax+by+c` gives

   ```text
   D3=ax+by+2c.
   ```

   Hence `D3=0` iff `a=b=c=0`, contradicting `f!=0`.

2. **Rank one with varying projective direction is exactly a target
   shear.**  Write

   ```text
   ell=ux+vy,   B1=a1 ell+c1,   B2=a2 ell+c2,
   delta=a1 c2-a2 c1.
   ```

   The projective direction `[B1:B2]` is nonconstant exactly when
   `delta!=0`.  Directly,

   ```text
   D3=delta (v f_x-u f_y).                                    (PC2)
   ```

   For affine `f`, equation (PC2) vanishes iff `f=m ell+n`: the gradient
   of `f` is parallel to the nonzero gradient of `ell` (the implication is
   immediate after choosing `ell` and one transverse affine coordinate).
   Since `delta!=0`, `(B1,B2)` is a basis of `span{ell,1}`.  The symmetric
   square change-of-basis from
   `(B1^2,B1B2,B2^2)` to `(ell^2,ell,1)` has determinant

   ```text
   delta^3 != 0.                                               (PC3)
   ```

   Therefore there is a unique homogeneous quadratic form `Q(S,T)` with
   `f=Q(B1,B2)`.  Conversely, `f/B1^2=Q(1,B2/B1)` makes (PC1) vanish,
   proving both directions.  The target map

   ```text
   (P1,P2,P3) -> (P1,P2,P3-Q(P1,P2))                          (PC4)
   ```

   is a polynomial automorphism, with inverse obtained by `+Q`; it preserves
   the Keller property, fibres, and field degree.  Its new `z^2` coefficient
   is `f-Q(B1,B2)=0`.  Thus this entire varying-direction face is genuinely
   Keller-equivalent to a `z`-affine map, although (PC4) may raise the
   `(x,y)` degrees of its lower coefficients.

   In fact the Keller equations close this normalized `z`-affine face
   completely.  A constant target change in the first two coordinates and
   an affine source change put the sheared map into

   ```text
   G=(P1,P2,P3)=(z+C1, xz+C2, a z+C3).                        (PC5)
   ```

   Define `h=C2-xC1` and `k=C3-aC1`.  Subtracting `x` times the first
   Jacobian row from the second and `a` times it from the third gives

   ```text
   det JG=(P1+h_x)(a_y P1+k_y)-h_y(a_x P1+k_x).               (PC6)
   ```

   The quadratic coefficient forces `a_y=0`, so `a=a(x)`.  The linear
   coefficient then gives `k_y=a'(x)h_y`, hence

   ```text
   k=a' h+b(x),             det JG=-h_y(a''h+b')=lambda in C*. (PC7)
   ```

   Both factors in the last product are units.  Thus `h_y=d in C*`, so
   `h=d y+e(x)`.  Its nonzero `y` coefficient forces `a''=0`, and then
   `b'=gamma in C*`.  For constants `alpha,beta,delta`, every such map has

   ```text
   a=alpha x+beta,    h=d y+e(x),
   k=alpha h+gamma x+delta,       lambda=-d gamma.             (PC8)
   ```

   It is a polynomial automorphism: from a target point `(P1,P2,P3)`,

   ```text
   x=(P3-beta P1-alpha P2-delta)/gamma,
   y=(P2-x P1-e(x))/d,
   z=P1-C1(x,y).                                               (PC9)
   ```

   Hence the entire rank-one varying-projective branch has field degree
   exactly one.  This conclusion uses the special top frame `(1,x,a)` and
   does not settle arbitrary `z`-affine degree-four maps.

3. **Constant direction: the rank-zero scale is exactly planar.**  If
   `delta=0`, `B1,B2` are proportional as affine polynomials wherever the
   projective direction is defined; then `L=M=N=0`, so `D3=0` imposes no
   condition on `f`.  There are nevertheless two different scale ranks.

   First suppose `rho=0` and `(B1,B2)!=(0,0)`.  A constant target change
   makes `(B1,B2)=(1,0)`.  Put

   ```text
   s=P1=z+C1,       G=B3-2fC1,
   R=C3-B3C1+fC1^2.                                      (PC10)
   ```

   The triangular source coordinate `z -> s` writes the map as

   ```text
   (s, C2, f s^2+G s+R),
   ```

   and its complete Keller system is

   ```text
   Jac(C2,f)=0,      Jac(C2,G)=0,      Jac(C2,R)=lambda in C*. (PC11)
   ```

   If the affine `f` is nonconstant, take `f=x`.  Then (PC11) forces

   ```text
   C2=d x+c,       G=G(x),       R=r y+E(x),       dr=lambda, (PC12)
   ```

   with `d,r` nonzero.  The inverse is explicit:

   ```text
   x=(P2-c)/d,
   y=(P3-x P1^2-G(x)P1-E(x))/r,
   z=P1-C1(x,y).                                              (PC13)
   ```

   If `f` is a nonzero constant, subtract `f P1^2` from the third target.
   Since `(C2,R)` is a planar Keller pair and `Jac(C2,G)=0`, THM-2230's
   Cheng--McKay--Wang centralizer theorem gives a unique `H` with
   `G=H(C2)`.  Subtracting `H(P2)P1` then leaves

   ```text
   (P1,C2,R),                                                   (PC14)
   ```

   which, after the triangular source coordinate `P1=z+C1`, is exactly the
   direct suspension of the planar Keller map `(C2,R)`.  Conversely every
   planar Keller map lifts to this genuine point-cap face by adding
   `f P1^2` to its suspended third target.  Hence field degree is preserved
   in both directions: **rank-zero constant-direction degree four is exactly
   the degree-four planar Jacobian problem**, not a new G1 species.

   The earlier minimal control

   ```text
   F0=(z,y,x+y z^2),             det JF0=-1,                  (PC15)
   ```

   lies in the nonconstant-`f` automorphism row; its inverse is
   `(a,b,c)->(c-ba^2,b,a)`.  No homogeneous quadratic `Q(B1,B2)` removes its
   cap, while the cubic target shear `P3->P3-P2 P1^2` does.  It therefore
   marks the failure of the quadratic shear only, not a G1 obstruction.

4. **The rank-one scale-varying constant direction is empty.**  Normalize
   `(B1,B2)=(x,0)`.  Direct expansion gives

   ```text
   det JF=(2fz+B3)(z C2_y+Jac(C1,C2))
           +x[z^2 Jac(C2,f)+z Jac(C2,B3)+Jac(C2,C3)].          (PC16)
   ```

   Write the nonzero affine cap as `f=ax+by+c`.  The quadratic equation is

   ```text
   b x C2_x+(a x+2b y+2c)C2_y=0.                              (PC17)
   ```

   If `b!=0`, the affine coordinate
   `Y=y+c/b+(a/b)x` changes its derivation to

   ```text
   b(x partial_x+2Y partial_Y).
   ```

   Every nonconstant monomial has positive weight, so its polynomial kernel
   consists only of constants.  This would make `C2` constant and
   `det JF=0`, impossible.  Hence `b=0`; (PC17) then forces `C2=C2(x)`.
   The constant equation in (PC16) factors as

   ```text
   C2'(x)[x C3_y-B3 C1_y]=lambda.                             (PC18)
   ```

   Thus `C2'=d in C*` and the bracket is `r=lambda/d in C*`.
   The linear equation is

   ```text
   x B3_y=2(ax+c)C1_y.                                        (PC19)
   ```

   At `x=0`, a nonzero `c` would make `C1_y(0,y)=0`, contradicting
   (PC18), so `c=0` and `a!=0`.  Equation (PC19) now gives
   `B3=2aC1+h(x)`, and integrating (PC18) in `y` gives

   ```text
   x C3=a C1^2+h(x)C1+r y+k(x).                               (PC20)
   ```

   Specializing (PC20) at `x=0`, with `u(y)=C1(0,y)`, would require

   ```text
   a u(y)^2+h(0)u(y)+r y+k(0)=0.                              (PC21)
   ```

   If `u` is nonconstant, the square has degree `2 deg u`, strictly larger
   than every other `u`-term and than `r y`; if `u` is constant, the nonzero
   linear term `r y` remains.  Both cases contradict (PC21).  Thus no Keller
   map exists in this last scale-varying branch.  The already-proved `b=0`
   exclusion handles `(B1,B2)=(0,0)`.

Consequently every varying-projective affine point-cap branch is an
automorphism, and the rank-zero constant-projective branch is exactly either
an automorphism or a planar suspension; the rank-one constant-direction
scale is empty.  Therefore **every affine point-cap Keller map in this
degree regime is polynomially equivalent either to an automorphism or to
the suspension of a planar Keller map**.  No arbitrary `z`-affine
field-degree-four map is excluded.

- Point-cap (`A = f v0`) with `b == 0`: empty of Keller maps
  (`det J = jac(C1, C2) (2 f z + B3)`).
- Line-cap doubly-degenerate face (`B3 == 0` and `j(Abar, Bbar) == 0`):
  empty (collapses into U1).
- Conic-cap `A`: empty of Keller maps of ANY field degree in the
  boxes `(A=C1, B<=1, C<=2)`, `(A=C2, B<=1, C<=2)`,
  `(A=C1, B<=2, C<=2)`, `(A=C1, B<=2, C<=3)` (normal forms
  `C1 = (x^2, xy, y^2)`, `C2 = (x^2, x, 1)` exhaust the cap modulo
  gauge); plus all 47 enumerated B-rays of the minimal quartic box
  and their gauge orbits, and the (2,4,4) kernel directions. The
  exact kernel law `ker D4|_{C1} = {hA + fA_x + gA_y} + {linear B}`
  is part of the record.
- Box (1,2,3): the `r2a` stratum forces field degree exactly 1 (the
  Euler-eigenvalue kill `D4 = 0 <=> B3 = 0`, then `C3` constant,
  then a polynomially z-solvable shape which is an automorphism);
  the cusp engine needs `deg B >= 3`.

## 3. Realized degrees and the hardness floor (PROVED)

`W1 = (F1, F2, F3 + F1^2)` (point-cap) and
`W2 = (F1 + F2 F3, F2 + F3^2, F3)` (line-cap), both built on the
THM-1310 wild map, are 2-jet Keller maps of field degree 3; degree 2
is Galois and hence impossible by the cited Campbell criterion.  The affine
point-cap classification in §2.1 proves that the `W1` mechanism is exhaustive
there: modulo polynomial source/target automorphisms, every potentially
non-automorphic map on that face is the direct suspension of a planar Keller
pair.  In particular,
an affine point-cap field-degree-four map exists **if and only if** a planar
field-degree-four Keller map exists.

Hence field degrees `{1, 3}` are realized on both residual strata, and **an
unconditional degree-4 exclusion there would prove the degree-4 case of the
z-affine order-{1,3} conjecture**.  THM-3438 settles global G1 positively on
a z-quadratic stratum but supplies no z-affine witness.  The residual arena
for the order conjecture is the point/line-cap z-affine quartic-engine problem plus the
named pockets with floors (r2b, r1, staircase branches, conic-cap
off-ray kernel families, deg >= 5 tails, and the pending GB boxes).

## 4. Tower constraints on any witness (PROVED)

Purity (Zariski-Nagata) plus `pi_1^et(A^3) = 1` re-proves that a
wild Keller map has nonempty Jelonek set `A_F`, and forces the
quartic field, its nontrivial matching quotient, and the Galois closure
to be unramified in codimension one off `A_F` while that nontrivial
quotient must ramify over some Jelonek component.  On the two live lanes the
quotient is cyclic cubic for `A4` and full-`S3` cubic for `S4` (THM-2598).
`S4` witnesses have odd
`v_D(Delta_4)` along a Jelonek component (the HYP-9027 shape), while
`A4` witnesses have square `Delta_4` with order-3 inertia.

For historical comparison, the now-excluded `D4` row has a further
conditional boundary invoice.  Its
generic deck involution always extends over the finite Zariski-main
normalization.  It can evade the polynomial-deck fixed-point exclusion only
by failing to preserve the open copy of `A^3`.  Because normalization
ramification is invariant and lies in the missing boundary, THM-2598 shows
that a `D4` witness must contain an **unramified missing divisor** exchanged
with an included divisor over the same target component.  At its generic
point the deck `C2` pair has opposite omitted/present owner bits; it is the
antipodal pair in the canonical square selected by the rational resolvent
matching.  This is an upstairs condition; even
set-theoretic equality of branch and Jelonek images would not detect it.

A witness carrying THM-1310's map as its resolvent
layer must satisfy: (N1) `G = S4`; (N2) `Delta_4 == -L o iota` mod
squares; (N3) its Jelonek set contains a copy of `{L = 0}`;
(N4) `v(Delta_4)` odd along it; (N5) its Galois closure contains
`K3(sqrt(-L))` as the degree-6 layer.

## 5. Detection floors and pendings (honest negatives)

Sampled-only (floors, not proofs): `r2b` (3/3 inconsistent), `r1`
(3/3 degenerate), staircase branches (8/8), conic-cap B off the 47
rays inside the 13/15-dimensional kernel families, `deg C >= 5`,
`deg B >= 5`, `deg A >= 3` caps, and the `(C2, B<=2, C<=2)` joint
box (watchdog kill, NO VERDICT). The GB wildcert boxes
(`(1,1,1)+fiber4`, `(1,2,2)+fiber4`, `(2,2,2)+alpha3`) are built
but unfinished: **no verdict is claimed from them**; the decisive
next computation is finishing `(1,2,2)+fiber4` -- `GB = [1]` there
would be the first stratum-free G1 exclusion in a box not killed by
the quadratic-injectivity argument (BKK ceiling 19 >= 4).

## 6. Reproduction

```bash
python 04-computation/jacobian_g1_stratum_exclusions_kps_S134.py   # 28/28 asserts
python 04-computation/jacobian_g1_degree_arithmetic_kps_S134.py
python 04-computation/jacobian_g1_conic_cap_hunt_kps_S134.py       # long
python 04-computation/jacobian_g1_gb_wildcert_kps_S134.py box122   # pending box
```

Outputs in 05-knowledge/results with matching names; the stored
outs are the workflow transcripts whose load-bearing claims were
adversarially re-derived cross-prime in-session.
