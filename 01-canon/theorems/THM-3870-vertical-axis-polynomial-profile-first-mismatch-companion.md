---
id: THM-3870
title: "Vertical-axis polynomial profiles force a projective companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The vertical line
  A=0 is a component of the
  depressed-cubic branch exactly for b=C^2/6+A Q(A,C).  For every polynomial
  Q, comparison with the unique formal binomial profile at the first A-adic
  mismatch produces a distinct reduced companion whose affine normalization
  misses at least two projective points.  The only degree resonance has a
  nonzero perfect-square leading coefficient.  Reverse polynomial graphs
  A=h(C) are independently classified as the vertical line or invertible
  affine graphs over A.
source: jc_sparse_direct_search / vertical first-mismatch completion of THM-3866, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit rederived the
  binomial indexing and exact quadratic response, including N=0; checked all
  three degree regimes and the noncancellation of the resonant square;
  verified that a factor selected at finite A and infinite C has positive
  C-degree and therefore dominates the A-line; and separately reconstructed
  the reverse-graph UFD argument, including its unit boundary.  The
  assertion-free companion checks
  the vertical classification, the universal nonlinear response, the exact
  truncation residuals for N=0,...,5, and every degree regime
  d=0,N+1,N+2,N+3,N+4, including the resonant square.  It also checks the
  nonvertical reverse-graph affine boundary.  Normal and optimized runs
  byte-match the frozen 206-gate transcript.  A structurally separate
  862-gate audit additionally tests A-dependent generic leading terms,
  arbitrary shared-q_N cancellation, actual-degree factorization with
  multiplicities, a nonreduced residual packet, dominance, and the two
  distinct normalization endpoints.  Its normal, optimized and frozen
  transcripts also byte-match.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3863-finite-binomial-hensel-peels-force-projective-companion-contact
  - THM-3866-all-polynomial-graph-branches-force-projective-companion
script: 04-computation/jc2_vertical_axis_first_mismatch_thm3870.py
output: 05-knowledge/results/jc2_vertical_axis_first_mismatch_thm3870.out
script_sha256: d3ff51e8021590f30ac453d8e0b37bfe455047091d16ae14c214d9d1003274a1
output_sha256: fbeb2f787e93323d1ba1e3e13a8cdc6dce9ae4c3a22162fa52462394c6b55acd
semantic_sha256: 9987de8723700c1944d67ebc50631afdab3952db4cc2bf8108f6ee5b562b90b0
independent_audit_script: 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py
independent_audit_output: 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
independent_audit_script_sha256: 0853470843c68d302466016e8bbc1c444d5497c75ca10097a05cb37df11ed401
independent_audit_output_sha256: 0a9c9ea793027f79e578271b5f53f396c705161fffe90c6f74cd11a8f426410a
independent_audit_semantic_sha256: 48aa2b48290aa53595bf9bbea0f0a7015fde2ea13515543cf62b6a941b857079
hash_basis: raw LF bytes
---

# THM-3870 -- the vertical axis has no polynomial escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an algebraically closed field `k` of
characteristic zero.  For `b=b(A,C) in k[A,C]`, put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                         (1)
```

The vertical line `L=V(A)` is a component of `V(Delta_b)` if and only if

```text
b=C^2/6+A Q(A,C),                    Q in k[A,C].               (2)
```

For every polynomial `Q` in `(2)`, the reduced branch packet contains an
irreducible component `Xi` different from `L` such that the affine
normalization of `Xi` omits at least two distinct points of its smooth
projective normalization.  One omitted point lies over

```text
A=0,                         C=infinity,                         (3)
```

and a second lies over `A=infinity`.  There is no degree bound on `Q`.

There is also an exact reverse-graph classification.  If a component is a
polynomial graph `A=h(C)`, then either `h=0`, or `h` is affine of nonzero
slope and the same component is an invertible affine graph `C=g(A)`.  Thus,
combining with proved THM-3866, every polynomial graph over either coordinate
axis in this carrier is closed by a projective-companion obstruction.
Genuinely non-graph polynomial `A1` branches remain open.

## 1. Exact vertical classification

Specializing `(1)` gives

```text
Delta_b(0,C)=9C^2-54b(0,C).                                    (4)
```

Consequently `A | Delta_b` if and only if `b(0,C)=C^2/6`, which is
equivalent to the unique form `(2)`.  Direct substitution records the full
first residual:

```text
Delta_(C^2/6+A Q)=A V_Q,

V_Q=-C^3-54Q-(3/4)AC^4-54ACQ-9A^2C^2Q-27A^3Q^2.               (5)
```

This identity is a useful low-order control, but repeated vertical
multiplicity requires the all-depth comparison below.

## 2. The unique formal vertical profile

The cusp identity is

```text
A^2 Delta_b=27(P^3-u^2),
P=1+(2/3)AC,                    u=1+AC+A^2b.                     (6)
```

In `k[C][[A]]`, the unique solution of `Delta_b=0` with `u=1 mod A` is

```text
b_*={(1+(2/3)AC)^(3/2)-1-AC}/A^2
   =sum_(j>=0) beta_j A^j C^(j+2),                              (7)

beta_j=binom(3/2,j+2)(2/3)^(j+2) !=0.                          (8)
```

Its constant term is `beta_0 C^2=C^2/6`.  Therefore the unique formal
vertical quotient is

```text
Q_*=(b_*-C^2/6)/A
   =sum_(j>=0) beta_(j+1) A^j C^(j+3).                          (9)
```

All coefficients in `(8)` are nonzero, so `Q_*` is not a polynomial.  Every
polynomial `Q` in `(2)` consequently has a unique first mismatch `N>=0`:

```text
Q=Q_<(N)+A^N T,
Q_<(N)=sum_(j=0)^(N-1) beta_(j+1)A^j C^(j+3),
t(C)=T(0,C) != beta_(N+1)C^(N+3).                              (10)
```

The empty sum is used at `N=0`.  The polynomial `T` may be zero; this is the
exact finite-truncation boundary.  Equality with the full formal profile is
not a polynomial boundary and therefore never evades `(10)`.

## 3. Truncation residual and exact nonlinear response

Put

```text
b_N=sum_(j=0)^N beta_j A^j C^(j+2),
u_N=1+AC+A^2b_N.                                                 (11)
```

The formal tail in `(7)` starts at order `A^(N+1)`.  Since
`Delta_(b_*)=0`, the exact quadratic response in the profile variable gives

```text
Delta_(b_N)=A^(N+1)R_N,                                        (12)

R_N(0,C)=54 beta_(N+1)C^(N+3).                                 (13)
```

For clarity, no completed-to-polynomial descent is hidden here.  The
polynomial `Delta_(b_N)` has its first `N+1` coefficients in `A` equal to
those of the zero series `Delta_(b_*)`; hence it is divisible by
`A^(N+1)` already in `k[A,C]`.  Moreover

```text
deg_C b_N=N+2,          lc_C(b_N)=beta_N A^N,
deg_C u_N=N+2,          lc_C(u_N)=beta_N A^(N+2).               (14)
```

The unique top term of `Delta_(b_N)` is the square term in `(1)`.  Thus

```text
deg_C R_N=2N+4,
lc_C(R_N)=-27 beta_N^2 A^(N+1).                                (15)
```

These formulas include `N=0`; there is no exceptional base case.

The exact response identity, valid for every perturbation `delta b`, is

```text
Delta_(b_N+delta b)-Delta_(b_N)
=-54u_N delta b-27A^2(delta b)^2.                               (16)
```

Equation `(10)` means `b=b_N+A^(N+1)T`.  Dividing `(16)` by
`A^(N+1)` gives

```text
Delta_b=A^(N+1)S,

S=R_N-54u_NT-27A^(N+3)T^2,                                   (17)

S(0,C)=54[beta_(N+1)C^(N+3)-t(C)] !=0.                        (18)
```

In particular the multiplicity of the marked vertical component is exactly
`N+1`, and no component of `S` is vertical.

## 4. All degree regimes and the resonant square

If `T=0`, equations `(13)` and `(15)` already give

```text
deg_C S=2N+4>N+3=deg_C S(0,C).                                 (19)
```

Now suppose `T!=0` and put

```text
d=deg_C T,                         tau(A)=lc_C(T).               (20)
```

The special degree is bounded by

```text
m:=deg_C S(0,C) <= max(N+3,d).                                  (21)
```

There are three exhaustive cases.

### 4.1. Lower degree d<N+2

The three degrees in `(17)` are respectively

```text
2N+4,                 at most N+2+d<=2N+3,
                       and 2d<=2N+2.                            (22)
```

Hence

```text
deg_C S=2N+4>max(N+3,d)>=m.                                    (23)
```

### 4.2. Higher degree d>N+2

The last term has unique degree `2d`; the other degrees are `2N+4` and at
most `N+2+d<2d`.  Therefore

```text
deg_C S=2d>max(N+3,d)>=m.                                      (24)
```

### 4.3. Resonance d=N+2

All three terms can have degree `2N+4`.  Their exact combined leading
coefficient is

```text
-27A^(N+1)[beta_N^2+2A beta_N tau+A^2tau^2]
=-27A^(N+1)(beta_N+A tau)^2.                                  (25)
```

The bracket has nonzero constant term `beta_N`, so it cannot vanish as a
polynomial in `A`.  Consequently

```text
deg_C S=2N+4>N+3>=m.                                           (26)
```

Together `(19)` and `(23)--(26)` give the strict finite-base drop

```text
D:=deg_C S > deg_C S(0,C)=m                                   (27)
```

for every polynomial profile and every first-mismatch boundary.

## 5. Multiplicity-safe projective companion

Homogenize `S` to its generic `C`-degree `D`, using `[C:Z]` on `P1_C`.
If `S_0^h` is the degree-`m` homogenization of the nonzero polynomial
`S(0,C)`, then

```text
mathcal S(0;C,Z)=S_0^h(C,Z)Z^(D-m).                             (28)
```

Since `S_0^h(1,0)` is its nonzero leading coefficient, `(28)` places the
total projective residual through

```text
P_0=(A=0,[C:Z]=[1:0]).                                         (29)
```

Factor with all multiplicities retained:

```text
S=c product_i G_i^(e_i),                                       (30)
```

where the `G_i` are distinct irreducibles.  Homogenizing each factor at its
actual `C`-degree shows that at least one reduced component
`Xi=V(G_i)` has projective closure through `P_0`.  It is not vertical by
`(18)`, and it is not a component contained in `Z=0`, because `G_i` is an
affine polynomial homogenized at its actual degree.  In fact `G_i` has
positive `C`-degree: a degree-zero irreducible factor vanishing at `A=0`
would be `A` up to a unit, already excluded by `(18)`.  Hence `Xi` dominates
the `A`-line and is distinct from `L=V(A)`.

Close `Xi` in `P1_A x P1_C` and normalize.  A point `p_0` above `(29)` is
absent from the affine normalization because `C/Z` has a pole there.  The
nonconstant projective morphism to `P1_A` is surjective, so there is also a
point `p_infinity` above `A=infinity`; it is absent because `A` has a pole.
Their `A`-images are zero and infinity, respectively, so they are distinct.
Repeated factors in `(30)` change only contact multiplicity, not the chosen
reduced component or either deleted place.

## 6. Reverse polynomial graphs

Suppose now that `V(Delta_b)` contains a graph

```text
A=h(C),                         h in k[C].                       (31)
```

Set `B(C)=b(h(C),C)`.  In the UFD `k[C]`, the restriction of `(6)` gives
`P^3=u^2`.  Absorbing the constant units gives a polynomial `z(C)` with

```text
P=z^2,                         u=z^3.                            (32)
```

Thus

```text
hC=(3/2)(z-1)(z+1),
h^2B=(1/2)(z-1)^2(2z+1).                                      (33)
```

If `h=0`, this is the vertical case already proved.  Assume `h!=0`; then
`z-1` is not the zero polynomial.  Eliminating `h` from `(33)` gives

```text
9(z+1)^2B=2C^2(2z+1).                                         (34)
```

Because `gcd(z+1,2z+1)=1`, the UFD law forces `z+1 | C`.  If `z+1` were a
unit, polynomiality of `h=(3/2)(z^2-1)/C` would make `z^2=1` and hence
`h=0`, contrary to the present case.  Therefore

```text
z+1=lambda C,                    lambda in k^*,
h=-3lambda+(3/2)lambda^2 C,                                    (35)

B=(2/(9lambda^2))(-1+2lambda C).                               (36)
```

In particular `(31)` is the invertible affine graph

```text
C=(2/(3lambda^2))(A+3lambda).                                  (37)
```

This proves the reverse-graph classification independently of THM-3866.
Applying proved all-polynomial transverse-quotient theorem THM-3866 to
`(37)` then gives the same projective-companion obstruction.

## 7. Frontier after the two graph closures

The vertical theorem closes the last coordinate component omitted by the
forward-graph first-mismatch proof.  Equations `(31)--(37)` also show that
there is no separate high-degree reverse-graph family: every such graph is
vertical or invertibly forward.

The exact remaining branch-design frontier in the depressed-cubic carrier is

```text
- polynomial A1 components that are not graphs over either coordinate;
- more general polynomial-curve normalizations;
- constructions outside this carrier.                          (38)
```

No planar Jacobian counterexample and no all-branch-packet theorem is
claimed.

## 8. Exact replay

Run

```bash
python3 04-computation/jc2_vertical_axis_first_mismatch_thm3870.py
python3 -O 04-computation/jc2_vertical_axis_first_mismatch_thm3870.py
python3 -B 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
python3 -B -O 04-computation/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.py --verify-frozen 05-knowledge/results/jc2_vertical_axis_first_mismatch_independent_audit_thm3870.out
```

Both commands must byte-match
`05-knowledge/results/jc2_vertical_axis_first_mismatch_thm3870.out`.  The
assertion-free companion performs 206 exact gates.  It replays every
truncation `N=0,...,5` and the hostile degree set
`d=0,N+1,N+2,N+3,N+4`, checks the universal perturbation and vertical
classification, and verifies the affine reverse-graph boundary.  The
all-degree proof is `(7)--(27)`; the finite grid is an adversarial control,
not its logical substitute.  Raw-LF SHA-256 hashes are recorded in the
metadata.  The last two commands replay the independent 862-gate hostile
audit described in the frontmatter.
