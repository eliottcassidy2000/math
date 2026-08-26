---
id: THM-4143
title: "Two-term collision-wall critical boundary and monodromy exclusion"
status: >
  PROVED RELATIVE TO THM-3827/4053/4103/4120/4122/4130/4138/4141 +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. On the live exact-M=8 two-term wall
  delta*theta!=0 and delta+theta=0, the Phi=0 stratum is impossible. The
  K!=0 and K=0 strata have affine critical lengths 19 and 17, normalized
  infinity packets (6,6,3,2,2,1) and (6,6,5,1), and no compatible
  transitive monodromy response. Hence the entire wall is empty. Together
  with THM-4130/4138/4141, this empties THM-4053's exact-M=8 trichotomy on
  the live b=d=0 reduced (2,3) seam. M>=9, entry, other cells, JC(2), and
  DC(2) remain OPEN.
source: codex-planar-jacobian-cycle6-63-20260825
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
  - THM-4141-delta-d-collision-wall-boundary-monodromy-exclusion
related:
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
script: 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143.py
output: 05-knowledge/results/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143.out
independent_audit_script: 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143_independent_audit.out
script_sha256: c5a711a8852ad0d1a9e15c794c7f9b513ee1fb1238d7d2dac2686ab9ea39d03b
output_sha256: ffd874eab873bbbaf989b83fd90eef9b6a4f1d7ded3f1512b5463414029329fd
independent_audit_script_sha256: 493529db444b7b0665c55549929ebaf80101cc20f194428086e7e3476ebb76a2
independent_audit_output_sha256: 9b35bed2a4cdf1458bbbc10a1d26ae6d0d6cea0487b57bd707553bd65b0a85db
semantic_sha256: b91c8b2d1205fb15b7d5c7a77e81d30a4bb0fc0442c6756ea8afdd98f1b15d77
independent_semantic_sha256: 1d2e2cd6d6fde167cdf32e8c27b273b1e5aae82a9067471ff6ad6369f9049d63
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone SymPy calculation reconstructs the complete normalized
  wall polynomial, excludes Phi=0 by four exact critical-value eliminations,
  factors both live critical resultants, rebuilds both Newton polygons and
  the ordinary-node normalization, records the exact response ledger, and
  checks every monodromy identity case plus sharp hostiles.
  Normal, optimized, and hash-seeded replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A clean-room SymPy referee imports no primary certificate. It uses
  a separately derived source equation and valuation charts, independently
  obtains the two critical lengths, packets, residue-field labels, response
  degrees, and monodromy exclusions in 113 live checks. Its seven-node
  collision model is corroborating local evidence only and is not a proof
  dependency. Normal, optimized, and two hash-seeded replays byte-match.
---

# THM-4143 -- the two-term collision wall is empty

**PROVED RELATIVE TO THM-3827/4053/4103/4120/4122/4130/4138/4141 +
VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2) REMAINS OPEN.** Work over
`C` and retain THM-4053's live `b=d=0` reduced `(2,3)` cell at exact total
residual weight eight.

## 1. Theorem and inheritance

> **Theorem.** No nonautomorphic planar Keller pair lies on the two-term
> collision wall
>
> ```text
> delta*theta!=0,                 delta+theta=0.        (1)
> ```

The inheritance pass is deliberately wall-sensitive:

- the closest proved mechanism is THM-4141's critical-length versus
  boundary-response monodromy obstruction;
- the canonical hostile is THM-4134, where a repeated boundary root changes
  both genus and response degree;
- the corrected near miss is an arbitrary-subset response argument, which
  loses the residue-field labels and is not valid;
- the least-used sidecar is the exact contraction
  `P^4-PY^2=TP^3`, which lowers the critical projection by one degree.

The repair is to determine every response over `C(q)` before doing any
permutation arithmetic. This leaves exactly the degrees `20`, `16`, and
`18`, each of which fails for a different sharp reason.

## 2. The contracted normalized polynomial

Use the normalization of THM-4130 and put

```text
P=T+X^2T^2,                    Y=XTP,
Delta=a^17 delta,              Theta=a^17 theta,
Phi=a^(29/2) phi,
K=a^12 kappa=2848/45-(7/6)Delta.                       (2)
```

Wall `(1)` is `Theta=-Delta` with `Delta!=0`. Its structural identity is

```text
P^4-PY^2=P^3(P-X^2T^2)=TP^3.                           (3)
```

Consequently the **entire** normalized source polynomial is

```text
G=-X^2T/2-3P+(8/3)P^2-(1376/135)P^3
  +K Y^2+Phi P^2Y+Delta TP^3.                          (4)
```

There is no omitted lower or top term. Split the proof into `Phi=0`, then
`Phi!=0,K!=0`, and finally `Phi!=0,K=0`. By `(2)`, the last equality is

```text
K=0                         iff Delta=5696/105.         (5)
```

## 3. The `Phi=0` subwall is impossible

On the line `X=0`, equation `(4)` becomes

```text
g(T)=-3T+(8/3)T^2-(1376/135)T^3+Delta T^4.             (6)
```

Every root of `g'` is a critical point of `G`. In a Keller realization
`G=E(A,C)`, the Hessian congruence at the two target nodes makes all source
critical points Morse, and their values belong to `{0,1/2}`. Thus the monic
critical-value eliminant

```text
Res_T(g'(T),Z-g(T))                                     (7)
```

would have to equal one of

```text
Z^j(Z-1/2)^(3-j),                   0<=j<=3.            (8)
```

For each of the four choices, exact coefficient comparison gives three
polynomials in `Delta` whose gcd is `1`. The primary and independent audits
derive this separately, including all repeated-root specializations. Hence
no value of `Delta` satisfies `(8)`, and

```text
Phi=0 is empty.                                         (9)
```

## 4. Exact affine critical lengths for `Phi!=0`

Since `G_X` is divisible by `T`, set

```text
f=G_X/T,                         h=G_T.                 (10)
```

Both have common leading row `7Phi T^6` as polynomials in `X`, of degrees
six and seven. For `K!=0`, exact elimination gives

```text
Res_X(f,h)=7T^30 Phi(6T+1)^2 Q_15(T),
deg Q_15=15,
[T^15]Q_15=576 Delta^4K^4,
Q_15(0)=-(7203/32)Phi^4.                               (11)
```

For `K=0`, hence `(5)`, it gives

```text
Res_X(f,h)=7T^30 Phi(6T+1)^2 Q_13(T),
deg Q_13=13,
[T^13]Q_13=
 2185746832794987941330944/16544390625,
Q_13(0)=-(7203/32)Phi^4.                               (12)
```

All displayed endpoints are nonzero on their stated strata. The factor
`T^30` is a Sylvester degree-drop artifact, because

```text
f(X,0)=-X,                    h(X,0)=-(X^2+6)/2.        (13)
```

For `T!=0`, the common nonzero leading row excludes a critical point at
`X=infinity`. The factor `(6T+1)^2` is the universal Morse pair

```text
T=-1/6,                       X^2=6,       G=1/2,       (14)
```

while the actual critical ideal `(Tf,h)` restores the pair omitted by
`(13)`:

```text
T=0,                          X^2=-6,      G=0.         (15)
```

Their Hessian determinants are respectively `-6` and `6`. In a hypothetical
Keller realization,

```text
grad(E o F)=DF^t grad(E),
Hess(E o F)=DF^t Hess(E)DF                              (16)
```

at a critical point. Since `det DF=1`, every point counted by `(11),(12)`
is Morse and the schemes are reduced. If `r_i` is the number of affine
preimages of the target node of value `i/2`, then

```text
K!=0:       r_0+r_1=19,
K=0:        r_0+r_1=17,
both:       r_0,r_1>=2.                                (17)
```

## 5. Geometric connectedness

The source-fibre transitivity below cannot be inferred merely from a
primitive equation over `C(q)`. Use instead the closed-polynomial factor
theorem cited and audited in THM-3827. If the geometric generic fibre of the
polynomial `J=E(A,C)` were disconnected, then

```text
J=R(J_0),                         deg R>1.              (18)
```

A root `c` of `R'` would make the nonempty curve `J_0=c` lie in `Crit(J)`.
That contradicts the finite critical schemes in `(17)`. Thus the geometric
generic source fibre is connected on both live `Phi!=0` strata.

## 6. Complete boundary and normalization

Put

```text
s=XT,              p=T+s^2,              t=p-s^2,
gamma=-1/2,        lambda=-3,             alpha=8/3,
epsilon=-1376/135.                                      (19)
```

Then

```text
G=gamma s^2/t+H(s,p),
H=lambda p+alpha p^2+epsilon p^3
  +K s^2p^2+Phi sp^3+Delta p^3(p-s^2).                 (20)
```

For `Q=q^-1`, a generic source fibre has equation

```text
F_Q=(s^2-p)(1-QH)+gamma Qs^2=0.                        (21)
```

After coincident monomials are combined, its Newton polygons are

```text
K!=0: conv{(0,1),(2,0),(4,2),(4,3),(0,5)},
       2Area=26, boundary=10, arithmetic genus=9;

K=0:  conv{(0,1),(2,0),(4,3),(0,5)},
       2Area=24, boundary=8,  arithmetic genus=9.       (22)
```

In the generic polygon, call these vertices `A,B,C,D,E` in the displayed
order. The simple boundary inventory is

| stratum | edge/place | residue label | ramification index |
|---|---|---|---:|
| `K!=0` | `AB` | rational | `1` |
| `K!=0` | `BC` | one quadratic closed point | `2,2` |
| `K!=0` | `CD` | rational | `3` |
| `K=0` | replacement `BD` | rational | `5` |

The four roots on the vertical edge `EA` have `s=0,p!=0`; they are affine
points in the original source chart, not punctures at source infinity. The
only nonsquarefree infinity face is the repeated `DE` edge, handled next.

Set

```text
s=1/z,                  p=(1-a)/z^2,                  r=1-a,
L(a,z)=z^10 F_Q(1/z,r/z^2).                            (23)
```

Exact expansion gives

```text
L=a[z^8+QDelta r^3a-QPhi r^3z
      -Q(epsilon r^3+Kr^2)z^2-Qalpha r^2z^4
      -Qlambda rz^6]+gamma Qz^8.                      (24)
```

At `(a,z)=(0,0)` its tangent cone is

```text
Qa(Delta a-Phi z).                                     (25)
```

Because `Delta Phi!=0`, this is an ordinary node. Its two normalization
branches begin

```text
a=(Phi/Delta)z+O(z^2),
a=(gamma/Phi)z^7+O(z^8).                               (26)
```

The Keller residue identity of THM-4103 becomes, up to a nonzero sign,

```text
omega=Q ds/(F_Q)_p=Qz^6 dz/L_a.                        (27)
```

On the two branches in `(26)`, `L_a=+/-QPhi z+O(z^2)`.
Thus `omega` has order five on each: the node normalizes into two rational
punctures of index six. Since the arithmetic genus in `(22)` is nine, this
node gives the upper bound `g<=8` for the smooth projective normalization.

At the other corner put `u=1-a`. When `K!=0`, there is a rational branch

```text
u=(K/Delta)z^2+O(z^3)                                  (28)
```

of index three, and two branches `u=wz^3+O(z^4)` with

```text
K w^2=q+gamma=q-1/2,                                   (29)
```

each of index two. Since `q-1/2` has odd valuation at `q=1/2`, `(29)` is
one irreducible quadratic closed point over `C(q)`. When `K=0`, the primitive
edge equation is

```text
Delta u^3+(q-1/2)z^8=0,                 gcd(3,8)=1.    (30)
```

Its toric coordinate `u^3/z^8` satisfies a linear equation over `C(q)`, so
it is one rational place. Writing locally `z=tau^3,u=c tau^8` in a splitting
field and applying `(27)` gives differential order four, hence index five;
no rationality of the auxiliary coefficient `c` is being assumed.

The complete packets are therefore

```text
K!=0:             (6,6,3,2,2,1),
                   labels rational (6,6,3,1)
                          + quadratic (2,2);

K=0:              (6,6,5,1), all rational.             (31)
```

Both displayed packets have ramification defect fourteen. Riemann--Hurwitz
for the nonconstant map to the elliptic target gives

```text
14<=deg R=2g-2.                                        (31a)
```

Thus `g>=8`, while the ordinary node above gave `g<=8`. Hence `g=8`, equality
holds in `(31a)`, and `(31)` is the whole ramification divisor. In
particular, there is no additional singularity lowering the normalization
genus and no hidden affine branch point.

## 7. Exact response degrees

THM-4120 proves for the target generic fibre that

```text
E_q(C(q))={O}.                                         (32)
```

Every rational puncture in `(31)` therefore maps to the target origin. The
quadratic pair in `(29)` responds together. If its finite images coincided,
that image would be Galois invariant and hence a finite `C(q)`-point,
contradicting `(32)`. Thus a finite response consists of two distinct
conjugate points. Since the polynomial target coordinates have no affine
pole, no affine source point maps to target infinity. The exhaustive list is

```text
K!=0, BC->O:       n=20, O-packet (6,6,3,2,2,1);
K!=0, BC finite:  n=16, O-packet (6,6,3,1),
                         two distinct index-two carriers;
K=0:               n=18, O-packet (6,6,5,1).           (33)
```

This Galois-locked list is the step missing from the discarded
arbitrary-subset argument.

## 8. Fixed-sheet transport with the finite carrier

The only additional issue in `(33)` occurs for `n=16`. Choose a constant
square root of `K` and put `rho=sqrt(K)w`; equation `(29)` becomes

```text
q=1/2+rho^2.                                           (34)
```

This is the same quadratic target base change audited in THM-4138. Its
source-wall-independent argument is re-established here: THM-4122 makes the
normalization of a horizontal nonproperness component `A1`; degree
factorization and `(32)` force degree two over the `q`-line; and the exact
rank-one Mordell--Weil calculation on

```text
y^2=x^3-3x+2+rho^2                                     (35)
```

leaves the polynomial sections `+-P,+-2P,+-3P`. Their coordinate degree
pairs are `(0,1),(0,1),(2,3)`. THM-4122 requires the positive intrinsic pole
pair `(2rho_C,3rho_C)`, so only the `+-3P` carrier survives. Hence the two
finite branch values are exactly the marked conjugate points `Q_+,Q_-` of
THM-4138.

For an affine inverse of target node `o_i`, Keller etaleness gives local
coordinates

```text
E-q_i=uv.                                              (36)
```

A nearby Milnor core therefore has one closed degree-one lift for every such
inverse. Resolve the compactified pencil graph and delete the finite extra
Hurwitz discriminant. After shrinking the base, the morphism of proper
smooth curve families is quasifinite and therefore finite. Deleting the
complete marked branch divisor `O` and, for `n=16`, `Q_+,Q_-`, makes it
finite etale. THM-4138's parallel-core construction avoids a carrier/node
collision, and Ehresmann transport only conjugates sheet labels. Vertical
nonproper components are confined to the deleted special fibres by the
`A1`-normalization gate.

Let `X,Y` be the transported vanishing permutations around the two nodal
fibres. Equations `(17),(36)` give

```text
#Fix(X)>=r_0,                    #Fix(Y)>=r_1.           (37)
```

Removing finitely many points from the connected projective source does not
disconnect it. Hence the resulting punctured cover has transitive monodromy.

## 9. Full-boundary degrees `20` and `18`

For a permutation `sigma`, write `supp(sigma)` for its moved letters. From
`(17),(37)`, the two full-boundary cases satisfy

```text
n=20: |supp(X)|+|supp(Y)|<=2n-19=n+1,
n=18: |supp(X)|+|supp(Y)|<=2n-17=n+1.                  (38)
```

The once-punctured target torus is generated by `X,Y`, and its origin
meridian is their commutator up to inversion. If one handle permutation is
the identity, transitivity would force the other to be an `n`-cycle, but
`(17),(37)` give it at least two fixed sheets.

If both are nonidentity, transitivity and `(38)` force their supports to
cover all `n` letters and meet in exactly one pivot. Any second nontrivial
cycle of either permutation would form a separate orbit, so each is one
cycle on its support. The commutator of two cycles meeting at one pivot has
type

```text
n=20: (3,1^17),
n=18: (3,1^15).                                       (39)
```

Neither is the corresponding origin packet in `(33)`. This excludes both
full-boundary degrees.

## 10. Finite-carrier degree `16`

Now the punctured-torus monodromy is generated by `X,Y` and the two distinct
carrier meridians `Z_+,Z_-`. Each `Z` is a transposition. Equations
`(17),(37)` give

```text
|supp(X)|+|supp(Y)|<=2*16-19=13.                       (40)
```

For every nonidentity permutation,

```text
ind(sigma)=16-#Cycles(sigma)<=|supp(sigma)|-1.          (41)
```

If both `X,Y` are nonidentity, the total generator index is at most

```text
(13-2)+ind(Z_+)+ind(Z_-)=13.                           (42)
```

If exactly one is the identity, it is at most `12+2=14`; if both are
identities, it is two. But a transitive permutation system on sixteen
letters has total generator index at least fifteen: starting from sixteen
singleton orbits requires at least fifteen orbit mergers. Every case
contradicts transitivity. Thus `n=16` is impossible.

This bound is sharp. If the critical length in `(40)` were `18`, a
fourteen-cycle and two attaching transpositions could generate all sixteen
sheets. The exact nineteenth critical point is essential.

## 11. Consequence and replay

Sections 3, 9, and 10 exhaust every coefficient and response stratum, so

```text
the two-term exact-M=8 collision wall is empty.         (43)
```

THM-4053's other exact-weight-eight alternatives are already excluded by
THM-4130, THM-4138, and THM-4141. Therefore the complete exact-`M=8`
trichotomy is empty on this live reduced seam. This does **not** treat
`M>=9`, entry into the seam, another reduced cell, `JC(2)`, or `DC(2)`.

Replay with

```text
python3 -B 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143.py
python3 -B -O 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143.py
PYTHONHASHSEED=29 python3 -B 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143.py
python3 -B 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143_independent_audit.py
python3 -B -O 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143_independent_audit.py
PYTHONHASHSEED=123 python3 -B 04-computation/jc23_two_term_collision_wall_critical_boundary_monodromy_thm4143_independent_audit.py
```

All corresponding streams byte-match their frozen outputs. **QED.**
