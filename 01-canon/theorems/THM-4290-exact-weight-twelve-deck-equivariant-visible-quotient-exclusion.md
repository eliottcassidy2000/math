---
id: THM-4290
title: "Exact-weight-twelve deck-equivariant visible-quotient exclusion"
status: >
  PROVED RELATIVE TO THM-4012/4230 + FINITE-EXACT CHARACTER AUDIT PASS. On
  the complete exact-M=12 gate U*Z*D*Lambda!=0, the 12-fold base-change deck
  action makes the specialized positive-genus component map equivariant for
  tau:(S,P)->(xi*S,xi^2*P) on the source and [-omega] on E_0. Its sixth power
  therefore factors the curve map through the genus-two quotient C/<tau^6>.
  THM-4230's saturated visible Hom lattice forces its degree to be
  4(N(alpha)+N(beta)), contradicting the only response degrees 34 and 42.
  Thus the entire exact-M=12 gate is excluded, including all hidden-Hom loci.
  Coefficient walls, seam entry, and JC(2) remain open.
source: root/planar-jc-higher-order-20260830
depends_on:
  - THM-4012-weighted-leading-face-good-elliptic-factor-observer
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
related:
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4260-w0-canonical-node-reciprocal-denominator-attachment-exclusion
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
primary_script: 04-computation/jc23_exact_weight12_deck_equivariance_thm4290.py
primary_output: 05-knowledge/results/jc23_exact_weight12_deck_equivariance_thm4290.out
primary_script_sha256: 697043b3eac9d27de03117c70b4835fba232a117560f53cd914a65d5fa07c7e5
primary_output_sha256: 72d69aab65b3fde8a5dbfe4d6f5387419db940dc550388bf771295ac6edd5d93
hash_basis: raw LF bytes
audit: >
  PASS. A dependency-free exact character audit freezes all source, target,
  face, quotient, elliptic-map, differential, and sixth-power characters
  modulo twelve, then checks the inherited visible degree and response
  residues modulo four. Ordinary and optimized outputs are byte-identical.
  Three independent adversarial derivations checked graph specialization,
  component stability, integral quotient factorization, translations, and
  finite extensions against THM-4012/4230/4249.
---

# THM-4290 -- Exact-weight-twelve deck-equivariant visible-quotient exclusion

**PROVED RELATIVE TO THM-4012/4230 + FINITE-EXACT CHARACTER AUDIT PASS. THE
COMPLETE EXACT-`M=12` GATE IS CLOSED. COEFFICIENT WALLS, SEAM ENTRY, AND
`JC(2)` REMAIN OPEN.**

## 1. Statement

Retain THM-4230's reduced `(2,3)` seam and put

```text
D=W^2-4UZ,                       Lambda=U+W+Z.           (1)
```

Assume the complete exact-weight-twelve gate

```text
U*Z*D*Lambda!=0.                                         (2)
```

Its unique positive-genus special-fibre component and good target are

```text
C: 1-U P^6-W S^2P^5-Z S^4P^4=0,
E_0: Y^2=X^3+1.                                         (3)
```

> **Theorem.** A hypothetical nonautomorphic planar Keller pair in `(2)`
> would specialize to a map `m:C->E_0` whose degree is divisible by four.
> The inherited origin-fibre/degree interface permits only degrees `34` and
> `42`, both congruent to two modulo four. Hence no such Keller pair exists on
> the complete gate `(2)`.

This removes THM-4230's countable hidden-Hom locus from the live frontier.
The later `W=0` lattice, projector, torsion-envelope, node-incidence, and
formal-contact theorems remain valid independent structure and hostile
audits, but are no longer needed to exclude the interior exact-`M=12` gate.

## 2. Exact deck equivariance

Write the inverse target parameter and the weight-twelve base change as

```text
Q=q^(-1)=sigma^12.                                      (4)
```

THM-4230 and THM-4012 use the exact scaled coordinates

```text
s=sigma^(-1)S,              p=sigma^(-2)P,
A=sigma^(-4)X,              C_target=sigma^(-6)Y.       (5)
```

Let `xi=xi_12`, choose `omega=xi^4`, and let the deck generator send
`sigma` to `xi*sigma`. Holding the original coordinates in `(5)` fixed gives

```text
delta_source:(S,P)->(xi*S,xi^2*P)=tau(S,P),
delta_target:(X,Y)->(xi^4*X,xi^6*Y)=(omega*X,-Y)
                                      =[-omega](X,Y).   (6)
```

This is an exact action on the families, not an accidental symmetry of the
special equations. For example,

```text
H_sigma=sigma^12 H(sigma^(-2)P,sigma^(-3)SP)           (7)
```

is invariant, while both terms in

```text
(S^2-P)(1-H_sigma)-sigma^12 S^2/2=0                   (8)
```

have character two. The target equation has the shape

```text
Y^2=X^3+1-c*sigma^8X-d*sigma^12,                       (9)
```

and every term has character zero. Because the original Keller morphism is
defined over `C((Q))`, its scaled base change satisfies

```text
Phi o delta_source=delta_target o Phi.                 (10)
```

On `sigma=0`, equations `(6)` and `(10)` give

```text
m o tau=[-omega] o m.                                  (11)
```

## 3. Why specialization preserves the relation

The needed equivariant-specialization lemma is elementary at the generic
point of the component.

Let a finite tame group act on a DVR model `X`, a proper separated target
`Y`, and an equivariant generic-fibre morphism. The schematic closure of its
graph in `X x Y` is group-stable. At the generic point of any invariant
special-fibre component `D`, assume `X` is normal there, as in the application
below (otherwise first normalize `X`). The local ring is then a DVR.
Properness of `Y` extends the generic morphism uniquely over that DVR, so the
rational map is defined at this generic point. Normalizing `D` therefore gives
a unique component map to `Y_0`, and separatedness specializes the generic
equivariance identity to that map.

THM-4230 proves that the component `C` in `(3)` is smooth, multiplicity one,
and the unique positive-genus component. It is visibly `tau`-stable. Thus the
lemma applies and proves `(11)`. Resolving the graph or the source by point
blowups only adds rational components and preserves the relation on the
strict transform. A later finite base extension merely base-changes the
already-defined equality.

There is also no translation loophole. The scaled map in `(10)` fixes the
target origin and gives `(11)` literally. Even if a different target
identification wrote

```text
m o tau=[-omega]m+T,                                   (12)
```

the induced Jacobian homomorphism loses `T`, and the sixth iterate of `(12)`
does too, because

```text
1+[-omega]+...+[-omega]^5=0.                           (13)
```

## 4. Sixth-power factorization and the mod-four contradiction

The target automorphism `[-omega]` has order six. Therefore `(11)` gives

```text
m o tau^6=m.                                           (14)
```

Now `tau^6:(S,P)->(-S,P)`. By the universal property of the geometric
quotient, `(14)` factors the curve map through the degree-two quotient

```text
pi:C->B=C/<tau^6>.                                     (15)
```

THM-4230 gives exact quotient coordinates

```text
u=1/P,                  x=S^2/P,                 v=W+2Zx,
B: v^2=D+4Z u^6,                                       (16)
```

and two elliptic maps

```text
f_a:B->E_a,       (X,Y)=(u^2,v),
f_b:B->E_b,       (X,Y)=(u^(-2),v/u^3).                (17)
```

After scaling the two `j=0` targets to `E_0`, THM-4230 proves the **saturated
integral** lattice and its degree form

```text
Hom(J(B),E_0)=O f_a direct-sum O f_b,
deg_B(alpha f_a+beta f_b)=2(N(alpha)+N(beta)),          (18)
```

where `O=Z[omega]`. Translation does not affect degree. Combining `(15)` and
`(18)` yields

```text
deg_C(m)=4(N(alpha)+N(beta)) in 4Z.                    (19)
```

All other components in THM-4230's regular special fibre are rational and
therefore map constantly to `E_0`; `C` has multiplicity one. Its existing
proper-flat degree-conservation interface identifies the degree in `(19)`
with the generic response. But THM-4230's exhaustive carrier calculation
gives only

```text
deg(Phi)=34 or 42,                 both =2 mod 4.       (20)
```

Equations `(19)--(20)` contradict each other and prove the theorem.

## 5. The sharper eigenline statement

The sixth-power argument already closes the gate. The full generator gives a
useful refinement. On `(16)`, `tau` sends

```text
u -> xi^10 u,                    v -> v.               (21)
```

Consequently the two maps in `(17)` transform by

```text
f_a o tau=[omega^2]f_a,
f_b o tau=[-omega]f_b.                                  (22)
```

If the based homomorphism induced by `m` is
`alpha f_a+beta f_b`, relation `(11)` and `(22)` give

```text
([omega^2]-[-omega])*alpha f_a=0.                      (23)
```

The coefficient is `omega^2+omega=-1`, a unit. Hence

```text
m-translation=beta f_b,
deg_C(m)=4N(beta).                                     (24)
```

Equivalently, the target differential character is the unique character ten
in THM-4230's multiset `{5,7,8,9,10,11,11}`. The hidden primitive Prym has no
such channel. This is the integral explanation of why a Keller specialization
cannot use the hidden-Hom locus, even where that locus is nonempty.

At `W=0`, THM-4249 supplies an independent lattice check: its integral
projector evaluates to the identity at `T=-omega` and forces the same
specialization onto its `u` eigenline, again of degree `4N`. The basis label
there is normalized differently from `f_b`; the degree conclusion agrees.

## 6. Differential character and the A23 firewall

The deck argument does not prove the ambient-Kahler bridge sought in
THM-4288/4289. Let

```text
eta_0=dX/(2Y).                                         (25)
```

Then

```text
[-omega]^*eta_0=-omega*eta_0=xi^10*eta_0.              (26)
```

The differential in the original target coordinate is instead

```text
dA/(2C_target)=sigma^2*dX/(2Y).                        (27)
```

The factor `sigma^2` has character two and makes `(27)` invariant, but it
also makes its specialization zero. Thus the nonzero limiting form `eta_0`
is a semi-invariant eigensection on the cover, not an ordinary nonzero
ambient differential descended to the raw `Q`-model. This is exactly
compatible with THM-4289's conductor/Jacobian-ideal obstruction. The cyclic
argument bypasses that local descent question rather than solving it.

## 7. Scope and proof-graph change

The conclusion closes precisely the complete coefficient gate `(2)` inside
the inherited exact-weight-twelve reduced seam. It supersedes the **open
status**, not the valid intermediate statements, in THM-4230/4241/4249 and
their `W=0` descendants. In particular:

- hidden Hom factors can exist, but an equivariant Keller specialization
  cannot use them;
- arbitrary degree-`34/42` Hom vectors can exist, but they are not in the
  required deck eigenspace;
- attachment and fat-contact hostiles remain correct for arbitrary resolved
  maps; and
- no raw `A_23` Cartesian square or ambient-Kahler form has been constructed.

The coefficient walls where one of `U,Z,D,Lambda` vanishes, exact-`M=12`
seam entry across all walls, and `JC(2)` remain open. In particular, this
theorem does not silently infer a statement at `Lambda=0` from the interior
gate; THM-4272's collision model needs its own degree-interface audit.

There is nevertheless a proved intrinsic wall survivor. The quotient curve,
deck action, and saturated formula `(18)--(19)` use only `U*Z*D!=0`, not
`Lambda!=0`. Hence, including at `Lambda=0`, every **actual deck-equivariant**
morphism `C->E_0` has degree divisible by four. What is open on the wall is
the Keller identification: current canon proves neither that the resolved
Keller restriction to `C` has degree `34/42` nor that `C` is the sole
degree-carrying component of the required regular specialization. This is
why the intrinsic corollary does not close THM-4272's `A_23` wall.

## 8. Exact audit

The optimization-safe script

```text
04-computation/jc23_exact_weight12_deck_equivariance_thm4290.py
```

independently freezes the exponent vector

```text
(sigma,S,P,X,Y)=(1,1,2,4,6) mod 12,                   (28)
```

checks every face and target character, the sixth powers, quotient
coordinates, the two visible-map characters, the invariant `sigma^2` twist
in `(27)`, and the final mod-four contradiction. Its output is

```text
05-knowledge/results/jc23_exact_weight12_deck_equivariance_thm4290.out
```

The script audits exact character arithmetic; the graph-specialization lemma
and the saturated lattice/degree theorem are proved textually above and in
THM-4230.

**QED.**
