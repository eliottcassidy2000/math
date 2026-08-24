---
id: THM-3882
title: "Immersed-dual one-place Wronskian projection criterion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  an everywhere-immersed plane normalization of arbitrary genus, a line section of its dual
  has divisor twice the base fibre of the corresponding point projection
  plus that projection's ramification divisor.  One-point support occurs
  exactly when genus zero is forced, the projection has degree one, and its
  base fibre is a single point of multiplicity d-1.  Immersion therefore
  leaves only the smooth-conic boundary.  Without immersion the exact reduced
  formula subtracts the common tangent-base divisor, and a sharp family shows
  one-place lines in every degree.  Consequently no irreducible sextic with
  exactly six A2 cusps and four A1 nodes has an affine chart with A1
  normalization.  This closes an equisingularity architecture, not JC(2).
source: root / post-THM-3879 dual-Wronskian reframe, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_quartic_c3_construct, 2026-08-23).  The
  audit rederived the homogeneous tangent-map basepoint criterion, the exact
  `2D+Ram` divisor identity, the Riemann--Hurwitz/local-index equivalence,
  and the immersed one-address boundary.  It independently checked the
  Pluecker/reflexivity inference for the full `6A2+4A1` packet and the
  THM-3879 `(3,3)` node specialization.  The proof identifies the line pullback
  with the Wronskian of a point projection, factors its base divisor with
  multiplicity two, invokes Riemann--Hurwitz only after recording the exact
  projection degree, and checks the local immersion boundary.  The exact
  companion verifies the Wronskian factorization, every degree ledger through
  d=40, the Pluecker packet, and the THM-3879 node/projection specialization.
  Normal and optimized runs byte-match the frozen transcript.  A second,
  import-independent 25,899-gate audit proves the arbitrary-genus extension,
  the corrected nonimmersed formula `2D+Ram-E_nu`, the all-degree sharp
  hostile, and a separate degree/reflexivity check of dual immersion.  Its
  normal and optimized streams also byte-match the frozen output.
depends_on:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
related:
  - THM-3841-deleted-ramification-three-puncture-jelonek-nonentry
  - THM-3851-tricuspidal-quartic-rank-two-two-place-tradeoff
script: 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
output: 05-knowledge/results/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.out
script_sha256: 94374be74db4c08049d18a06316f53a342b20340bdd3dcdee12d1810bb060854
output_sha256: 8c9fa0ede45e359ea08a928506ba227662ce32d8825c9b9f1cc3d392f0d6438b
semantic_sha256: 6452278a5afd71f244134a6871e247195c6af4f80f641c1b05c7f6b6dd6b8fc2
second_independent_script: 04-computation/jc2_rational_dual_one_place_wronskian_second_independent_audit_thm3882.py
second_independent_output: 05-knowledge/results/jc2_rational_dual_one_place_wronskian_second_independent_audit_thm3882.out
second_independent_script_sha256: d4de70e66c98d7e342adddc77d0c4461e04ee638de2fe272732b34c4f79a4943
second_independent_output_sha256: 47dedac2d6b134eaf03419261c36aadd60b90384496ee54557a3082cd7a58ce9
second_independent_semantic_sha256: f577f0ff2345f79c0ac654713fa0325376def81ea79c7857a75f37777bae7b5f
hash_basis: raw LF bytes
---

# THM-3882 -- an immersed one-place dual line is a degree-one projection

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
Work over an algebraically closed field `k` of characteristic zero.  Let

```text
nu:X -> C subset P2                                           (1)
```

be the normalization of an irreducible nonlinear plane curve of degree
`d>=2`, where the smooth projective curve `X` has genus `q`.  Assume that
`nu` is an immersion at every normalization point, and write
`nu_dual:X -> C_dual` for its tangent-line map.

For a point `P in P2`, write `H_P` for the line in the dual plane whose
points are the lines through `P`.  Choose two independent linear forms
`L_0,L_1` vanishing at `P`.  Their pullbacks are sections of
`L=nu^*O_P2(1)`.  Factor their complete common base divisor:

```text
L_0 o nu = h x,              L_1 o nu = h y,
D_P=div(h),                  m=deg D_P,
e=deg L(-D_P)=d-m.                                           (2)
```

Then `[x:y]` is the resolved projection from `P`, a morphism

```text
phi_P:X -> P1                        of degree e,              (3)
```

and the exact divisor identity is

```text
(nu_dual)^* H_P = 2D_P + Ram(phi_P).                          (4)
```

In particular,

```text
supp((nu_dual)^*H_P)={p}
iff
q=0,  D_P=(d-1)p,  and deg(phi_P)=1.                          (5)
```

Because `nu` is immersive, a base fibre supported at one point has
coefficient one there.  Thus `(5)` forces `d=2`: one-place support itself
forces rationality, and the smooth conic is the unique positive boundary.
For every immersed plane normalization of degree `d>=3`, every line on the
dual curve meets its normalization in at least two support points.

## 1. A dual line is a projection Wronskian

Move `P` to `[0:0:1]` and write the parametrization as `[X:Y:Z]`.  The line
`H_P` asks whether the tangent to `C` passes through `P`.  In a local
parameter `t`, its pullback is therefore the Wronskian

```text
W_P=X Y'-Y X'.                                                 (6)
```

In a local trivialization write `X=hx` and `Y=hy`.  The derivative of `h`
cancels exactly:

```text
W(hx,hy)=h^2 W(x,y).                                          (7)
```

The zero divisor of the coprime Wronskian `W(x,y)` is precisely the
ramification divisor of `[x:y]`.  Formula `(7)` proves `(4)` locally, hence
globally.  It also explains the otherwise easy-to-miss factor two on the
fibre over `P`.

Immersion is the exact hypothesis which prevents a further cancellation.
The three tangent coordinates are sections of `L^2 tensor K_X`; at a
nonimmersed normalization address they acquire a common factor, and the
reduced dual map divides it out.  No such common tangent factor exists here.

The degree ledger is consequently

```text
deg(2D_P)+deg Ram(phi_P)
 =2m+(2e+2q-2)
 =2d+2q-2
 =deg (nu_dual)^*H_P.                                        (8)
```

## 2. Riemann--Hurwitz makes the criterion exact

Suppose the left side of `(5)` is supported at one point `p`.  Both effective
summands in `(4)` must then be supported at `p`.  Riemann--Hurwitz gives

```text
deg Ram(phi_P)=2e+2q-2,                                       (9)
```

whereas the ramification coefficient at one point is at most `e-1`.  Hence
`2e+2q-2<=e-1`, or `e+2q-1<=0`.  With `e>=1` and `q>=0`, the unique possibility
is `(q,e)=(0,1)`.  The ramification divisor vanishes, `(2)` gives `m=d-1`,
and the support assumption says `D_P=(d-1)p`.

Conversely, if `D_P=(d-1)p`, then `e=1`, so `(4)` is simply

```text
(nu_dual)^*H_P=2(d-1)p.                                      (10)
```

This proves `(5)`, including the fact that rationality is a conclusion rather
than a hypothesis.  Finally, when the fibre of `nu` over `P` has only the
address `p`, immersion says that at least one of the two local coordinates
through `P` has a nonzero first derivative.  Thus the common vanishing order
in `(2)` is one.  Equality `d-1=1` forces `d=2`.  The conic/tangent-line case
is the sharp positive boundary.

## 3. The entire six-cusp/four-node sextic packet is two-place

Let `Gamma` be an irreducible plane sextic having exactly six ordinary `A2`
cusps and four ordinary `A1` nodes and no other singularities.  The classical
Pluecker formulas give

```text
genus(Gamma) = (5*4)/2 - 6 - 4 = 0,
deg(Gamma_dual) = 6*5 - 3*6 - 2*4 = 4,
inflection weight = 3*6*4 - 8*6 - 6*4 = 0.                  (11)
```

Thus `Gamma_dual` is a rational quartic, and the zero inflection weight says
that its normalization map is everywhere immersive.  Biduality identifies
`Gamma` with the dual of that immersed quartic.  Applying `(5)` with `d=4`
shows that every line section of `Gamma` has at least two normalization
support points.

Therefore no member of the whole equisingularity class

```text
degree 6;  singular packet 6A2+4A1                            (12)
```

has an affine chart whose normalization is `A1`.  This is independent of a
chosen torus equation, contact conic, or coordinates.  A one-place
counterexample design cannot preserve the exact THM-3879 singular packet.

## 4. Why the THM-3879 best line has type (3,3)

For the trinodal quartic in THM-3879, the dual line `C=0` corresponds to the
primal node `P=[0:0:1]`.  Its two projection coordinates have the factorization

```text
X=ST(S^2-2T^2),                  Y=ST(S^2-T^2).               (13)
```

Hence

```text
D_P=[S=0]+[T=0],              e=2.                            (14)
```

The residual degree-two projection is ramified once at each of those two
addresses.  Formula `(4)` gives

```text
(nu_dual)^*(C=0)=3[S=0]+3[T=0],                               (15)
```

exactly the pullback `2S^3T^3` found in THM-3879.  The line's two-place
optimality is therefore not an accidental sparse coefficient pattern: its
two node branches each carry one base unit twice and one ramification unit.

## 5. Exact nonimmersion correction and sharp boundary

If `nu` is only generically immersive, let `E_nu` be the common base divisor
of its three unreduced tangent coordinates.  Cancelling that divisor to form
the reduced dual map changes `(4)` to the exact identity

```text
(nu_dual,reduced)^*H_P=2D_P+Ram(phi_P)-E_nu.               (16)
```

The subtraction is effective coordinatewise because `E_nu` divides every
tangent coordinate.  Thus the immersion hypothesis is structural, not a
technical convenience.  It has a sharp hostile in every degree.  For `d>=2`,
put

```text
nu_d([S:T])=[S T^(d-1):T^d:S^d].                           (17)
```

This map is birational onto a degree-`d` curve because `X/Y=S/T` on the dense
torus.  Its raw tangent coordinates are

```text
(-d^2 S^(d-1)T^(d-1), d(d-1)S^dT^(d-2), dT^(2d-2)).       (18)
```

They share `T^(d-2)`.  After cancellation, the third coordinate is `dT^d`,
so its corresponding dual line is supported at one address.  The divisor
ledger is

```text
D_P=(d-1)[T=0],  Ram(phi_P)=0,  E_nu=(d-2)[T=0],
2D_P+Ram-E_nu=d[T=0].                                      (19)
```

The local branch `[u^(d-1):u^d:1]` is immersive only for `d=2`; for every
`d>=3` this is exactly the forbidden nonimmersion.  Formula `(16)`, rather
than `(4)`, is therefore the correct boundary statement.

## 6. Counterexample-search consequence and scope

This theorem constructs no Keller map and proves no case of `JC(2)`.  It
does remove a full geometric architecture from the search: varying the
contact conic or torus presentation inside the `6A2+4A1` sextic packet can
never repair the one-place defect.

There are now only three honest ways around this obstruction:

1. make the primal normalization nonimmersive, so a common tangent factor is
   cancelled as in `(16)`;
2. change the dual/singularity packet, necessarily paying a different
   Pluecker and Cardano divisor ledger; or
3. abandon the dual construction and build the one-place branch directly.

Any such move must still preserve the connected codimension-one-unramified
`C3` packet of THM-3879 while avoiding its globally monogenic cubic order.

Run

```text
python3 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
python3 -O 04-computation/jc2_rational_dual_one_place_wronskian_projection_criterion_thm3882.py
python3 04-computation/jc2_rational_dual_one_place_wronskian_second_independent_audit_thm3882.py
python3 -O 04-computation/jc2_rational_dual_one_place_wronskian_second_independent_audit_thm3882.py
```

and compare each normal/optimized pair byte-for-byte with its frozen output.
