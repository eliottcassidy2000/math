---
id: THM-2596
title: "Modular free factors, Farey children, and the Gram-owner cocycle"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  The Bass--Serre tree of PSL_2(Z)=C_2*C_3 is
  (2,3)-biregular; contracting its degree-two vertices gives the trivalent
  Farey-dual tree, whose rooted noninitial branching is binary.  The two
  positive Farey children are parabolic words in the torsion generators.
  On a THM-2056 flank, the exact faithful state is the Gram-owner pair
  (U^T U,U^T w), transforming by congruence under every GL_2(Z) basis
  move.  Equal endpoint defects do not determine the mediant: two acute
  unimodular flanks with endpoint values (-90,-89) have child defects -177
  and +1.  Active order-three motion can likewise change a nonnegative
  defect into a negative one, so the THM-2056 Euclidean defect certificate
  is not invariant under active modular motion.  The three Berggren
  branches form a disjoint PGL_2(Z) reduction cross-section of (0,1), not
  a C_3 action; their triple matrices are identity mod 2 and not order
  three.  This repairs the binary/ternary "one object" analogy without
  proving LRC(14) or identifying a V_4 torsor.
source: codex-2026-07-27-modular-transfer
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2467-bicycle-spaces-of-the-star-flip-split
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
script: 04-computation/modular_farey_gram_owner_cocycle_thm2596.py
output: 05-knowledge/results/modular_farey_gram_owner_cocycle_thm2596.out
script_sha256: 096ec9b646794d981ec02951fa6144efb9496b1d7c1f80f6dfb5028944bbdf12
output_sha256: 450d30bacdfdcc6b56353de1e608e5cec37f432abe4b3ccf36671fb8e66c831b
hash_basis: normalized repository blobs (LF)
---

# THM-2596 -- the modular tree needs a metric sidecar

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**  The group, matrix, and hostile statements are proved
below.  The dependency-free companion checks the displayed identities and
`582,340` exact integer/rational instances.  The LRC paragraph is a
faithful reformulation of the proved THM-2056 defect, not a new LRC gate.

The tempting slogan is that a binary fraction tree and a ternary triple
tree are the two free factors of

```text
PSL_2(Z) = C_2 * C_3.                                      (1)
```

There is a rigorous core, but the literal branch-count reading is false.
The faithful common object is a congruence state carried by integral basis
moves.  It retains a quadratic Gram coordinate which every static endpoint,
`V_4`, partial-cube, or parity shadow discards.

## 1. What the free product tree actually says

Use the integral representatives

```text
S = [ 0 -1 ]       C = [ 0 -1 ]
    [ 1  0 ],          [ 1  1 ].                           (2)
```

Then `S^2=-I` and `C^3=-I`, so their projective classes have orders two
and three.  The Bass--Serre tree has vertices

```text
PSL_2(Z)/<S>  disjoint_union  PSL_2(Z)/<C>,                (3)
```

and one edge for every group element, joining its two cosets.  The two
vertex types therefore have degrees two and three.  Freeness of the normal
form proves that this incidence graph is a tree.  Contracting every
degree-two vertex gives the trivalent dual tree of the Farey tessellation.
After one root edge is chosen, each nonroot trivalent vertex has one parent
and two children.  That is the exact origin of the binary rooted grammar.
It is not a second free factor of valence two.

The positive Farey/Stern--Brocot children are the parabolic matrices

```text
L = [1 1] = -S C,          R = [1 0] = S C^(-1)            (4)
    [0 1]                      [1 1]
```

in projective notation.  If `U=[u v]` is an ordered unimodular flank, then

```text
U L = [u,u+v],             U R = [u+v,v].                  (5)
```

Thus the two children are **words in** the torsion generators.  They are
not the order-two and order-three generators themselves.

## 2. The faithful THM-2056 state

Fix any integral owner vector `w` and any `U in GL_2(Z)`.  For a coordinate
column `z`, put

```text
d=U z,
G=U^T U,                   ell=U^T w,
F_w(d)=||d||^2-91 w dot d.                                (6)
```

Then exactly

```text
F_w(Uz)=z^T G z-91 ell^T z.                               (7)
```

For every `g in GL_2(Z)`, change the basis passively by `U'=Ug` and the
coordinates by `z'=g^(-1)z`.  The state update is

```text
G'   =g^T G g,
ell' =g^T ell.                                             (8)
```

Equations (7)--(8) give

```text
U'z'=Uz,                 F_w(U'z')=F_w(Uz).                (9)
```

This proves the Gram-owner covariance.  It is a genuine action law because
two successive changes multiply their matrices.  In coordinates

```text
G=[a b],                  ell=(r,s),
  [b c]
```

the two Farey child updates are

```text
L: G -> [a,       a+b;     a+b, a+2b+c], ell -> (r,r+s),
R: G -> [a+2b+c, b+c;      b+c, c      ], ell -> (r+s,s). (10)
```

On a fixed signed hull-owner cone of THM-2056, `F_w>=0` is precisely its
determinant certificate.  Equation (7) is therefore the exact
finite-dimensional state
on which a modular/Farey dynamic program may run.  The scalar endpoint
defects are only two evaluations of this quadratic polynomial.

The mediant law isolates the missing coordinate:

```text
F_w(u+v)=F_w(u)+F_w(v)+2 u dot v.                          (11)
```

It is the off-diagonal Gram entry `u dot v`, not a tree label, which decides
whether the child repairs the endpoint defects.

## 3. Sharp endpoint/V4/partial-cube hostile

Take `w=(1,0)` and compare the two acute unimodular flanks

```text
(u,v_0)=((1,0),(1,1)),
(u,v_1)=((1,0),(90,1)).                                   (12)
```

They have identical endpoint data, even before defect-gate pass/fail
scalarization:

```text
F_w(u)=-90,             F_w(v_0)=F_w(v_1)=-89.             (13)
```

Both determinants are one and both dot products are positive.  Yet

```text
F_w(u+v_0)=F_w(2,1)=-177,
F_w(u+v_1)=F_w(91,1)=1.                                  (14)
```

Consequently none of the following determines the Farey child:

```text
the two endpoint defect-gate bits (a four-state V4 toggle set);
the two exact endpoint defect values;
the abstract edge/path joining the endpoints;
any partial-cube or graceful label of that abstract carrier.              (15)
```

The same carrier and endpoint labels occur in (12), while the target
predicate differs at (14).  The Gram cross term is a load-bearing sidecar.

There is also no **active** modular symmetry of the Euclidean defect.  Let

```text
w=(1,0), d=(91,-45),
d'=C d=(45,46),             w'=C^(-T)w=(1,1).              (16)
```

The linear pairing is preserved:

```text
w dot d=w' dot d'=91.                                     (17)
```

But the Euclidean norm is not:

```text
F_w(d)=2025,                 F_(w')(d')=-4140.             (18)
```

Thus the order-three generator can change the sign of the displayed
Euclidean defect.  Passive basis covariance (8) is lawful; active invariance
of the THM-2056 defect certificate without the transformed Gram metric is
false.  No actual LRC row or full LRC-safe predicate is transported by this
hostile.

## 4. The ternary Pythagorean tree is a reduction cross-section

Parameterize a primitive Pythagorean triple by

```text
x=m/n in (0,1),
(a,b,c)=(n^2-m^2, 2mn, n^2+m^2).                          (19)
```

Its three standard positive children act on `x` by

```text
A(x)=1/(2-x),       matrix [ 0 1],     image (1/2,1),
                            [-1 2]

B(x)=1/(2+x),       matrix [0 1],      image (1/3,1/2),
                            [1 2]

C_0(x)=x/(2x+1),    matrix [1 0],      image (0,1/3).       (20)
                            [2 1]
```

The three open images are disjoint and partition `(0,1)` up to the two
seams.  Their inverses are respectively

```text
2-1/y,                   1/y-2,                  y/(1-2y). (21)
```

Hence every non-seam positive rational has one reduction branch.  This is
the positive meaning of the ternary tree “covering the same area”: it is a
three-cylinder PGL_2(Z) cross-section of the same rational line on which
Farey bases live.  The middle matrix has determinant `-1`; the other two
have determinant `+1`.  All three act lawfully on the same Gram-owner state
by (8), because (8) holds for all of `GL_2(Z)`.

On triples the corresponding Berggren matrices are

```text
B_A = [1 -2 2;  2 -1 2;  2 -2 3],
B_B = [1  2 2;  2  1 2;  2  2 3],
B_C = [-1 2 2; -2  1 2; -2  2 3].                        (22)
```

Direct substitution in (19) intertwines (20) with (22).  They preserve
`a^2+b^2-c^2`, have determinants `1,-1,1`, and are all congruent to the
identity modulo two.  Moreover `B_A-I` and `B_C-I` are nonzero nilpotents
of index three, while

```text
(B_B+I)(B_B^2-6B_B+I)=0.                                (23)
```

The trace of `B_B` is five, so it cannot be finite order over characteristic
zero (three roots of unity have trace of absolute value at most three).
The two nontrivial unipotents are also infinite order.  They are therefore
infinite-order reduction moves, not one order-three
action.  Their mod-two action on triple coordinates is trivial, so a
purported moving `V_4`/`S_3` label is especially unfaithful.

## 5. Two order-six quotients must not be conflated

There is a true abstract coincidence:

```text
PSL_2(Z) -> PSL_2(F_2) isomorphic to S_3,
S_4/V_4 isomorphic to S_3.                               (24)
```

The first is reduction of the modular matrices; the second is the quartic
root-pair quotient audited in candidate THM-2598.  Equation (24) identifies
an abstract quotient group, not a canonical action between the Farey,
Berggren, and quartic objects.  In particular the Berggren triple
representation is identity modulo two.  A map of targets without a map of
actions, predicates, and fibres proves no transfer.

THM-2597 supplies the other natural order-six quotient,

```text
PSL_2(Z)_ab = C_2 x C_3 isomorphic to C_6.                (25)
```

There the free factors commute and all commutator/word data are lost.  In
`S_3` they do not commute.  The fork `(25)` versus `(24)` is the exact
distinction hidden by cardinality six.

## 6. Transfer contract and stopping boundary

The useful connection is now precise.

| source | target | map | preserved | lost / required sidecar |
|---|---|---|---|---|
| Farey flank | THM-2056 defect | `U -> (U^TU,U^Tw)` | every defect value | hull-owner validity outside the fixed cone |
| modular word | basis recursion | `(G,ell)->(g^TGg,g^Tell)` | passive physical vector and defect | active Euclidean symmetry |
| Berggren branch | rational-line cylinder | (20)--(22) | rational parameter, Lorentz triple law | torsion/free-word interpretation |
| mod-two or endpoint quotient | finite label | reduction/scalarization | its stated residue or bits | Gram cross term, action, phase, owner word |

The cheapest next LRC test is consequently not another tree isomorphism.
It is an exact basis recursion whose state includes `(G,ell)` and whose
accepting predicate is THM-2056's owner-cone gate, followed by the
independent phase/owner sidecars already named there.  The entries of
`(G,ell)` are unbounded, so this is not a finite automaton without a height
or residual truncation; THM-2056's finite uncertified residual is the
natural finite table to test.  This candidate proves neither that such a
truncation closes LRC(14) nor that every useful move stays in one owner
cone.

## 7. Exact reproduction

```bash
python 04-computation/modular_farey_gram_owner_cocycle_thm2596.py
python -O 04-computation/modular_farey_gram_owner_cocycle_thm2596.py
```

Both modes reproduce the stored transcript byte-for-byte after LF
normalization.  The companion checks:

```text
S^2=-I, C^3=-I, L=-SC, R=SC^(-1);
580,000 passive Gram-owner evaluations;
the active C3 and matched-endpoint hostiles;
all Lorentz/determinant/mod-two/order checks in (22)--(23);
2,340 exact rational branch/intertwining cases.           (26)
```

All truth-bearing checks use explicit exceptions and remain active under
`python -O`.
