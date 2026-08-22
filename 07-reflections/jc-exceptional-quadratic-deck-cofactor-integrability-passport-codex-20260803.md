# The exceptional deck restores a root label but not a Keller cofactor

**Status:** structural synthesis for
[THM-3309](../01-canon/theorems/THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction.md),
on the fixed `C=c+x`, `d=k=1` slice in the two THM-3212 accessory fields.  The
quadratic base change, cofactor pullback, deck traces and norms, and
first-normal identities below are exact.  The theorem constructs
no polynomial mate, Keller inverse, `JC(2)`, or `DC(2)` consequence.

The matching companion is
[`jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.py`](../04-computation/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.py),
with frozen
[`output`](../05-knowledge/results/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.out).

## Inheritance pass and the question pulled onto the deck

The closest proved mechanism is
[THM-3306](../01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md).
It separates two PRS pivots that must remain distinct:

```text
linear_a(x):  degree 36, coefficient of the primitive linear row,
P2(x):        degree 32 before reduction, leading coefficient of the
              preceding quadratic row.                           (1)
```

At the transverse base ideal of the linear row, the
[exceptional-quadratic scout](jc-affine-c-exceptional-quadratic-c2-blowup-codex-20260803.md)
constructs

```text
A_i=K_i[x]/(linear_a),
F_0(t)=P2 t^2-Q2 t+R2,
B_i=A_i[t]/(F_0).                                      (2)
```

It proves that `F_0` is irreducible and separable.  Thus `B_i/A_i` is a
quadratic field extension.  Its degree is `2` relative to `A_i` and `72`
over `K_i`; after fixing one of the `36` embeddings of `A_i`, there are two
geometric directions, while over an algebraic closure of `K_i` there are
`72` points in `36` pairs.  This field naming is load bearing by
[MISTAKE-362](../01-canon/MISTAKES.md).

The closest positive sidecar is the
[critical inverse-graph and cofactor probe](jc-critical-inverse-cover-cofactor-jacobian-probe-agent-20260803.md).
Away from the base divisor it obtains a unique critical graph and an
elimination-cofactor identity.  Its mate-integrability class `mu(P)` is
defined only after the true gradient ideal is known to be the unit ideal.

The canonical hostile is the cofactor blindness of
[THM-3066](../01-canon/theorems/THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor.md):
symmetric cofactor shadows do not recover a sheetwise Keller equation.  The
least-used positive comparison is the inverse-different discipline of
[THM-3064](../01-canon/theorems/THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary.md),
where a genuine primitive chain-rule cofactor is supplied before its
different valuation is interpreted.

The live board was therefore:

| object | operation on `B_i` | decisive question |
|---|---|---|
| tautological `t` | apply the quadratic deck involution | does adjoining `t` restore a branch label or choose one over `A_i`? |
| elimination pair `(U,W)` | evaluate at `y=-t` | does its unit identity become a Keller Bezout row? |
| monic derivative `f_0'(t)` | invert in the etale deck | inverse different or genuine Keller cofactor? |
| true gradient `(P_x,P_z)` | pull back to `B_i` | is the first mate gate open? |
| first-normal velocity | conjugate and trace | does the deck label persist equivariantly to first order? |

The operative method card was **close a section under its next native
operation before trusting its scalar shadow**.  Here that operation is base
change to the smallest root-label carrier and then pullback of every alleged
cofactor to it.

## The root label is restored only in the pointed quadratic category

Write

```text
T=Q2/P2,                 N=R2/P2,
sigma(t)=T-t.                                               (3)
```

Then in `B_i[s]`,

```text
F_0(s)=P2(s-t)(s-sigma(t)),
t+sigma(t)=T,             t sigma(t)=N.                    (4)
```

Thus the pointed algebra `(B_i,t)` carries a tautological root and its
conjugate.  This is an exact restoration of the root coordinate erased by
the scalar resultant.  It is not an `A_i`-rational choice: forgetting the
pointing leaves the connected quadratic closed point, and `sigma` exchanges
the two labels.

Set

```text
delta=Q2^2-4P2R2,
Fprime=2P2t-Q2.                                           (5)
```

Exact reduction gives

```text
Fprime^2=delta,                  sigma(Fprime)=-Fprime.   (6)
```

The inherited nonsquare certificate makes `delta` nonzero and nonsquare in
`A_i`.  Hence `Fprime` is a unit in `B_i`, the two roots are distinct, and
there is no hidden fixed branch on which the deck action could collapse.

## Pullback of the elimination pair

Use `R_1,R_2` for the localized gradient cubics and reserve `R2` in `(2)`
for the constant coefficient of their quadratic PRS row.  Put

```text
h=P2 y^2+Q2 y+R2,
D=6P2-4Q2,
Qelim=4P2 y+D,
U=P2^2-V'Qelim,
W=-4Qelim.                                                (7)
```

The universal fraction-free identity is

```text
ell=U R_1+W R_2,
U-(V'/4)W=P2^2.                                          (8)
```

On the exceptional deck `y=-t`.  Since `Q2-P2=4V^2`, equations `(5)--(7)`
become the especially transparent invariant/anti-invariant splitting

```text
Qelim=-24V^2-2Fprime,
W=96V^2+8Fprime,
U=P2^2+24V'V^2+2V'Fprime.                               (9)
```

The script checks that `Qelim,U,W` are each nonzero, hence units, in both
quadratic fields.  Their displayed traces are

```text
Tr(Qelim)=-48V^2,
Tr(W)=192V^2,
Tr(U)=2P2^2+48V'V^2,                                   (10)
```

and their norms are

```text
N(Qelim)=4(144V^4-delta),
N(W)=16N(Qelim),
N(U)=(P2^2+24V'V^2)^2-4(V')^2delta.                    (11)
```

Every element in `(11)` is nonzero in the two exact cases.  Trace, norm,
and the invariant combination in `(8)` descend to `A_i`.  The cofactor pair
itself does not even descend projectively:

```text
U sigma(W)-sigma(U)W=-16P2^2Fprime !=0.                 (12)
```

Thus the pair remembers which exceptional root was selected.  This is a
genuine branch passport, stronger than its trace/norm shadow.

But `(8)` must be read in the right direction.  It says that the two
**coefficients** `U,W` generate the unit ideal when `P2` is a unit.  It does
not say that a combination of the gradients equals one.  On the deck,
`ell=R_1=R_2=0`, so the first identity in `(8)` reduces to `0=0`.  A
unimodular elimination-cofactor pair is not a gradient Bezout row.

## The inverse different is a unit, but it is not Keller data

For the monic exceptional polynomial

```text
f_0(s)=s^2-Ts+N,
f_0'(t)=Fprime/P2,                                      (13)
```

equation `(6)` gives

```text
Tr(f_0')=0,
N(f_0')=-delta/P2^2.                                    (14)
```

Consequently the quadratic inverse different has the explicit unit

```text
eta=1/f_0'(t)=P2/Fprime,
sigma(eta)=-eta,
N(eta)=-P2^2/delta.                                     (15)
```

This is the cleanest possible etale primitive-element packet.  It answers
one part of the pullback question positively: the deck derivative and its
inverse remain units, and their anti-invariant parity is exact.

It does **not** answer the Keller-cofactor question positively.  THM-3064's
primitive Keller cofactor is chain-rule data supplied by an inverse-spectral
pair, schematically `q/f'`.  Here there is no supplied `q`, inverse map, or
second coordinate.  Equation `(15)` is the special algebraic reciprocal
`1/f_0'`; calling it a Keller cofactor would silently replace missing
chain-rule data by the unit numerator `1`.

Three objects carrying the word “cofactor” are therefore sharply distinct:

```text
(U,W):       fraction-free elimination coefficients;
eta:         inverse-different generator of the quadratic root deck;
q/f':        primitive Keller cofactor, absent here.                   (16)
```

## The mate pipeline fails before the divergence class

The decisive calculation is not a trace or norm.  For

```text
P(x,z)=(V(x)z^2+z+C(x))^2+A(x)z+E(x),
y=Vz,                                                      (17)
```

the localized cubics satisfy the exact change of gradient coordinates

```text
R_1=V P_z,
R_2=V^3P_x-(V'y/2)R_1.                                  (18)
```

The script directly evaluates both cubics at `y=-t` and obtains

```text
R_1(-t)=R_2(-t)=0 in B_i.                               (19)
```

THM-3306's owner separation makes `V` a unit on this base locus, so `(18)`
can be inverted.  Hence

```text
P_x=P_z=0 in B_i.                                       (20)
```

Since `B_i` is a nonzero field, the evaluation map to `B_i` witnesses

```text
1 notin (P_x,P_z).                                      (21)
```

For every polynomial `Q`, independently of its coefficients,

```text
Jac(P,Q)=P_xQ_z-P_zQ_x=0 in B_i.                        (22)
```

This identifies the first failed implication exactly.  The canonical class

```text
mu(P)=[A_x+B_z] in coker(P_x partial_z-P_z partial_x)
```

from the inverse-graph reflection presupposes a row `AP_x+BP_z=1`.  Equation
`(21)` proves that no such row exists on this slice.  Therefore `mu(P)` is
**undefined here**, not merely uncomputed and not proved nonzero.  The mate
pipeline stops at gradient unimodularity, before the divergence-class gate.

This also answers the “does a primitive Keller cofactor remain unit?” query
directionally:

- the deck inverse-different unit `eta` exists;
- no primitive Keller cofactor was supplied to pull back; and
- any true physical Jacobian pulls back to zero by `(22)`, never to a unit.

## First-normal compatibility of the pointed deck

Let the fixed blow-up normalization give

```text
F(u,t)=F_0(t)+u(F12 t^2+F11 t+F10)+O(u^2),
dot(t)=-F_1(t)/Fprime.                                  (23)
```

If `T(u)=Q2(u)/P2(u)` is the monic quadratic trace, then

```text
dot(T)=(-P2 F11-Q2 F12)/P2^2.                           (24)
```

The exact quotient computation in both fields gives

```text
Tr_B/A(dot(t))=dot(T),
sigma(dot(t))=dot(T)-dot(t).                            (25)
```

Equation `(25)` is precisely the derivative of
`sigma_u(t(u))=T(u)-t(u)`.  Thus the tautological label and its conjugate do
not collide or lose their involution at first order.  They continue as an
equivariant pair.  The velocity digests agree with the independent
first-normal packet in the exceptional-quadratic scout.

This is local first-order compatibility in the fixed normalization.  It is
not a global parameter-space monodromy computation and does not classify
higher normal coefficients.

## Connection contract and procedural reframe

```text
source:      transverse linear-row base ideal plus quadratic PRS row;
target:      pointed etale deck (B_i,t) with elimination and different data;
map:         base change A_i -> B_i followed by y=-t;
preserved:   both common critical roots, deck exchange, first-normal motion;
restored:    one tautological root coordinate after adjoining the sidecar;
destroyed:   descent of an individual branch and of the projective (U,W) pair;
absent:      chain-rule numerator q, gradient Bezout row, mate Q;
test:        Fprime^2=delta, projective defect (12), gradient pullback (20).
```

The useful reframe is that a resolving cover can be a **microscope rather
than a repair**.  Passing to `B_i` closes the root-selection problem exactly,
but it exposes rather than removes the earlier Keller obstruction.  The
correct operation order is

```text
scalar critical shadow
 -> quadratic root deck
 -> branchwise elimination passport
 -> test the true gradient ideal
 -> only if unimodular, form mu(P)
 -> only if mu(P)=0, reconstruct a mate.                (26)
```

On this slice the first three arrows succeed and the fourth fails.  Trace or
norm calculations after that failure cannot reopen the mate pipeline.

## Reproduction and honest boundary

Run from the repository root:

```text
python3 04-computation/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.py
python3 -O 04-computation/jc_exceptional_quadratic_deck_cofactor_integrability_scout_20260803.py
```

The companion reconstructs both response fields internally, pins the exact
THM-3306 and two prior sidecars, uses explicit quadratic-extension arithmetic,
and has no Python assertions, floating literals, or imported execution path.
The two modes must agree with the stored transcript.

```text
script SHA-256:    c9809b5cafc5defca5aeda49bdb321f9c2d20b57405020c76bebc78e4e5dd2c6
LF output SHA-256: 5026f0f17d0ff110d7f916a5b12858c04e09e2949209c7a75f6ab3e9cdc37b58
```

The result is sharply bounded:

- only the two fixed accessory fields and `C=c+x,d=k=1` are covered;
- `linear_a` has degree `36`, while `P2` has degree `32` before reduction;
- the cover has relative degree `2` and total degree `72`, as required by
  MISTAKE-362;
- the anti-invariant packet is algebraic deck data, not global monodromy;
- the result classifies neither the degree-119 residual component nor a
  deformation in `d,k,V,A`;
- no primitive Keller cofactor, polynomial mate, inverse map, `JC(2)`, or
  `DC(2)` follows.

The cheapest genuinely new deformation is now clear: release one clutch
parameter while carrying all three typed packets in `(16)`.  Test whether
the critical deck disappears before spending effort on `mu(P)`; if it does
not, the branchwise passport may describe the obstruction more finely, but
it cannot become a mate certificate.
