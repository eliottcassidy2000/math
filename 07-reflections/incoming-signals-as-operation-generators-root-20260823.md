# Incoming signals as operation generators

**Research synthesis -- 2026-08-23.**  This session promoted two exact
downstream theorems, froze one support-two operation atlas, and recorded one
honest stopping boundary.  [THM-3825](../01-canon/theorems/THM-3825-prime-colour-valuation-two-cube-decoder.md)
and [THM-3824](../01-canon/theorems/THM-3824-rule30-fixed-division-tariff-and-physical-phase-separation.md)
are **PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED**.  LRC(14),
every Rule 30 prize, and `JC(2)/DC(2)` remain **OPEN**.

The organizing rule was to treat every incoming result as a generator of a
new operation column:

```text
incoming object
  -> choose the next native operation
  -> ask whether the old equality is a congruence
  -> retain the smallest missing sidecar
  -> freeze the cheapest hostile
  -> turn the first surviving obstruction into the next task.             (1)
```

This is stronger than collecting analogies.  It makes each incoming idea
produce either a theorem, a minimal counterexample, or a typed stopping
reason.

## 1. Inheritance and live board

The session began from three incoming signals.

| lane | inherited signal | next operation | outcome |
|---|---|---|---|
| Anchor | reserved quadratic-`r`/constant-`z^2` nodal carrier | change elimination order and audit the forbidden divisor | origin promoted THM-3815 during the session; a third elimination exposed only an excluded-axis ghost |
| Niche | THM-3743's `l1<=356` support-two atlas and THM-3793's 5,855 singleton cube fibres | Stern--Brocot branching versus inert dilation | exact tree atlas plus THM-3825 valuation decoder; the operations are different and do not share one scalar law |
| Wildcard | THM-3804's all-period Smith quotient and the incoming adaptive-observer reflection | apply fixed `2^e` division, then restrict to physical sibling blocks | THM-3824 exact tariff; the physical free defect strictly separates same-scale phase |

The closest mechanisms were THM-3813's finite invariant charts,
THM-463's Eisenstein split, and THM-3804's factorization through the odd
period core.  The principal hostiles were an invalid scale normalization,
the taxicab collision `1729`, the orientation collision at cube address `28`,
and the two-slot exact-gap split `(1,1)->(3,1)`.

## 2. Three different meanings of “map it to a natural number”

The user's odd-square convention is the right starting attitude: once a
carrier has been selected, its `n`th element may be called `n` rather than
dragging the ambient value `(2n-1)^2` through every formula.  But three roles
must remain separate.

### Dense ordinal

A dense ordinal is a round-trip scheduler

```text
rho:X_N <-> {1,...,|X_N|}.                              (2)
```

It is ideal for exhaustive computation.  Its successor normally has no
mathematical meaning.

### Sparse action code

An action code is allowed to waste natural numbers, but it pulls a native
operation back to a closed law:

```text
h(Fx)=phi_F(h(x)).                                     (3)
```

Binary heap addresses are the model.  The unused integers are validity
sentinels, not a defect.

### Arithmetic carrier

An arithmetic carrier uses unique factorization, valuation, or another
ambient operation to encode coordinates:

```text
A(cx)=psi_c(A(x)).                                     (4)
```

THM-3825's tagged cube values are of this type.  They carry inert dilation
exactly but not Stern--Brocot branching.

The address contract is therefore:

```text
round trip:       decode(encode(x))=x or name the fibre;
operation:        name every F for which (3) or (4) holds;
predicate:        prove the target predicate is constant on address fibres;
gauge:            record label order, orientation and symmetry transporters;
sidecars:         retain every coordinate needed by the next consumer.     (5)
```

Overlaps are harmless when `(5)` is explicit.  A scalar can be a perfectly
good task address even when it is not an invariant or a badness classifier.

The centered triangular identity

```text
T(z+2)-T(z-2)=4z+2                                   (6)
```

illustrates the boundary.  Equation `(6)` becomes structural only after
`z->z+/-2` has been proved to be a native operation on the selected carrier.
Otherwise it is a syntactic finite difference.  In the 165-profile LRC
control, the triangular ordinal is exact but its centered difference is blind
to the second profile coordinate, so it cannot recover owner or arrival.

## 3. Support-two atlas: dense rank, tree rank, and cube rank

Let

```text
U_N={(a,b):1<=a<b, gcd(a,b)=1, a+b<=N}.                (7)
```

The shell `a+b=s` has `phi(s)/2` members, hence

```text
|U_N|=(1/2)sum_(s=3)^N phi(s).                         (8)
```

The native branches

```text
L(a,b)=(a,a+b),             R(a,b)=(b,a+b)            (9)
```

make every coprime `a<b` a rooted binary tree with root `(1,2)`.  Euclidean
subtraction gives the unique parent.  If `h(1,2)=1`, then

```text
h(Lx)=2h(x),                 h(Rx)=2h(x)+1.            (10)
```

Thus `h` is the sparse operation address.  A shell-first totient rank is the
dense ordinal; its successor is not a branch.  At `N=356` the exact atlas has

```text
19,314 primitive ratios,
maximum edge depth 353 at (1,355),
12,877 L-edges and 6,436 R-edges inside the cap.       (11)
```

THM-3825 gives a third address on the inert/cube-free subatlas.  Prime colours
decode the common inert scale, primitive shell, Eisenstein cofactor, and pair;
inert scale multiplication sends the address to its cube multiple.  A fixed
3-adic label tag gives `456,690` unordered or `913,380` oriented primitive
placements.

These structures do not collapse to one preferred rank.  In the selected
5,855-node subatlas, the Stern--Brocot induced graph has

```text
1,526 parent edges,        4,329 forest roots,
largest component size 7.                              (12)
```

The first branch hostile is

```text
(1,4), shell 5 -> (1,5), shell 6 and (4,5), shell 9;  (13)
```

neither child shell is admissible.  The cube decoder is closed under inert
dilation, not branching.  Conversely the heap address is closed under
branching but retains common scale only as a sidecar.

The full atlas also has `28` double cube fibres, beginning with

```text
1729=1^3+12^3=9^3+10^3.                               (14)
```

Thus the bare cube value is not a full-atlas address.  On the selected
subatlas it becomes injective for the precise prime-colour reason proved in
THM-3825.  The safe LRC row `(1,...,13)` still has `30` selected ratio
placements, `27` with admissible inert scale; reversibility is not loneliness.

The complete exact atlas is reproducible from

```text
04-computation/lrc_support_two_operation_address_atlas_20260823.py
05-knowledge/results/lrc_support_two_operation_address_atlas_20260823.out
```

with script/output SHA-256
`23b02703e2a4c2d33f085f80f7d021853ee50ca3410a186499fef7f0f77f874d`
and
`b0d9f82b90252e36743dbd5b9dc2b09424b9d7bd334a9e773fff80b061e7ee34`.

## 4. Rule 30: one generated task closed the next

THM-3804 says equality modulo `L=T_n^r Z^n` is a congruence for raw `T_n`,
because `T_nL<=L`.  The physical operation next divides by an adaptively
selected power of two.  Freezing one branch `e` generated the lattice

```text
K=L intersect T_n^(-1)(2^eL).                         (15)
```

The question “what information must be added?” became the exact quotient
`L/K`.  THM-3824 closes it for every `n,r,e`, including the decomposition

```text
K <= J <= L,
L/J = integral-division defect,
J/K = post-integrality carry.                         (16)
```

The hostile audit then generated a scope correction: `K` is maximal for the
fixed linear map on its divisibility domain, not for the nonlinear exact-gap
stratum.  Translation can change gap one to gap two already at `n=2,r=0`.

The next procedural task was physical realization of a nonzero tariff class.
The initial finite search found no same-scale collision, but the absence was
not left as data.  Its cause is the all-phase theorem

```text
Delta_m(t)=U_m(t+2^m)-U_m(t),
Delta_m(t+1)>Delta_m(t)>0.                             (17)
```

Packed top-bit bounds prove `(17)`.  The exact Smith free coordinate therefore
separates the entire nonnegative phase ray at a fixed scale.  This is a useful
obstruction and a useful warning: the coordinate is an unbounded encoding of
chronology, not a finite signalizer.

The next physical task is now sharply typed: search cross-scale, off-ray, or
larger-block realizations of the ambient tariff while retaining the phase
owner and exact adaptive branch.  Repeating a same-scale two-slot scan would
be redundant.

## 5. Quadratic mixed carrier: incoming proof and the normalization boundary

While this session was deriving the reserved mixed carrier, origin promoted
[THM-3815](../01-canon/theorems/THM-3815-quadratic-r-profile-constant-z2-nodal-carriers-have-critical-points.md)
with two independent eliminations.  The session therefore switched from
duplicate closure to reconciliation.

A third target-first chart uses

```text
u=re,       K=1+2u,       U=u(1+u),
g=a0+a1e+a2e^2,                                           (18)

P=gu^2-K(2e^3+ueg'),
Q=e^2K^3-729g^3U^2-162gehUK^2-216eh^3UK^3.             (19)
```

The resultant has the form

```text
Res_u(P,Q)=g^3e^4H(e),       deg H=16.                 (20)
```

The honest finite-invariant chart on `a2*h!=0` is

```text
v=a0/a2^3,       x=a1/a2^2,
y=h/a2^5,        T=a2^7.                               (21)
```

One must **not** set `a2=1`: the `h^3` terms retain `T`.  The tempting
weighted-homogeneous normalization was detected as invalid and its scratch
artifact was deleted before promotion.  In the honest chart the only common
resultant survivor lies on `y=0`, exactly the already excluded `h=0` axis.
This supplies a third elimination-order explanation of the canonical theorem,
but no new dependency was promoted because THM-3815 was already proved and
independently audited.

The reusable lesson is that boundary ghosts are often useful signals: a
resultant factor supported entirely on a previously closed parameter axis
explains why a provisional packet survives without representing a point in
the genuine cell.

## 6. Procedurally generated frontier

The session leaves the following tasks, ranked by structural value rather
than by the numerical size of their addresses.

1. **THM-3818 orientation/facet packet.**  Start from THM-3825's oriented
   `156`-channel tag.  Reconstruct the exposed facet pair and the complete
   pair-sum residue schedule.  Retain ambient row, owner, root, phase, and
   arrival; the unordered `78`-channel tag is a mandatory hostile.
2. **First-return branch grammar.**  For each of the 4,329 selected forest
   roots, find the first later admissible Stern--Brocot descendant and retain
   the skipped `L/R` word.  Test whether return gaps are bounded or form a
   renewal language.
3. **Rule 30 cross-scale tariff realization.**  Same-scale exact two-slot
   collisions are impossible by THM-3824.  Search cross-scale pairs and
   larger/off-ray blocks, keeping the exact gap branch and phase owner.
4. **Bounded observer plus section sidecar.**  Combine the `D+B` adaptive
   observer with the depth-`D+B` portrait compiler and determine the least
   off-ray section state that recovers the next selected bank, not only the
   zero-ray parent.
5. **Boundary-ghost compiler.**  On the degree-at-least-six nodal residual
   frontier, label every residual factor by the parameter axis it represents
   before launching another Groebner computation.  An excluded-axis factor is
   a stopping reason, not a new resonance.
6. **Typed task ordinal.**  Shortlex-rank only realized tuples

   ```text
   (object, representation, invariant, operation, quotient,
    symmetry gauge, scale, hostile, required sidecar).                 (22)
   ```

   Treat this rank as a scheduler.  Promote it to structure only after the
   chosen operation has a proved pulled-back law and the target predicate
   descends through every address fibre.

That is the durable interpretation of “the `n`th odd square is just `n`”:
normalize the selected carrier freely, but never confuse the normalized index
with an operation, an invariant, or a semantic certificate until those extra
laws have been proved.
