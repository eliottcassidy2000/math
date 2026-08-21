# Multiplicity, conductor, and phase debt: one sidecar principle on three frontiers

Status: **TYPED SYNTHESIS + CORRECTED FINITE-EXACT LRC EXTENSION**
(2026-08-21).
This is not a theorem and does not create a dependency.  It specializes the
existing meta-patterns “controlled forgetting requires a sidecar” and “the
same representation is not the same carrier.”

## 1. The common diagnostic

Suppose an exact map `q:X->Y` retains a quotient observable but the desired
predicate asks for a distinguished object inside a fibre of `q`.  Three debts
must be separated:

```text
multiplicity debt: several lines carry the same quotient character;
conductor debt:    separate completed branches do not descend to one function;
phase debt:        a normalized scalar forgets the intercept/cocycle of a word.
```

The useful next object is not another polynomial in the quotient observable.
It is the smallest typed sidecar which acts inside the unresolved fibre.  The
cheap tests are correspondingly different:

1. compute the commutant and repeated-character spaces;
2. compute the evaluation/conductor equalizer across branches; or
3. write an exact defect invoice before fitting a limiting scalar.

The sidecar must preserve the target predicate and its lost information must
be stated.  A larger algebra, a formal interpolation, or a better numerical
fit alone is not a completion.

## 2. LRC: a noncommuting arc involution pays multiplicity debt

THM-3625 reconstructs the thirteen point-by-branch operators `A_r` on the
six-point module over one pinned finite field.  They generate

```text
mathcal A isomorphic to F_p^4,
character multiplicities (2,1,1,2).                    (1)
```

The inherited four-space `W=P6sharp intersect Q6sharp` is the sum of the two
doubled characters, while `K2sharp` takes one line in each.  Every element of
`mathcal A` is scalar on a character block, so no further polynomial in the
branch addresses can select `K2sharp` as a kernel or image.

The least-used sidecar is already native to the pointed `P4` carrier: reverse
each directed arc,

```text
S=(01)(23)(45).                                         (2)
```

An exact independent probe (pending theorem-level hostile audit) gives, in the
pinned point frame,

```text
W = span(e0,e1,e4,e5),
K2sharp = span(e0+e1,e4+e5),
K2sharp = W intersect ker(S-I).                         (3)
```

Moreover the address algebra already contains the projector

```text
Pi_W=diag(1,1,0,0,1,1),                                (4)
```

whereas adjoining `S` enlarges it to a noncommutative algebra of provisional
dimension eight.  Hence

```text
Pi_K=Pi_W (I+S)/2                                      (5)
```

is the candidate rank-two projector onto `K2sharp`.  The gain is precise:
`S` splits the multiplicity-two blocks.  The loss is also precise:
`(I+S)/2` forgets directed-arc orientation, and `(4)` discards the middle arc.
Neither operation supplies chronology, a current, characteristic-zero
transfer, or an LRC(14) row exclusion.

An older, independently audited third-digit carrier gives a striking hostile
alignment without licensing a transfer.  There the same arc reversal fails on
the full `78`-coordinate parent, and its entire rank-two defect is localized on
the middle point pair `(2,3)`; the endpoint pairs `(0,1)` and `(4,5)` have zero
defect.  Stabilizing by the involution yields an eight-dimensional quotient
with balanced arc-even/odd characters, but none of the thirteen child spaces
is itself arc-invariant.  Thus `(4)` removes exactly the historical defect
location, while the child failure proves that the static projector is not yet
a chronological transition.  This is evidence of the mechanism, not an
identification of the two carriers.

The first tempting quotient map was the endpoint-arc imbalance on the centered
coefficient frame,

```text
delta_arc:R4sharp -> F_p^2,
delta_arc(c)=(c0-c1,c4-c5).                             (6)
```

The hostile reconstruction **refutes** this interpretation.  In the same
centered point frame,

```text
K2sharp=span(e0+e1,e4+e5),
R4sharp=span(e0+e1,e2,e3,e4+e5).                       (7)
```

Hence `delta_arc` vanishes on all of `R4sharp`; its kernel is four-dimensional,
not `K2sharp`, and the augmentation cannot factor through it.  The corrected
map is the middle-arc coordinate

```text
delta_mid(c)=(c2,c3).                                   (8)
```

It has rank two and kernel `K2sharp`.  Better, it is already address-theoretic:

```text
Pi_M=I-Pi_W=diag(0,0,1,1,0,0) in mathcal A,
ker(Pi_M|R4sharp)=K2sharp,
im(Pi_M|R4sharp)=span(e2,e3).                           (9)
```

The THM-3615 augmentation now has the exact factorization

```text
LambdaSharp|R4sharp=T_mid o delta_mid,                  (10)
```

where, in the pinned middle-coordinate and canonical augmentation bases,

```text
T_mid=
[126498113787680818370196  646561219255993169342961]
[384727693242765231857657  452865069702498262026030],

det(T_mid)=275730381649850587765623 !=0 mod p.          (11)
```

Thus the intrinsic quotient isomorphism
`R4sharp/K2sharp -> A2sharp` is real, but the endpoint story was false.  The
displayed matrix depends on the pinned bases; the quotient map does not.  It
still needs a temporal-entry sidecar before any physical current claim.

## 3. Jacobian fold: completion pays jets but not diagonal constants

THM-3630 constructs a decomposable formal pair on a product of three completed
source germs and polynomial approximants through every prescribed finite
cutoff.  Newton interpolation succeeds because each disconnected branch may
carry its own integration constant:

```text
f_i=Cx+h_i(w).                                          (7)
```

THM-3632 identifies the conductor debt.  If one target-local algebraic
function `F` paired with the stable coordinate `w`, its pullback on the
connected source would satisfy

```text
partial_x(Phi^*F)=C,
Phi^*F=Cx+H(w)                                          (8)
```

with one common `H`.  A retained same-`w` collision `x_1!=x_2` maps both
points to the same target value, contradicting `(8)` when `C!=0`.  The formal
construction is not false: completion replaces one diagonal constant by three
branch constants.  It is precisely nonalgebraic in the target local ring.

Thus the relevant quotient is

```text
one connected rational function
    -> product of completed branch functions,           (9)
```

and the missing sidecar is the diagonal residue/evaluation condition on the
collision fibre.  More collision jets cannot recover it uniformly: THM-3630
proves survival at every finite cutoff.

Incoming THM-3635/3637 sharpen this boundary twice.  For the minimal Hermite
curve, the actual restriction ring has global conductor degree `82` and does
contain a genuine rank-two/rank-two zero-stable determinant witness with fixed
pair `(U,V)=(c,e)`; the least symmetric coefficient cutoff is `94`.  That
witness even admits an exact first-stable lift `J_1=0`.  Nevertheless every
unbounded determinant witness in the same restriction classes which reaches
`J_1=0` has the identical second-stable retained debt

```text
lambda(J_2)=-2072/81.                                  (9a)
```

The value is independent of the full affine first-jet gauge, higher lifts,
and representatives modulo the restriction kernel.  The vertical compiler
terms are load-bearing: deleting them changes the quotient rather than merely
changing coordinates.  Thus completion debt, actual-ring conductor debt, and
stable-order debt are three different obstructions.  The fixed `(c,e)` cell
is now closed at second order; other zero-stable pairs and the tangent-rank-one
stratum remain open.

## 4. AMM: the horizon ratio hides a three-term phase cocycle

At a failed Rule-A trace the adjoint horizon `q` is exact, but the scalar
`q/R` forgets how the Sturmian word meets the clamp front.  For the audited
dyadic scales define

```text
h=5R/8-d0,
b=2d0-R+2,
m=j-b,
ell=j-s,
e=ell-((sqrt(5)-1)/8)R,
theta=(3-sqrt(5))/8.                                   (10)
```

Then the exact invoice is

```text
q-theta R=m+1-2h-e.                                    (11)
```

For the fixed failed traces at `R=8192,16384`, promoted THM-3633 gives

```text
(h,m,ell)_8192  =(31,57,1271),
(h,m,ell)_16384 =(43,43,2538),

q_16384-2q_8192=-30
                  =(-71-1)-2(-19)-(-4).                (12)
```

The normalized error deteriorates by exactly `15/8192` at the larger fixed
trace.  Equation `(12)` explains why a one-coordinate golden fit can rebound
and then worsen: the margin collapse, headroom defect, and Sturmian depth phase
move separately.  The sidecar is the phase invoice, not a replacement
asymptotic law.  The `R=16384,D0=400` trace is known to fail but is not known
terminal in this package.

## 5. Transfer rule and counterindications

When an exact quotient almost selects the desired object:

1. name the quotient fibre and compute its dimension/multiplicity;
2. test whether the available algebra is scalar on that fibre;
3. seek a native involution, residue, intercept, or conductor map which acts
   nontrivially there;
4. express the target as an intersection, kernel, image, or diagonal equalizer;
5. audit what the new sidecar destroys and whether the physical action exists.

This move is counterindicated when the target predicate is already invariant
on quotient fibres, or when the proposed sidecar is only a coordinate gauge.
In particular, a static arc projector is not a current, a completed Newton
interpolant is not algebraic, and an exact finite-scale phase invoice is not a
limit theorem.

The productive research question is therefore:

```text
What is the smallest noncommuting or nonlocal observable that separates the
remaining fibre while preserving the actual target predicate?              (13)
```

That question is now concrete on all three frontiers: temporal entry after
the corrected middle projector `(9)` for LRC, a zero-stable pair outside the
second-order-closed `(c,e)` Russell cell (especially the tangent-rank-one
stratum), and the terminal-offset/front-phase transition at `R=16384` for AMM.
