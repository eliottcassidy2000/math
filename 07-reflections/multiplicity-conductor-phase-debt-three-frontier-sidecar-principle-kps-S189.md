# Multiplicity, conductor, and phase debt: one sidecar principle on three frontiers

Status: **TYPED SYNTHESIS + CORRECTED FINITE-EXACT CROSS-FRONTIER EXTENSION**
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

THM-3647 compresses the static observable further.  Point-order reversal

```text
J=(05)(14)(23)                                          (11a)
```

pairs the branch addresses by `r -> -r-1`.  For every one of its seven orbit
packets `B_r`, including the fixed singleton `B_6=A_6`, the exact spectrum is

```text
B_r=alpha_r Pi_W+beta_r Pi_M,       alpha_r!=beta_r.    (11b)
```

Hence any single packet gives

```text
Pi_W=(B_r-beta_r I)/(alpha_r-beta_r),
Pi_K=Pi_W(I+S)/2.                                       (11c)
```

The fixed branch `A_6` alone therefore recovers the endpoint/middle split;
native directed-arc reversal `S=(01)(23)(45)` then recovers rigidity.  The
distinction between `J` and `S` is load-bearing.  This removes thirteen-way
spectral overhead, but not the missing physical type: an address operator is
still not a child transition or clock edge.

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
changing coordinates.

THM-3638/3639 then remove the apparent cell dependence.  For **every**
zero-stable/first-stable restriction cell on the fixed minimal Hermite fold,
including both retained tangent ranks, the division-free identity is

```text
lambda(J_2)=-2072/81+(5/27)J_0'(-1)-(13/27)J_0'(1).    (9b)
```

Thus the single polynomial identity `J_0=1` already forces the nonzero debt;
`J_1=0` is not needed.  The live frontier has moved from another cell on the
same fold to another collision polynomial or compiler.

An exact `Q_h` restriction-ring probe makes that move nonvacuous.  Its actual
restriction ring has generator degrees `24,27,39,119,134,214`, quotient
dimension `155`, semigroup conductor `299`, and global conductor degree `310`.
For `(U,V)=(c,e)` it nevertheless contains a unique first actual-ring Bezout
witness at symmetric source cutoff `334`, with `deg(T,A,B)=(308,334,331)`.
THM-3627 still closes the same hostile at source degree six for arbitrary
target two-forms.  Hence completion debt, conductor debt, zero-stable
determinant debt, and transverse stable-order debt are genuinely distinct;
the compiler must retain them as separate coordinates.

THM-3641 now supplies the missing curvature coordinate.  On every ordinary
retained cell of the full projective first-jet atlas it proves

```text
lambda(J_2)=D_Q+mu_-J_0'(-1)+mu_0J_0'(0)+mu_+J_0'(1). (9c)
```

In the principal chart,

```text
D_Q=-2((9-4u)r_-+(9+4u)r_++243)/81.                  (9d)
```

For `Q_h` (the degree-six polynomial called `Q_6` in THM-3642),
`(u,r_-,r_+)=(1,0,-243/13)`, so the retained second-order debt is exactly
zero.  This is a qualitative boundary rather than a smaller nonzero invoice:
the universal three-value obstruction that closed the minimal Hermite fold
has disappeared.  Promoted and independently audited THM-3642 reaches the
stronger global identity.  After a global actual-ring `J_1=0` lift, `J_2=0`
reduces to

```text
c'W-e'Z=r_2.                                            (9e)
```

In the `310`-dimensional quotient target, the correction map has rank `309`
and its augmentation by `r_2` still has rank `309`; exact lifts have
`deg(T_2,Z,W)=(308,1044,1041)`.  The unique cokernel is

```text
Lambda([P],[Q])=lambda(c'Q-e'P),
lambda(f)=(5f(-1)-18f(0)+13f(1))/18.                  (9f)
```

Since every representative gauge from `ker gamma_h` has vertical derivative
vanishing at the three retained points, it lies in `ker Lambda=im M`.  Thus
the global `J_2` solve is gauge-complete for every actual-ring witness in this
cell that reaches `J_1=0`.

The apparent `J_3` frontier is now superseded.  THM-3642 proves, for arbitrary
target two-forms on this fixed compiler, the universal fourth-stable identity

```text
J_0=1 and J_2=0  ==>  lambda(J_4)=365888/6561 !=0.    (9g)
```

Neither `J_1` nor `J_3` is needed.  A symmetric-curvature degree-seven
control on the same zero-second-debt wall similarly forces
`lambda(J_4)=5440/81`.  Thus zero second debt is a delay, not a passage:
the fixed degree-six actual lift reaches `J_1=J_2=0` but cannot have constant
source Jacobian.  The lawful frontier is no longer another jet on these two
folds; it is a different collision polynomial, nonquadratic fold, or compiler.

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
asymptotic law.

THM-3644 now resolves the immediate local offset question without inventing a
global threshold.  At `R=16384` the complete consecutive bracket is

```text
D0=401,...,412 : DIE,
D0=413,...,416 : CLOSED,
OPEN_RESIDUAL   : none.                                (13)
```

The adjacent event jumps from death at row `4116` for `D0=412` to closure at
capture row `10116` for `D0=413`; an independent FLINT/fmpz reconstruction
matches the load-bearing pair.  Thus `413` is the first closure **within this
bracket**, not a proved global least offset: Rule-A closure monotonicity below
`401` is absent from the proof graph.  The discontinuous event coordinate in
`(13)` is another instance of the sidecar principle.  An offset records the
input word, but the terminal branch and capture/death row record different
fibres of the exact automaton.

## 5. Two-cube and Mordell arithmetic: support and Selmer sidecars

The Berggren/two-distinct-cubes lane exhibits the same failure of scalar
compression in arithmetic form.  For a primitive parity-correct slope set

```text
U=n-m,                V=2m-n,                T=2m+n.   (14)
```

THM-3645, now independently hostile-audited, proves that every odd mod-prime
obstruction to the associated conic is supported on `nVT`.  At a support
prime its entire mod-prime content is one of two characters:

```text
p|n or p|V : (6/p)=1,                 p|T : (-2/p)=1. (15)
```

Together with the `p=3` condition, a screen through `n` is complete because
`V<n` and `T/3<n`.  It forces exactly the two cones

```text
(n,V)=(5,23) or (23,5) mod 24.                         (16)
```

This replaces a long blind prime screen by a typed support divisor plus two
quadratic characters.  It also exposed a concrete correction: the first
slope passing the fixed `p<=997` screen but failing the complete support test
is `(m,n)=(512,1019)`, already obstructed at `p=1019`; the formerly reported
`(1012,1039)` is only a later escape.  Neither `(15)` nor `(16)` supplies
`p`-adic solubility, an integral cube decomposition, or an admissible Pell
orbit.

The exact S190 continuation counts the complete mod-prime gate through
`n<=49999`: among `126,652,918` primitive parity candidates there are
`633,416` survivors.  The normalized values
`A(N)(log N)^(3/2)/N^2` remain near `0.009`, motivating—but not proving—the
Frobenian three-semigroup law

```text
A(N) ~ kappa N^2/(log N)^(3/2).                        (16a)
```

This is a fractional-sieve conjecture, not a Pell or global-solubility count.

On the fixed-`107` Mordell curve

```text
E: y^2=x^3+1225041,                                    (17)
```

THM-3643 now independently certifies the ordinary real-quadratic class number
`h(Q(sqrt(1225041)))=1`.  Independently audited THM-3646 proves rank exactly
three: a rational point with denominator `14723^2` adds a third independent
mod-`3` class, while explicit `3`-isogeny support bounds the rank above by
three.  The conceptual correction is important:
`3`-saturation of the old two-point subgroup did not imply completeness of
the Mordell--Weil group.  The missing sidecar was the complementary isogeny
Selmer support, not a denser search inside the old lattice.

An exact bounded continuation scans all `68,920` nonzero combinations
`aP+bQ+cR` with `|a|,|b|,|c|<=20`.  Its only integral absolute points are the
two known points `P,Q`; hence it finds no new fixed-`107` two-cube depth.  The
companion paths and hashes are

```text
04-computation/berggren_fixed107_rank3_integral_lattice_box20_kps.py
fcda4363faa7d14275e692ffc7e78c3f80b0a811b3652cf0a2ab5013fe8374f8

05-knowledge/results/berggren_fixed107_rank3_integral_lattice_box20_kps.out
4e1b6847467c62ab12f30293386e902cb934f06787b40514b9dca821b4a01efa.
```

This is only a coefficient box in `<P,Q,R>`.  Until the subgroup index and a
height reduction are known, it is not an integral-point classification.

## 6. Transfer rule and counterindications

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
interpolant is not algebraic, an exact finite-scale phase invoice is not a
limit theorem, and a mod-prime support classification is not a rational or
integral point.

The productive research question is therefore:

```text
What is the smallest noncommuting or nonlocal observable that separates the
remaining fibre while preserving the actual target predicate?              (18)
```

That question is now concrete on five frontiers: a lawful digit/entry action
realizing the fixed packet `A_6` and arc involution for LRC; a new
collision polynomial, nonquadratic fold, or compiler beyond THM-3642's
universal fourth-order wall; offset monotonicity versus terminal-event type
at `R=16384` for AMM; a
`p`-adic/global sidecar beyond the two-cube Frobenian support characters; and
the precise index/height structure behind the new third Mordell point.
