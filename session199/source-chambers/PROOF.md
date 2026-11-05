# Three source chambers recover the endpoint coordinate

Session199. **PROVED (elementary, pending independent audit):** on actual
nontrivial sources, E is an explicit piecewise-linear function of D,U,d.
The complete O/Z/J source grammar therefore has an exact nonlinear
three-coordinate realization. Away from the balanced plane, its physical
inverse has at most one candidate predecessor. Balanced resets are the only
step that needs a deletion-coordinate search. **OPEN:** all-target scalar
coverage. A triple satisfying the chamber inequalities need not be physical.

## 1. Inherited source grammar

Use the complete physical positive source grammar from Session198
`bilinear-descent/PROOF.md`, published at
`6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9`, and the canonical balanced inverse
from `balanced-source-normal-form/PROOF.md`, published at
`cf2263bb4449776d1ce13852528d5d55141c9252`. Frozen copies are in dependencies/.

For a nontrivial start-one rooted source let

    (D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)),
    (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E).

All positive coordinates are nonnegative, D and d are positive odd integers,
U>=2d, and d=e+p. The seed is tau=(1,0,2,1), with positive vector(0,0,0,1).
The complete physical moves are

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).                           (1)

Equivalently in raw coordinates,

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).                         (2)

Every source has such a word; every word is physical. This is not a grammar
on arbitrary integer points in a cone.

## 2. A strict gap in the positive f coordinate

**Lemma1.** Every actual nontrivial source satisfies

    f=0 or f>d.                                         (3)

**Proof.** The seed has f=0. A J resets f to0. A Z preserves both f and d.
At O, equation(1) gives

    f_new-d_new = a+2f+p >=1.

Indeed D=a+2e+p is odd, so a+p is a nonnegative odd integer and is at least1.
This proves the strict gap by induction on the physical word. QED.

The strictness matters. It excludes every formal band 1<=f<=d rather than
merely providing nonnegativity. No positive f value in that band is physical.

There is a second strict separator: every strong source with f=0 satisfies

    e>p, equivalently E>U-D.                            (3a)

Such a source ends in Z or J. A Z output has p=0<e=d. A J output has
 e_new-p_new=a+p>=1 by the same oddness argument. This proves(3a).

A source with E>0 ends in Z or J. Equation(1) gives f<=a for either output.
Consequently

    U<=2(D-E)<2D whenever E>0.                           (4)

If E=0, the source is the seed or ends in O, and U>=2D because O doubles a
predecessor cover count U_old>=D. Therefore

    E=0 iff U>=2D.                                      (5)

This also recovers the weak/strong distinction in the source class without
reading the endpoint coordinate.

## 3. Exact endpoint recovery

**Theorem2.** On every nontrivial physical source,

    E = 0                 if U>=2D;
    E = d                 if D+d<U<2D;
    E = D+d-U             if U<2D and U<=D+d.            (6)

**Proof.** The first case is(5). In the other cases E>0. Every generated
source lies on at least one of E=0, p=0, f=0 according to its last O,Z,J
letter. If f>0 and E>0, its last letter must be Z, hence p=0 and E=d.
Lemma1 then reads f=U-D>d, so U>D+d. If f=0, its defining equation gives
E=D+d-U, forcing U<=D+d. These alternatives are disjoint and exhaustive.
QED.

The equality U=D+d in the third case would give E=0. Thus it cannot occur
for a physical strong source in that chamber; retaining it in the formula
makes the boundary convention unambiguous on formal inputs. It is not a new
physical equality case. For example the edge has U=2D=D+d and belongs to
the first case. The balanced plane U=D has E=d.

The endpoint is therefore not an independent coordinate of the physical
profile set. This is a nonlinear, chamber-dependent recovery statement, not
a linear identity between D,E,U,d.

## 4. Exact three-coordinate physical language

Let E(D,U,d) be the piecewise expression(6). Starting from(1,2,1), define

    O(D,U,d)=(D,2U,d),
    Z(D,U,d)=(D-E(D,U,d)+U,2U-d,d),
    J(D,U,d)=(D-E(D,U,d)+U,D+U,D).                        (7)

Their orbit is exactly the projection of all nontrivial physical sources
onto(D,U,d). Recovery by(6) is a bijection between this orbit and the actual
four-coordinate profile set, with each construction word retained.

This follows by induction using(2) and Theorem2. It does not extend the
recovery theorem to all integer triples as a realizability assertion.
For instance(1,3,1) has E=0 under(6), but is not physical: D=1 forces a
transitive source, whose U is a power of two. The formal triple(13,15,7)
is another failed small inverse control.

This does not contradict the four-dimensional minimal *linear* realization
of a fixed-alphabet source readout. Formula(7) explicitly branches at the
planes U=2D and U=D+d. It is also different from the earlier infinite-indexed
three-coordinate grammar on(D,E,U), which eliminated d from its transitions.
Here the finite three-letter alphabet is retained and d is not discarded.

On the strong subset there is a complementary recovery in the other direction:

    d=E              if E<U-D;
    d=U-D+E          if E>U-D.                          (7a)

Equality cannot occur. If f>0, E=d and U-D=f>d. If f=0, (3a) gives
E>U-D and the defining equation of f gives the second line. Thus on strong
sources either of E and d can be recovered from the other three raw
coordinates; this does not identify them or give a global linear relation.

The deletion sidecar remains genuine on the full source language. There are actual profiles(9,0,18,1)
and(9,0,18,3), with words O O J O and Z O J O from tau, respectively.
They share(D,E,U) but have distinct deletion counts.

## 5. A deterministic inverse away from balance

The following is a complete arithmetic membership inverse, always requiring
recursive physical membership. Begin with positive integer D,U,d, D and d
odd, U>=2, and recovered E. Reject any failed inherited inequality

    0<=E<=d, D-E>=d, U>=D+d-E, U>=2d.

Accept tau separately. Each candidate predecessor below must satisfy the
same guards, and applying its indicated physical letter must return the
full proposed profile. These checks prevent a chamber formula from being
mistaken for a physical certificate.

### Weak chamber: one candidate

If E=0, the last letter is O. U must be even, and its only possible physical
predecessor triple is

    (D,U/2,d).                                          (8)

E of that predecessor is uniquely recovered by(6). For a source other than
tau this predecessor is nontrivial. Thus there is at most one branch.

### Strong nonbalanced f>0 chamber: one parity-selected candidate

Here E=d, U>D+d, and the last letter is Z. Put

    u=(U+d)/2, K=D-u.

Reject if u is not an integer. The only possible predecessor triples are

    (K,u,d), (K+d,u,d).                                 (9)

Indeed the predecessor satisfies D_old-E_old=K and U_old=u. If it is weak,
E_old=0 gives the first candidate. If its f>0, E_old=d gives the second.
If f_old=0, then E_old=D_old+d-u, so K=u-d, which would give D=U and
contradict the current nonbalanced hypothesis. Since d is odd, K and K+d
have opposite parity. Physical D_old must be odd. Therefore K odd selects
the first candidate, and K even selects the second. Recursive membership
and the forward check decide whether that unique candidate is physical.

### Strong nonbalanced f=0 chamber: one chamber-selected candidate

Here p=d-E>0 and the last letter is J. The predecessor has

    D_old=d, E_old=p, U_old=U-d.

It is strong since p>0. By(6), its deletion count must be either p, when
E_old=d_old, or U_old-D_old+E_old, when f_old=0. Thus the only predecessor
triples are

    (d,U-d,d-E), (d,U-d,U-d-E).                         (10)

The two alternatives are selected uniquely by the sign of D-2d. To see
this, the forward J equation gives

    D_new-2d_new=U_old-D_old-E_old.

If the old strong source has f_old>0, then p_old=0 and the strict gap gives
U_old-D_old=f_old>d_old=E_old. Thus D_new>2d_new and the first candidate
in(10) is mandatory. If the old source has f_old=0, (3a) gives the opposite
strict inequality; D_new<2d_new and the second candidate is mandatory.
Equality is also ruled out by D_new odd. Both forward strong J branches
are unimodular on their applicable triples:

    J(D,U,d)=(D+U-d,D+U,D)  when E=d,
    J(D,U,d)=(2U-d,D+U,D)   when f=0.

The determinant of each displayed 3-by-3 integer map is 1. Their applicable
physical domains and the selector above are essential; neither is a global
linear representation. As always, discard an invalid or nonphysical result.

### Balanced chamber: the only deletion search

Here U=D and E=d. The inherited canonical balanced inverse says the source
is J(O(W)) for an actual source W, possibly the singleton, with

    D_W=d, U_W=(D-d)/2.                                (11)

The case D=3,d=1 gives the singleton W and is accepted as the triangle.
Otherwise W is nontrivial, and enumerate its odd deletion count d_W over

    1<=d_W<=min(d,floor(U_W/2)).                         (12)

Recover E_W by(6) and require recursive membership of(d,U_W,d_W).
This is complete because a nontrivial source always has d_W<=D_W and
2d_W<=U_W. It is sound by the actual JO construction. Also d<=D/3 is
necessary because U_W>=D_W; violation may be rejected first.

Every accepted predecessor has smaller U. For(8) this is immediate; for Z
use U_old=(U+d)/2<U since U>=2d; for J use U_old=U-d; and(11) gives
U_W<D/2=U/2. Thus the recursion terminates. Completeness follows from the
last-letter alternatives and the canonical balanced inverse. Every nonbalanced inverse step is deterministic. The only
unbounded list of candidate deletion coordinates occurs at a balanced reset.
No polynomial input-bit-time bound or positive answer for every target follows.
The canonical balanced OJ reset is used only for membership, not to preserve
a minimum J-letter budget. For example tau.Z^k=(2^k+1,1,2^k+1,1) has a zero-J
recipe although its canonical transitive-precursor OJ recipe uses one J.

## 6. Scope and exact checks

The verifier compares the full four-coordinate physical orbit through U=2000
with the projected piecewise orbit, checks every endpoint recovery and strict
gap, and tests the inverse against the complete orbit on every guarded triple
in its declared small box. It also records both nonphysical controls and the
lost-deletion witness. Normal and optimized runs use explicit exceptions and
must agree. The all-order proof is the induction and inverse classification,
not an extrapolation of these checks.

The result makes source membership and inverse branching more explicit.
It does not classify all physical triples in a static region, prove a
balanced conductor, or establish the {7,21} spectrum conjecture.
