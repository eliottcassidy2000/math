# Closed area predicates and the rotation constraint on distance-graph colourings

2026-10-20. **PROVED:** a boundary-complete observer-area criterion and an
explicit counterexample to the recalled fractional-colouring identity in
[HYP-9350](../../../05-knowledge/hypotheses/HYP-9350-polygonal-order-kepler-flow-and-the-hadwiger-nelson-bridge-of-the-lonely-runner-conjecture.md).
**FINITE-EXACT:** the independent rational controls in
[boundary_bridge_check.py](boundary_bridge_check.py).
These repairs preserve the useful Kepler and graph constructions while
keeping the missing common-clock condition explicit.

## 1. One-sided area rates retain tight lonely times

On the regular n-gon, n>=3, with an observer at vertex0, let Phi_i be the
swept-area clock about the observer of a positive-speed runner. The rate
is zero on the two incident sides and positive on every other side.
At a vertex use its two one-sided time derivatives L_i^- and L_i^+.
Then the **closed** loneliness predicate is exactly

    for every i, max(L_i^-(t),L_i^+(t))>0.             (1)

At either far endpoint of an incident side, one derivative is zero and
the other positive. At the observer both are zero. In each side interior
the two coincide. This proves (1), including all endpoints. Multiplying
the positive area rate by2 gives angular momentum; positivity is the
same, but the rate itself is half that momentum.

The right-derivative-only criterion fails for speeds D={1,2}, n=3, at
t=1/3: both runners are exactly distance1/3 from the observer, but the
second enters the last incident side and has zero right area rate.
The same happens to the other runner at t=2/3. These are the only lonely
times modulo1, so discarding a measure-zero boundary discards every
witness. Random interior tests cannot check this implication.

For integer speeds, the measure of the closed lonely set agrees with
the measure after deleting vertex times. It is positive iff kappa(D)>1/n.
Zero measure means kappa(D)<=1/n: the strict inequality is the hypothetical
counterexample case and cannot be removed by assuming LRC. When equality
holds, the nonempty lonely set is finite; a finite positive-speed lower
envelope has no positive-length constant segment at a positive height.

## 2. A four-colour graph whose best rotation needs five

Let D={1,3,4,7} and let G have integer vertices, with an edge whenever
the absolute difference is in D. The vertices0,3,4,7 form a K4: their six
differences are3,4,7,1,4,3. Conversely the periodic colouring

    c(z)=2 floor(r/4)+(r mod2),       r=z mod8,

uses four colours. Differences1,3,7 change parity, while difference4
changes the first bit. Thus

    chi_f(G)=chi(G)=4.                                 (2)

For the rotation parameter kappa(D)=max_t min_(d in D)||dt||, symmetry
and period one restrict t to[0,1/2]. These intervals cover it:

| Interval | A speed whose distance is at most1/5 |
|---|---|
| [0,1/5] | 1 |
| [1/5,3/10] | 4 |
| [3/10,2/5] | 3 |
| [2/5,9/20] | 7 |
| [9/20,1/2] | 4 |

For example, in the fourth interval7t lies in[14/5,63/20], entirely
within distance1/5 of3. The other rows follow immediately in the same
way. At t=1/5 the minimum is1/5. Therefore

    kappa({1,3,4,7})=1/5,       1/kappa=5>4=chi_f(G).  (3)

This refutes the recalled equality in HYP-9350(c); it is not a refutation
of a checked theorem attributed to its named literature source. The
incoming note explicitly marked that citation as memory-only.

## 3. The exact survivor retains the group action

A circular r-colouring places adjacent vertices on a circle of circumference
r at circular separation at least1. A rotation colouring is the restricted
form z -> r{zt}, for one common t. It is valid exactly when
||dt||>=1/r for all d in D. Thus the optimum among these rotation colourings
is exactly1/kappa(D). Ordinary circular colourings need not be rotations.
The valid direction is

    chi_f(G)<=chi_c(G)<=1/kappa(D).                    (4)

For completeness, a circular r-colouring gives the fractional bound by
randomly rotating an open unit arc: the vertices it contains form an
independent set and every vertex has inclusion probability1/r. This
defines the bound on every finite subgraph, hence the usual fractional
chromatic number of the countable graph. An ordinary q-colouring is a
circular q-colouring by placing its labels at integer points of that
circle. These constructions give every direction used in(2)–(4).

In the counterexample, c(1)=1 but c(2)=0 modulo4, so the four-colouring
does not satisfy c(z+w)=c(z)+c(w). Forgetting this additive equivariance
is exactly the loss that invalidated the equality. The map from a shared
runner clock to a graph colouring is sound; its inverse needs a separate
rotation constraint.

Indeed every finite subgraph of G(Z,D) is |D|-degenerate: its smallest
integer vertex has at most |D| larger neighbours. Greedy colouring and
the usual finite-to-infinite compactness theorem give chi(G)<=|D|+1.
That elementary unrestricted-colouring bound cannot by itself prove LRC,
which asks for a rotation colouring at the same scale.

## 4. Conditional denominator and real-orbit guards

For integer D define N_min(D) as the least denominator of a lonely
rational time, or infinity when no such time exists. The periodic
piecewise-linear objective attains its maximum at a rational time, so

    N_min(D)<infinity iff kappa(D)>=1/(|D|+1).

Rational attainment alone establishes neither side. The vertex criterion
at time1/n remains valid: it works exactly when n divides no speed.

There is nevertheless an explicit finite-but-unbounded denominator family.
Fix k>=2, L>=3, put m=lcm(2,...,L), and take

    D={1,m,2m,...,(k-1)m},       t=1/2+1/(km).

At any rational time with denominator at most L, the speed m collides
with the observer. But m is even, so at the displayed t the other speeds
have phases j/k, j=1,...,k-1, whose circular distances are at least1/k.
The speed1 has distance1/2-1/(km)>=5/12>1/(k+1). Thus

    L<N_min(D)<=km

for every fixed k>=2. This elementary construction settles unboundedness
without assuming universal LRC; it does not settle a growth bound in the
largest speed. It also shows why a universal finite search cannot cap its
denominators using only the number of runners.

For real speeds, a dense orbit's closure meeting a closed safe box proves
the supremum threshold, not actual hitting without an additional argument.
An interior safe point supplies that argument by density; a closed orbit
also does by attainment. In THM-4737's rank-one case, n>=4 gives interval
length at least3(1-2/n)>1, so an integer lies in its interior and the safe
box has an interior point on the relation torus. For n=3, rank one makes
the two speeds commensurable and the orbit closed. Thus the stated
rank-one survivor remains valid after the quantifiers are made explicit.

The shared lesson with timed Collatz blocks is precise: an existence or
colouring witness is usable only for the retained source, group action,
boundary convention and clock. Each of these is a mathematical constraint,
not presentation metadata.

## 5. Reproduction

    python -B 05-knowledge/results/runner_area_clock_20261020j/boundary_bridge_check.py
    python -B -O 05-knowledge/results/runner_area_clock_20261020j/boundary_bridge_check.py --verify

The checker tests every sampled polygon endpoint exactly, compares the
one-sided predicate with rational circular distance, verifies the literal
four-colouring and clique, and independently computes kappa from all
affine-envelope candidate times. It is not a general LRC computation.
The LCM controls cover k=2,...,8 and L=3,...,12; the all-parameter argument
is the proof above. Independent proof review by `family_algebra` confirmed
the boundary, colouring, conditional-denominator, explicit LCM-family
and real-orbit repairs. Normal and optimized runs both pass11,860 exact
controls; the LF script hash and finite ranges are in the saved receipt.
