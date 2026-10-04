# Collatz word geometry needs a marked axis, an ordered carry, and a source

2026-10-04. **PROVED** elementary affine, word, and rational-score statements
below; **FINITE-EXACT** for the declared controls. Collatz convergence and
classification of all signed cycles remain **OPEN**. This is a carrier and
obstruction result, not new integer coverage or a triangle-group conjugacy.

The useful geometric object is the actual affine word map with its finite
fixed point and signed multiplier. Its finite fixed point specifies a
vertical hyperbolic axis. Two nonempty word maps have the same axis exactly
when their marked words repeat one common primitive word. Different route
patterns therefore require an ordered carry when they are composed.

## 1. Recovery and the live boundary

The [September 21 affine blueprint audit](collatz_blueprint_20260921_affine.md),
sections 2--3, already proves that the unrestricted inverse-branch group is
an affine, torsion-free, nondiscrete group; the proposed triangle-group
identification is false. It proves a free guarded word semigroup and gives
equal-spectrum maps fixing 1 and 5/7 as the carry-loss hostile. The matching
`01-canon/MISTAKES.md` entry is titled **2026-09-21 Collatz blueprint -- invalid
geometry and missing trajectory quantifiers**. It overrides the old proposed
elliptic/parabolic classifications.

[The mediant-tree note](collatz_mediant_tree_K_20260930.md), Theorems M1--M2,
already identifies the prefix fixed point with the descent threshold and
gives the source-displacement identity. The
[periodic-chart separation proof](periodic_chart_separation_20261004.md)
already uses marked primitive periods and elementary period cancellation;
rotation and repetition are different quotients. The
[Pillai/clock census](collatz_mod6_20260921_pillai_convergents_cycle_gates.md)
already finds the seven integer anchors on the 11/7 clock. No priority claim
is made for these facts or their geometric packaging here.

The newly inherited [boundary compiler](collatz_boundary_compiler_20261004.md)
proves its sharp all-length carry bound and resolves the old one-copy guard
problem. This note imports that theorem as stated; its elementary geometric
results do not require the compiler's cited linear-forms input.

Closest mechanism: affine source displacement. Canonical hostile: moving
checkpoint `27 -> 41 -> 31`, where the last decrease does not repay 27.
Corrected near miss: projective spectral type does not classify integer
closure or signed drift. Least-used sidecar: the position of the finite
fixed point. The board is **word / axis / signed multiplier / source guard /
ordered composition / curvature score**. The anchor is a faithful Collatz
certificate interface; the niche is exact geometric information loss; the
wildcard is the proposed 5--6--7 curvature comparison.

## 2. Every nonempty positive valuation word is hyperbolic

Let `w=(a_1,...,a_j)` be a nonempty word of positive integers, with

    A=sum a_i, P=3^j, Q=2^A,
    B=sum_(i=0)^(j-1) 3^(j-1-i) 2^(a_1+...+a_i).

The formal affine map is

    F_w(z)=(Pz+B)/Q = r_w z+t_w,
    r_w=P/Q>0, t_w=B/Q>0.

It preserves the upper-half-plane metric `|dz|/Im(z)`: numerator and
imaginary coordinate both multiply by `r_w`. Unique prime factorization
gives `P!=Q`. Thus its two boundary fixed points are infinity and

    a_w=B/(Q-P),

and its invariant axis is the vertical line `Re(z)=a_w`. In the usual
determinant-normalized projective classification,

    J=tr(M_w)^2/det(M_w)
     =(P+Q)^2/(PQ)=4+(P-Q)^2/(PQ)>4,
    M_w=[[P,B],[0,Q]].                              (1)

Both expanding and contracting words are therefore hyperbolic. No nonempty
positive valuation word is parabolic or elliptic. On its axis the signed
vertical translation is `log r_w=j log3-A log2`; its length is the absolute
value. The sign is a drift coordinate, not the elliptic/parabolic/hyperbolic
type. There is no nonempty exactly neutral word, since `3^j!=2^A`.

For an **exactly guarded integer source** n,

    F_w(n)-n=(r_w-1)(n-a_w).                        (2)

For growth words `r_w>1`, the axis endpoint is negative, so every positive
legal source grows. For contracting words `r_w<1`, the endpoint is positive:
a source descends exactly when it lies above `a_w`. The metric length alone
forgets that position. The exact source cylinder is

    n=(Q-B)P^(-1) mod2Q.                            (3)

For example `(1,2)` and `(2,1)` have the same `P=9,Q=8` and spectral
invariant, but have carries 5 and 7 and fixed points -5 and -7. These are
different marked phases of the known negative cycle, not two legal choices
at one integer. Their exact cylinders must stay separate.

An arbitrary contracting word can still increase a small source. The
incoming compiler's exact hostile is

    w=(4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3),
    165 -> 167 after its 17 odd steps,
    a_w=1106233681/(134217728-129140163)>165.

Its first step is already `165 ->31`, so this is not a counterexample to
the open first-stopping coincidence. It tests precisely the missing
source-versus-axis inequality in a general word.

## 3. Same axis means one primitive pattern; different axes retain order

Write words chronologically: `uv` means perform u and then v. Direct affine
composition gives

    F_(uv)=F_v o F_u,
    (P,Q,B)_(uv)=(P_v P_u, Q_v Q_u, P_v B_u+B_v Q_u).

**PROVED order defect.** For any nonempty words u,v,

    F_v(F_u(z))-F_u(F_v(z))
      =(r_u-1)(r_v-1)(a_v-a_u).                    (4)

The right side is constant. Substituting `t_w=(1-r_w)a_w` proves (4).
Since neither multiplier equals one, the maps commute if and only if their
**finite fixed points** agree. Merely sharing infinity says nothing: every
map considered here shares it.

The representation of positive valuation words by affine maps is injective.
The reduced multiplier `3^j/2^A` recovers j and A. It remains to decode B.
When j>=2, split off the first letter:

    B=3^(j-1)+2^a_1 B_tail,

where B_tail is the positive odd carry of the remaining word. Consequently
`a_1=v2(B-3^(j-1))`, after which B_tail is recovered by division. Repeat,
and recover the final exponent from A. The final carry is 1. This is a
direct decoder for the inherited free-semigroup fact.

It follows that the following are equivalent:

1. `a_u=a_v`;
2. `F_u` and `F_v` commute;
3. the literal marked words satisfy `uv=vu`;
4. `u=s^m, v=s^n` for one common marked primitive word s and positive m,n.

For the last step, if `|u|<=|v|`, word equality forces v to begin with u,
say `v=ut`. Cancel the first u from `uv=vu` to obtain `ut=tu`. Induct on
the sum of lengths; stripping powers gives a common primitive root. Equal
lengths force equality. The converse follows by concatenating powers.

Thus **axis plus multiplier** is a lossless geometric encoding of a supplied
word. Axis alone forgets its positive repetition count. Passing to a cyclic
necklace loses a different coordinate: its marked phase. In fact if
`w=(a_1,...,a_j)`, then

    F_(a_1)(a_w)=a_(a_2,...,a_j,a_1),                (5)

as follows by conjugating the word map by its first letter. These rational
axis endpoints need not be integer cycle points.

For letters u=(1),v=(2), their axes are -1,+1. Equation (4) gives -1/4.
In the unrestricted real affine group,

    (F_v o F_u) o (F_u o F_v)^(-1): z -> z-1/4

is parabolic. This does not contradict (1): this commutator uses inverses
and is not a nonempty positive word, nor a legal integer Collatz path. It
shows exactly where a parabolic map can appear when one changes the type
of permitted operations.

The practical rule is narrow: macros about one fixed axis are repetitions
of one primitive pattern. Genuinely different patterns require the ordered
carry or equivalent axis-change data. This is not an obstruction to mixing
macros; the existing affine compiler already retains the required data.

## 4. Why fixed rational curvature weights cannot encode all drift signs

In an equilateral triangulation the familiar degree-five and degree-seven
vertex defects have opposite equal weights, with degree six flat. To test
an explicit candidate transfer, assign weights +1 and -1 to valuation
letters 1 and 2. This is a **formal score** on words; no claim is made that
every such word is the degree sequence of a closed triangulated surface.

The shortest flat score is `(1,2)`, but its multiplier is `9/8>1`; the
exact legal source `11 ->17 ->13` grows. The first strict wrong sign among
all binary-letter counts occurs at `j=7,A=11`: three 1-letters and four
2-letters give score -1, while

    3^7-2^11=139>0.

This is the same clock as the negative 17-cycle, not a new cycle theorem.
The exact census has 210 words with `j=7,A=11`, 210 different marked axes,
and just seven integer axes: `-91,-61,-55,-41,-37,-25,-17`. The other axes
are rational nonintegers. The clock and a curvature sign cannot choose
among them.

**PROVED general obstruction.** No additive rational score on letters 1 and
2 has the correct signed drift for all words. If its letter weights do not
satisfy `c_1>0>c_2`, a one-letter control already fails. Otherwise multiply
by a common denominator to obtain integer weights. With
`g=gcd(c_1,-c_2)`, counts `n_1=-c_2/g,n_2=c_1/g` make the score zero, while
the drift cannot be zero by unique prime factorization.

Even allowing an arbitrary fixed tie convention does not repair the score.
At a proportion x of 2-letters its sign changes at the rational threshold

    theta=c_1/(c_1-c_2),

whereas the genuine sign changes at `beta=log_2(3)-1`, which is irrational.
Choose a rational proportion strictly between theta and beta and clear its
denominator. The resulting word has strictly opposite score and drift
signs, so ties play no role. Word order is irrelevant to this particular
drift test, though it remains essential to the carry and legality.

This excludes fixed rational local curvature scores as complete drift
classifiers, not geometry in general. A carrier with the actual irrational
weights `log3-a log2`, or the exact integer pair `(3^j,2^A)`, retains the
drift perfectly. It still needs the axis/carry and source guard for descent.

### An all-height hostile consisting of actual positive sources

The inherited negative-17 mechanism in
[the return-family note](collatz_minus17_return_20261003.md) makes the
failure concrete at arbitrary finite length. Extend the formal score to
all positive valuation letters by `c(a)=3-2a`, retaining the weights +1,-1
at a=1,2. For the marked negative-cycle word

    W=(1,1,1,2,1,1,4),
    (P,Q,B)=(2187,2048,2363), a_W=-17,

the score is -1 but every nonempty prefix has expanding coefficient. For
each m>=1, define the actual positive odd integer

    n_m=2*2048^m-17.

It has exactly the valuation prefix `W^m` and endpoint

    U^(7m)(n_m)=2*2187^m-17>n_m.                     (6)

Every intermediate odd iterate after the source is also strictly above
the immutable source n_m. At the completed blocks the scores are
-1,-2,...,-m, although every positive-source prefix grows. These are
arbitrarily long finite examples. At m=1 this is `4079 ->4357` after seven
odd steps.

For exactness, let c_i be the ith point of the repeated negative cycle
and A_i its accumulated valuation. The formal ith state is

    c_i+2*3^i*2^(11m-A_i), 0<=i<=7m.

Each is an odd integer, since c_i is odd and the added term is even. This
proves that every prescribed division has exactly its stated valuation.
Every repeated-word prefix has coefficient greater than one, and its
positive carry then proves it exceeds n_m. The endpoint is (6).

This uses a previously proved shadow family; only its role as a curvature
hostile is added here. The source changes with m. There is no conclusion
that one positive source follows the negative cycle forever, nor a claim
that these growing prefixes by themselves provide a completed route home
or new coverage.

## 5. Join to the incoming boundary compiler

For a contracting coarse cylinder `n=r+Qk`, the incoming sharp theorem is

    a_w/Q=B/[Q(Q-P)] <= 1121/3328 < 1/2.

Its geometric reading is exact: the fixed-point threshold lies below Q/2,
while every lift with k>=1 lies above Q. Therefore all such positive lifts
descend; only the least positive member r can require a separate certificate.
Coarse cylinders may have extra final halvings, so their actual endpoint is
the odd part of the formal affine endpoint. Formal descent is sufficient;
the converse applies only to the exact cylinder (3).

Enlarging only the final valuation changes Q but leaves P and B fixed. The
axis moves from the negative side to the positive side exactly when the
coefficient crosses from growth to contraction. The finite interface test
uses every prefix-growing word of length at most six, with one, two, or
three repetitions: 165 rows, 164 descending boundaries and root 1. It checks
the same-axis repetition rule, the guard's axis-side change, the sharp bound,
and a positive lift against the inherited compiler. No new coverage delta
is asserted.

The user's root rays with H-exponents `5+6b,3+6b,7+6b` have exact odd
valuations `6+6b,4+6b,8+6b`. All these single-letter forward maps contract.
Thus the three residue rows are not spherical/Euclidean/hyperbolic cases
of the word drift. Their legitimate ternary-address and completed-route
structure remains in the [inverse-ray codec](inverse_ray_ternary_addresses_20261004.md).

The remaining global question is still how to choose a successful guarded
word or completed route for every source. Neither a metric name, a spectral
class, nor a local curvature sign supplies that selection theorem.

## Reproduction and audit scope

```text
python -X utf8 -B 04-computation/experiments/geometry_collatz_drift_carriers_20261004.py
python -X utf8 -B -O 04-computation/experiments/geometry_collatz_drift_carriers_20261004.py
```

[Script](../../04-computation/experiments/geometry_collatz_drift_carriers_20261004.py)
and [output](geometry_collatz_drift_carriers_20261004.out). The explicit universe:

- all 4,095 positive valuation words with total at most 12, comparing a
  chronological carry update with the independent sum, decoding the word,
  checking cyclic marking and J, and replaying three exact integer lifts;
- all 65,025 ordered pairs whose two totals are at most eight, checking the
  order defect, commutation and common primitive root independently;
- all 210 words on the `j=7,A=11` clock, retaining all their distinct axes;
- 144 rational-score controls `c_1=1,...,12`, `c_2=-12,...,-1`, finding each
  first strict wrong sign by exact powers; the latest first error is at
  length 53 for weights `(7,-5)`;
- the inherited negative-17 shadow sources for m=1,...,24, with every
  prescribed valuation and every source-relative growth comparison replayed;
- the 165 incoming compiler rows and the two source/axis hostiles above.

The normal and optimized runs use the same explicit checks. The finite
counts are controls for the proofs; they are not evidence of universal
Collatz convergence. No new geometric literature theorem is imported.
