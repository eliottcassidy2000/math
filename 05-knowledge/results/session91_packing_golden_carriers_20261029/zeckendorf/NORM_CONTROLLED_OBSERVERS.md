# The golden norm controls exact modular observation memory

**Status:** PROVED elementary finite-state and Smith-form statements below;
FINITE-EXACT audits are separately scoped. The ordinary divisibility minimum is
classical. This session supplies a self-contained constructive proof, a general
linear-observer formulation, and exact interfaces to the repository's golden
quotients and joint guards. No external-priority or major-open-problem claim.

Source snapshot: `90e704a8f35614d2085778d1b518f03ff0d750ca`.

## 1. Inheritance and a precise change of question

The relevant source files are:

- `05-knowledge/results/duck_zeckendorf_20260925.md`: unique canonical words,
  two-coordinate charges, and the carry when repeated atoms are added.
- `05-knowledge/results/zeckendorf_guard_automaton_20261003.md`: the exact
  value/shift register, the previous-digit guard, and the upper bound2m² on
  live modular states. Its modulus3 example separates integers2 and8 after
  appending a zero. The supplied historical color stripe is not this charge.
- `01-canon/theorems/THM-4592-the-eleven-squares-fibonacci-torus-is-the-period-five-collatz-cycle-set.md`,
  with the Session74 all-period correction: golden cyclic words map to
  `Z[phi]/(phi^n−1)`, bijectively at odd period; at even period only the three
  zero words collapse. Those are inherited results, not proved anew here. Odd period is not the
  same hypothesis as odd group order: the period3 golden quotient already has
  four elements. Its even-period exception comes from the alternating
  primitive-period2 words, not from tournament automorphism parity.
- Session90's `trijection/TREE_COMPRESSION.md` and
  `carriers/ITERATED_TRIANGLE_FAMILY.md`: exact carrier changes need their frame;
  a named root does not universally recover an anti-involution.

A primary literature check found that the exact ordinary Fibonacci divisibility
minimum is already known: Charlier, Rampersad, Rigo and Waxweiler,
[*State complexity of testing divisibility*](https://orbi.uliege.be/handle/2268/64877),
EPTCS31 (2010),48–57, gives2m² states in the trim minimal automaton of greedy
representations of multiples of m. The complete automaton here has one extra
rejecting sink. We do not relabel this known theorem as a new discovery.

The useful extension asks: **which linear measurement of the two registers is
being preserved under every future legal digit?** Its exact answer is an
integer quadratic form, the same indefinite golden norm found by the packing
lane. This is an operational bridge, not a claim that a digit automaton is a
physical square packing.

## 2. Native digit dynamics

Use Fibonacci weights f_i=F_(i+2), beginning1,2,3,5,8,... . Read finite words
high to low. Empty words and leading zeros are allowed. A word is legal when it
has no adjacent11. For a legal prefix retain

    X = sum_i d_i f_i,
    Y = sum_i d_i f_(i+1),
    p = last digit, with p=0 for the empty prefix.

Appending digit d is legal unless p=d=1, and gives

    (X,Y,p) -> (Y+d, X+Y+2d, d).

Let A=[[0,1],[1,1]] and b=(1,2)^t. On the two registers the update is
v -> A v+d b. The matrix is invertible over every Z/m because det A=−1.
The previous-digit bit is separate from the registers: both values occur at
every register pair modulo m.

### Explicit reachability, without an irrational-density argument

Let h be the order of A modulo m. For m>=2, h>=3, since neither A nor A² is the
identity modulo m. Set L=h; for m=1 set L=3. Starting with p=0, the legal blocks

    U = 1 followed by L−1 zeros,
    V = 01 followed by L−2 zeros

end with p=0 and act on the register vector by

    U: v -> v + A^(-1)b = v+(1,1),
    V: v -> v + A^(-2)b = v+(0,1).

These two translations generate all of (Z/m)². Thus every register pair with
p=0 is reachable by nonnegative repetitions of U and V. Appending1 then reaches
every pair with p=1 because v -> A v+b is bijective. Appending11 reaches the
rejecting sink. This is a constructive finite reachability proof; the period
h controls one witness length, not the number of minimal states.

## 3. General observer theorem

Fix modulus m>=1, integers a,b, and target residue r. On legal words ask whether

    a X+b Y = r mod m.

Write ell=(a,b), and define the observation matrix

    O = [ ell ; ell A ] = [[a,b],[b,a+b]].

Let K=ker(O:(Z/m)² -> (Z/m)²) and J=|im O|. The complete deterministic automaton
for this language has the following exact minimum:

- if r is not in the image of ell, the language is empty and its minimum is1;
- otherwise its minimum is **2J+1**: two legality states for each observable
  register class, plus the rejecting sink.

The image condition is `gcd(a,b,m) divides r`. It is essential when the
observer is not surjective. In particular the identically zero observer has
three states for target0, and one for every impossible target.

**Proof: exact quotient.** Two register vectors differing by z have future
output differences ell A^j z after j appended digits; the digit contributions
cancel. Since A²=A+I, all these differences vanish iff ell z=ell A z=0, that
is, z belongs to K. Moreover A K=K, so both digit updates descend to V/K.
Appending digit1 changes the common affine offset but not this criterion.
Thus V/K together with p is a closed sufficient carrier.

**Proof: minimality for one Boolean target, not just all residue queries.**
Every quotient state is reachable by the preceding translation blocks. If two
register classes differ, at least one of ell A z or ell A² z is nonzero:
otherwise ell z=ell A² z−ell A z would vanish as well. Append respectively
one or two zeros. This is legal from either p and resets both legality bits to0.
Now use U,V translations to send the first output to r. Such a translation
exists because r and every current output lie in im ell. The other output
still differs, so exactly one word is accepted.

If the register classes agree but p differs, append1. One path is dead; the
other is live. Append0 and then translation blocks to reach the target. A live
state is distinguished from the dead state by the same reset-and-translate
argument. Thus every one of the2J+1 states is necessary. QED.

A concrete quotient coordinate is `(u,v)=(ell v,ell A v)` in im O. Its update is

    (u,v,p) -> (v+d(a+2b), u+v+d(2a+3b), d) mod m.

Not every pair (u,v) need occur. Using all m² pairs would conceal precisely the
compression measured by the theorem.

## 4. Smith form turns memory into a golden norm

Suppose (a,b) is not the zero integer pair. Set

    g = gcd(|a|,|b|),
    Delta = a²+ab−b².

Delta is nonzero: a nonzero rational ratio cannot solve the irrational golden
quadratic. The integer Smith invariants of O are

    g, |Delta|/g.

Indeed g is the gcd of all entries, and their product is |det O|=|Delta|.
An integer unimodular row or column operation is invertible modulo every m.
Multiplication by an integer d on Z/m has image size m/gcd(m,d). Therefore

    J = m² / [gcd(m,g) gcd(m,|Delta|/g)].

For a reachable residue target, the exact complete-DFA size is

    2m² / [gcd(m,g) gcd(m,|Delta|/g)] + 1.

For the ordinary value observer ell=(1,0), Delta=1, recovering the classical
2m²+1 complete-state count. For a primitive observer, only primes dividing its
norm Delta can reduce the two-register memory. This is a proved arithmetic
criterion; it is not a generic analogy between prime numbers and information.

Let phi²=phi+1 and phi'=1−phi. Then

    Delta = (a+b phi)(a+b phi').

With N=[[2,1],[1,−2]], we also have

    Delta = (1/2) ell N ell^t,
    A^t N A = −N,
    Delta(ell A) = −Delta(ell).

Consequently appending zeros changes the observer's golden norm only by sign,
and leaves the exact modular carrier size unchanged. The packing lane's
indefinite anti-isometry and this automaton's Smith obstruction are literally
the same quadratic form. No Euclidean metric preservation follows: N is split,
not positive definite. See the companion [metric analysis](../packing/GOLDEN_METRIC_AND_PALEY_TORI.md).

## 5. The eleven-cell quotient and its missing channel

For any digit word its algebraic golden reading z=sum_i d_i phi^i is

    z = (2Y−3X) + (2X−Y) phi.

Check the initial weights and recurrence, or note that the matrix changing
(X,Y) to these coefficients has determinant−1. It is an integral change of
coordinates, not a loss of information.

At phi=4 modulo11,

    z = 5X−2Y = 5(X+4Y) mod11.

Thus the eleven-cell golden channel is, up to an invertible scalar, precisely
the observer ell=(1,4). Its norm is−11. It needs only11 observable register
classes, or23 states when legality and the sink are included. By contrast,
ordinary integer value X modulo11 needs121 register classes, or243 complete
states. The smaller golden quotient is not a faithful replacement for every
integer congruence guard. Explicitly3 (word100) and15 (word100010) both have
golden residue5 and last digit0, but their integer residues modulo11 are3 and4.

There is an exact second channel. The golden polynomial splits modulo11 as
(t−4)(t−8). The two eigen-observers

    X+4Y,       X+8Y

recover (X,Y) modulo11 jointly, since their coefficient determinant4 is a unit.
Each channel alone has11 register classes; together they have121. This is an
actual CRT/splitting statement for Z[phi]/(11), not a claim that two independent
isomorphism classes automatically recover their shared alignment.

### Changing the modulus exposes the carry defect

Keeping ell=(1,4) while changing m from11 to121 does **not** retain a scalar
one-channel recurrence:4²−4−1=11 is nonzero modulo121. The norm formula gives
J=121²/11=1331, hence2663 complete states.

The correct Hensel lift is37:37²−37−1=1331=11³. With ell=(1,37) modulo121 the
norm is−1331, so J=121 and the complete automaton has243 states. Thus a genuine
lift of the golden eigenchannel restores scalar memory; merely copying the
old representative4 does not. This is a concrete exact-carry test when moving
between moduli, not a floating-point approximation to phi.

The period5 cyclic-word theorem and the11-cell packing quotient are inherited
from THM-4592/Session74. This observer calculation identifies which arithmetic
queries that quotient can carry. It does not infer positive-integer Collatz
convergence, packing optimality, or an LRC conclusion from the common ring.

### Split, inert, and ramified precision have different exact costs

At a prime-power modulus p^k, call an observer p-primitive if a,b are not both
divisible by p. The norm formula gives J=p^(2k−min(k,v_p(Delta))). Therefore:

- If t²−t−1 is irreducible modulo p, every p-primitive norm is a unit and
  **J=p^(2k)** for every such observer. Examples are p=2,3,7.
- If that polynomial has two distinct roots modulo p, either root lifts
  recursively to a root lambda modulo p^k: solve the next coefficient using
  the unit derivative2lambda−1. The observer(1,lambda) then attains
  **J=p^k**, the smallest possible value for a p-primitive observer. Example:p=11.
- For p=5, a primitive norm divisible by5 has valuation **exactly1**. Indeed
  Delta=(a+3b)² mod5; vanishing forces a=2b+5c with b a unit modulo5, and
  Delta=5(b²+5bc+5c²). Hence the minimum is **J=5^(2k−1)**. Beyond the first
  precision level, no primitive scalar eigenreader exists modulo5^k.

The lower bound J>=p^k in the split case follows already from the norm formula:
its denominator can remove at most p^k for a p-primitive observer. All three
costs count register classes; multiply by2 and add the sink for the complete
DFA. The splitting hypothesis, not a chosen decimal, determines the result.

This gives an exact interpretation of the old {2,3,11} inspiration:2 and3
require both registers at every power, whereas11 admits a one-channel reader
at every power if its eigenvalue is lifted correctly. Prime5 is a separate
ramified precision defect. It does not identify these primes with roles in
an unrelated tournament or packing theorem. In particular the11-state
complete DFA at modulus5 and the11-cell golden quotient at period5 are
*different counts of different objects*.

## 6. Joint moduli: an exact case of recovering joint data

For ordinary value guards at m and n, let L=lcm(m,n), d=gcd(m,n). A pair of
live modular states comes from one word iff:

1. their previous-digit bits are equal;
2. their X registers agree modulo d;
3. their Y registers agree modulo d.

CRT then recovers a unique register pair modulo L. Conversely every such pair
and either legality bit is reachable, by the translation construction. Hence
the exact joint live carrier has2L² states, not(2m²)(2n²). The redundancy factor
of the naive independent product is

    (2m²)(2n²)/(2L²) = 2d².

For coprime moduli this saves exactly the duplicated legality bit. For4 and6,
the product has2304 states while the joint carrier has288, a factor8. If the
Boolean targets are CRT-compatible, their conjunction is one residue guard
modulo L and has2L²+1 complete states; incompatible targets give the empty
language and one state.

More generally use finitely many linear observers modulo possibly different
moduli, with common modulus L. Intersect their observable kernels O_i in
(Z/L)². A reachable joint target has minimum complete carrier2[V:intersection K_i]+1,
by the same reset/translation proof. The target must lie in the joint image;
individual target reachability is insufficient.

This is a precise positive counterpart to the repository's joint-port hostiles.
The recovery works because explicit residue compatibility equations determine
one common source state. It does not transfer to arbitrary tournament orbit
marginals, whose missing transport can remain undetermined.

## 7. Boundaries of normalization and comparison with tree compression

Zeckendorf normalization chooses a representative of numerical value. The
register/observer theorem instead classifies **all future legal digit tests**.
Those are distinct tasks. Two canonical prefixes can have the same current
value residue but different shifted residues, so the latter must be retained.
Two golden residue channels can jointly recover the full register, while either
alone deliberately forgets information.

Session90's tree-edge compression is also invertible only with its retained
oriented tree. The exact common lesson is to specify the source frame and the
future operation, then compute the kernel of the actual observation map.
Here the kernel is O; in the tree map there is no kernel until the frame is
forgotten. We do not identify their operations or claim an isomorphism of the
underlying combinatorial objects.

The same distinction limits LRC/unit-distance transfers: a modular digit guard
is a congruence predicate, not a simultaneous metric separation condition.
The golden norm and its primes supply a concrete arithmetic interface. A new
geometric conclusion still needs a map preserving distances, inequalities,
endpoints, and the quantifiers of the target conjecture.

## 8. Reproduction

Run `python3 check_norm_observers.py` and `python3 -O check_norm_observers.py`.
The standard-library checker performs353,234 active checks in each mode,
with identical output. It uses exact integers and explicit runtime checks.
It compares direct Fibonacci values with the register reader; verifies legal
translation blocks over all register states for moduli1–40; minimizes every
observer/target language for all(a,b,r) modulo m<=8; independently builds the
observable-image automata for the11/121 and other hostile cases; and checks
all joint-modulus state compatibility for m,n<=12. It also exhausts primitive
norms for powers of2,3,5,7,11 in the stated output ranges and verifies eight
successive11-adic Hensel lifts. Generated JSON is a local
receipt, not a required input file.

The minimality proofs hold for every modulus; finite minimization is a control,
not the proof. The general norm/Smith formula is an elementary application of
observability and finite abelian groups. No claim of complete literature
priority, a polynomial-time tournament classification, or a major conjecture
closure is made.
