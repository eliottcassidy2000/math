# Full rational image and a six-label arithmetic carrier

Session97 second phase, 2026-10-29. New files only; the first-phase
`README.md` and `check.py` are unchanged in this phase. Companion:
`full_image.py` (stdlib exact arithmetic; imports `check.py`).

## Status and external input

**CITED INPUT:** `E(Q)=Z G ⊕ Z/6 T`, where

    E: y²=x³+109x²+224x,
    G=(-4,28), T=(56,728).

This is the generator assertion in [Bremner–MacLeod, Remark2.2 and Lemma2.1](https://ami.uni-eszterhazy.hu/uploads/papers/finalpdf/AMI_43_from29to41.pdf).
The prior repo audit correctly did not infer full saturation from a
bounded saturation run. Here we explicitly use the published generator
assertion, not that finite computation. Everything below is a deduction
from this input, standard isogeny/reduction theorems, and displayed exact
certificates. It is not an independently reexecuted descent proof of
that input.

The standard dual-isogeny theorem gives `phihat phi=[3]` and
`phi phihat=[3]`, over the same ground field; see
[Sutherland, Lecture7, Theorem7.8](https://math.mit.edu/classes/18.783/2017/LectureNotes7.pdf).
Good-reduction injection on prime-to-p torsion and the three-division
polynomial are standard elliptic-curve group-law facts. No BSD or Sha
finiteness assumption is used.

## 1. Rational torsion of the quotient curve, exactly

Write

    E': Y²=Z³+421Z²+43904Z,
    H=(-196,196), U=(0,0).

Both11 and17 are primes of good reduction. Direct enumeration gives
`#E'(F_11)=12`, `#E'(F_17)=18`, so rational torsion has order dividing6.
For completeness the11-primary part injects at17 and hence vanishes;
the17-primary part injects at11 and vanishes; every other primary part
injects at both. This avoids falsely assuming unrestricted injection
at the residue characteristic.

The three-division polynomial factors as

    psi_3(Z)=(3Z+784)(Z³+300Z²+9408Z−2458624).

The cubic factor is monic and has no root modulo19, so it has no rational
root (a rational root of a monic integral polynomial is integral).
At the other possible root `Z=−784/3`, the right hand side of E' is
`−15366400/27<0`, so there is no rational Y. Thus there is no nonzero
rational3-torsion. Since U has order2,

    E'(Q)_tors = {O,U}.                            (7)

The checker supplies every count, coefficient and mod19 test.

## 2. H is a generator—not an assumption

The first-phase isogeny is

    phi(x,y)=(x(x−56)²/(x−4)²,
              (x−56)(x²+44x+224)y/(x−4)³).

Exact addition and substitution give

    3H=(-225,210)=phi(G),       phi(T)=U.           (8)

The isogenous curves have equal rank, hence E' has rank1 from the CITED
input. Choose a free generator K of E'(Q), and write H=mK modulo torsion.
Then phi(G)=3mK modulo torsion. The dual-isogeny relation says the free
image index `|3m|` divides3: on the two rank-one free groups the integer
multipliers of phi and phihat multiply to3. Therefore `m=±1`.
In particular there is no hidden saturation index for H. We conclude

    E'(Q)=Z H ⊕ Z/2 U,
    phi(E(Q))=3Z H ⊕ Z/2 U,
    E'(Q)/phi(E(Q)) ≅ Z/3.                        (9)

This closes the first-phase restriction to `<H>` without assuming a
primitive H or searching for generators at bounded height.

## 3. One local test suffices on rational target points

H corresponds to q=0 in the symmetric root polynomial

    f_q(t)=t³−t²+qt−(6q−1)/7.

For q=0 its reduction modulo2 is `t³+t²+1`, rootless over F_2.
A root over Q_2 would be integral because f_q is monic over Z_2, so
it would reduce to a root. Consequently H is not in `phi(E(Q_2))`.
This implication uses the same explicit cubic/fruit lift as the first
phase; H is not a boundary point where that chart fails.

By(8), H has exact order3 also in `E'(Q_2)/phi(E(Q_2))`. Since U is in
the image, the classification(9) gives, for **every rational** Q on E',

    Q in phi(E(Q))  iff  Q in phi(E(Q_2)).         (10)

More explicitly Q=nH+eU lifts iff `n=0 mod3`. This is a local-global
statement for this particular rational isogeny image. It does not
classify all torsors, all Selmer classes, all elliptic curves, or Sha.
It also does not assert that merely reducing every arbitrary q modulo2
is a complete algorithm: denominators and repeated reductions require
retaining the actual Q_2 equation. The simple rootless polynomial is
used to certify the one nonzero generating coset.

## 4. Complete reconstruction with six named labels

For a rational target Q=nH+eU, `e in {0,1}`, the exact lift controller is

    reject if 3 does not divide n;
    if n=3m, return P_j=mG+(e+2j)T, j=0,1,2.      (11)

These are exactly all three rational phi-preimages. To reconstruct the
full ordered triple from the unordered q quotient, also retain the
orientation label `sigma in {+1,−1}` and use `sigma P_j`. The map
`P -> -P` swaps a and b; the kernel `{O,2T,4T}` cyclically permutes the
three coordinates. These six actions generate S3. Their order is not
a free informational bit: the cubic splitting condition must be met
before these six rational labels even exist.

The fruit map is

    [a:b:c]=[56−x+y:56−x−y:−56−12x].

Clear its exact denominators, divide the integer gcd, and choose positive
sum a+b+c. For m nonzero the point P_j has infinite order, so it cannot
be any of the six exceptional torsion points. Hence all pair sums are
nonzero and the original rational expression—not just its cleared
cubic—is valid. For m=0 every lift is torsion and every candidate has a
zero pair sum: the controller must reject that boundary if asked for a
fruit solution.

No two nonboundary rational coordinates coincide. Setting a=b=t in
sum-one coordinates gives `(t−1)(14t²−11t+1)=0`; t=1 is the excluded
point `[1:1:-1]`, and the quadratic has nonsquare discriminant65.
Thus the six labels give six distinct ordered primitive triples.

It follows from the CITED input that **all unordered rational fruit
triples** are indexed uniquely by

    (m,e),  m a positive integer, e in {0,1}.      (12)

Here quotienting by negation removes the sign of m, and quotienting by
kernel translation removes j. This is an exact enumeration grammar,
not a claim of polynomial-time computation of the discrete coordinate
m from an arbitrary large rational point.

## 5. The two unordered partners cannot both be positive

Translation by U on E' sends `Z -> 43904/Z`. Since
`q=−Z/196−1`, its quotient action is the rational involution

    J(q)=(1−7q)/(7(1+q)).                         (13)

The positive-real chamber requires q>1/6. In that chamber the numerator
of J(q) is negative and its denominator positive, so J(q)<0. Thus of
the two e-partners in(12), at most one is positive. This recovers the
prior pair-sum reciprocal positivity obstruction by a short quotient
calculation while retaining the exact source index and labels.
The two partners need not include any positive one at a given m.
Ternary updates of the integer source n remain exact, but neither a
residue nor an unlabelled quotient certifies eventual positivity.

## Reproduction and scope

Run `python3 full_image.py` and `python3 -O full_image.py` in this folder.
Besides the torsion and primitive-generator certificates, the script
checks40 unordered families (m=1..20, both e), all six ordered labels,
original nonzero-denominator equations, and the reciprocal guard.
The all-order proofs above are not extrapolations from that finite bank.
The next genuinely open computational obligation is an efficient
coordinate extraction/local test for arbitrary large target input;
repeating bounded multiples is not a substitute for that algorithm.

## 6. A constant-size ternary decision guard, without discrete logarithms

The coordinate-extraction obligation above is unnecessary for the
**membership decision**. Prime3 is good for E': its displayed
Weierstrass discriminant is2 modulo3. The six points of E'(F_3) are

    O,(0,0),(1,1),(1,2),(2,1),(2,2).

Reduction of H is `(2,1)`, of order6; reduction of U is `(0,0)=3Hbar`.
For `Q=nH+eU` its reduction is `(n+3e)Hbar`. Therefore(9) yields

    Q in phi(E(Q)) iff Qbar is O or (0,0).        (14)

This is a finite-field classifier for **every rational point of E'**,
not merely for a bounded bank. Unlike a search for n, it only needs exact
rational coordinates and reduction modulo3.

In symmetric coordinates `Z=−196(q+1)`, let `q=a/d` in lowest terms,
d>0. If3 divides d then Z is nonintegral at3 and the elliptic point
reduces to O. Otherwise it reduces to `(0,0)` precisely when
`q=−1 mod3`. We obtain the explicit arithmetic guard

    rational lift iff 3|d OR 3|(a+d).            (15)

The projective reduction of q in P1(F_3) has four values: infinity and−1
accept, 0 and1 reject. **Input guard:** q must already have a rational
W satisfying the quotient curve(2); arbitrary rational q with a passing
residue is not promised a split cubic. The boundary q=−1 passes but its
fruit denominators vanish; that separate boundary guard must remain.
For nonboundary rational quotient points, (15) is equivalent to rational
splitting of the fruit cubic.

This is a precise version of “a decision may need very little information
while reconstruction needs the source”: two accepting residue states
certify existence, but they do not reconstruct n, primitive integer
heights, the three cyclic owners, or the orientation. Recovering these
continues to require exact arithmetic data. The theorem does not imply
that denominators may be reduced before the curve/quotient input guard.

The finite checker tests(15) on both torsion translates of nH for
`−40<=n<=40`, excluding only O from the affine q chart, against the
all-order source criterion `n=0 mod3`. This finite bank audits the
formula; the proof is the rational group classification plus good
reduction, not extrapolation.
