# Prime arrival clocks, 290/227, and the orientation bit in MRB

**PROVED elementary clock and CRT identities; FINITE-EXACT census of all
even 4<=m<=10000; INHERITED periodic-sum theorem for MRB. No probabilistic
law for deterministic first primes is asserted.**

## 1. The seven numerical rows are right; their interpretation needs a clock

For even m>2 let p_*(m) be the least prime congruent to ±1 modulo m.
Let tau(m) count primes coprime to m up to and including that prime, and
let omega(m) count distinct prime factors. Then

    p_*(m)>=m-1,
    tau(m)=pi(p_*(m))-omega(m),
    L(m)=pi(m-2)-omega(m),
    delta(m)=tau(m)-L(m)=pi(p_*(m))-pi(m-2).             (1)

All prime divisors of m are below m-1, which proves the formulas.
The first L(m) eligible primes have exactly zero chance of a hit in this
deterministic experiment. The special modulus 2 is excluded: its two
nominal target residues coincide.

| m | phi(m)/2 | tau | p_* | zero-hazard prefix L | remaining prime rank delta |
|---:|---:|---:|---:|---:|---:|
|120|16|49|239|27|22|
|144|24|81|431|32|49|
|208|48|173|1039|44|129|
|262|65|54|263|53|1|
|278|69|57|277|56|1|
|290|56|227|1451|58|169|
|292|72|60|293|59|1|

The quotient U(m)/{±1} has phi(m)/2 elements. Independent uniform draws
from this finite group, waiting for the identity, have geometric mean
phi(m)/2. That is a specified probability model; Dirichlet's theorem or
equidistribution does not identify the increasing prime sequence with it.
Even adding L to that geometric mean would require a new stochastic
assumption about the remaining sequence. Equation (1) is exact instead.

The neighbor-prime examples have delta=1. They need not have a low raw
ratio: m=210 has p_*=211, tau=43, phi/2=24, ratio 43/24>1.79 despite
the immediate post-boundary hit. For m=290, L=58 already exceeds the
suggested iid expectation 56.

The exhaustive even-modulus census through 10000 gives a largest ratio
tau/(phi/2) of 6905/996, about 6.933, at m=9970, p_*=149551, tau=13810.
The first new strict ratio record after m=290 is m=1688, with p_*=15193,
tau=1772 and ratio 443/105. These are finite-universe records only.

Chebyshev bias compares cumulative prime-count races, typically with
logarithmic averaging and additional hypotheses in precise distribution
theorems. It does not by itself provide a first-hit waiting law. Both
1 and -1 are quadratic residues modulo 290, with four square roots each;
so even their quadratic character should be checked rather than guessed
from the fact that 1 is a square. The linked [2016 study](https://arxiv.org/html/1608.08647v2)
investigates Legendre-symbol walks, including Gaussian primes; it proves
no waiting theorem for 290. The linked Reddit post contains no 290 table.
For the distinction between races and sampling measures, see the primary
[Lichtman–Martin–Pomerance paper](https://arxiv.org/abs/1809.03033).

## 2. A deterministic infinite family behind the delay

For m=290 the candidates before 1451 are

    289=17²,       291=3*97,
    579=3*193,     581=7*83,
    869=11*79,     871=13*67,
    1159=19*61,    1161=27*43,
    1449=9*7*23.

Take M=2*3*7*11*13*17*19=1,939,938. The same proper small divisors
certify all these affine candidates whenever

    m=290+M t, t>=0.

Consequently p_*(m)>=5m+1 on this entire progression. Equality is not
asserted: 5m+1 may itself be composite. This is a finite covering of
candidate slots which persists under CRT, not a fitted probability.

**General CRT delay theorem.** For every K>=1 there is an infinite
arithmetic progression of even m for which every km±1, 1<=k<=K, is
composite, hence p_*(m)>Km. Choose distinct odd primes q_(k,s)>K and
impose m=-s/k modulo q_(k,s), together with m=0 modulo 2. CRT solves
these congruences. Discard the finite beginning where a candidate could
equal its divisor. Every later progression member works. Dirichlet's
theorem ensures p_*(m) still exists for each fixed modulus. ∎

This proves unbounded p_*(m)/m. It does **not** prove that tau/(phi/2)
is unbounded, since the CRT modulus itself grows and prime ranks use
another scale. Keeping these clocks distinct is part of the theorem.

## 3. The less obvious 227 connection is also a covering certificate

Prior [Stern-prime work](../special_primes_curling_genus_zero_langlands_20261015.md)
records that 227 cannot be written q+2a² with q prime, a>=1. The complete
finite obstruction is

    227-2a² = 225,219,209,195,177,155,129,99,65,27
    for a=1,...,10.

Every entry is covered by one of 3,5,11. Thus N=227+165t preserves this
obstruction for the first ten a values. It does not preserve being a
Stern prime: 557=227+2*165 is prime, and

    557=269+2*12².

The first eleven candidate differences for 557 are composite, so the new
successful branch is precisely a=12. This is the exact boundary of the
analogy: both constructions retain a fixed finite list of polynomial
divisibility obstructions; the admissible Stern range grows with N,
whereas the claimed prime-delay range was explicitly fixed at K.

The same integer 227 is also the retained-prime sum after deleting 11 in
[duck_primes_20260925](../duck_primes_20260925.md). That earlier work
distinguishes an original prime-index cutoff from a retained-prime count.
It supplies the correct *operation* to reuse here, not an equality of the
Stern predicate and a prime-waiting rank.

## 4. A 290/56 coincidence with a useful hostile test

The [THM-1430 tiling/metagraph dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md)
has 56 tournament classes and 290 class-metagraph edges at order six.
Meanwhile phi(290)/2=56. The finite group is explicitly

    U(290)/{±1} = (C4 x C28)/<(2,14)> ~= C2 x C28.

The Smith invariants are 2 and 28: the gcd of the presentation entries
is 2 and the gcd of its maximal minors is 56. Its element-order census
is 1,3,4,6,18,24 elements of orders 1,2,4,7,14,28 respectively.

This gives a candidate matching of 56 objects but rules out the simplest
structural identification: a Cayley graph on 56 vertices is regular,
whereas a graph with 290 edges has average degree 580/56=145/14,
which is not an integer. Thus the class metagraph cannot be that regular
Cayley graph. A useful replacement must retain its nonuniform incidence
or the H/Aut presentation weights; a bare 56-element group erases them.

Other recovered 290 occurrences are a discriminant valuation in
[THM-2469](../../../01-canon/theorems/THM-2469-degree-twenty-two-CD-plane-coprime-weight-ramification-closure.md)
and the slope (290,431) in
[THM-3645](../../../01-canon/theorems/THM-3645-berggren-two-cube-local-prime-support-and-mod24-cones.md).
The latter has an actual CRT/local-support mechanism worth comparing.
Neither occurrence makes the numerical prime-waiting rank invariant.
Gao's five-dimensional example also has determinant -290. Multiplying
one target coordinate by -1/290 normalizes that determinant to 1,
without changing fibres or injectivity. Thus its particular determinant
value is a coordinate choice, unlike the arithmetic modulus in (1).

## 5. MRB: merging the two orientations destroys cancellation

Let f(n)=n^(1/n)-1. The repository's
[periodic-filter theorem](../session52_saturation_residues_metagraph_20261024/mrb-supersingular/RESEARCH.md)
says that a periodic coefficient sequence a_n of mean mu satisfies

    sum_(n<=N) a_n f(n) = (mu/2)(log N)² + O(1).         (2)

For zero mean the sum converges. Indeed f(n)=log(n)/n+
O(log²(n)/n²); bounded partial sums of a_n-mu and summation by parts
give convergence of the zero-mean term, while the constant term supplies
one half log² N. This recalls the proof rather than claiming novelty.

For even m, every unit is odd. Therefore the alternating MRB sum filtered
to units has coefficients of mean -phi(m)/m. For m=290,

    sum_(n<=N, gcd(n,290)=1) (-1)^n f(n)
       = -(28/145)(log N)²+O(1).                        (3)

Filtering instead to the two target residues gives

    sum_(n<=N, n=±1 mod 290) (-1)^n f(n)
       = -(1/290)(log N)²+O(1).                         (4)

Both diverge. A density multiplier of the original MRB constant cannot
represent either filtered series. Retaining the orientation gives the
convergent replacement

    sum_(n>=1) [1_(n=1 mod290)-1_(n=-1 mod290)] f(n).    (5)

Character orthogonality makes the information loss precise. The sum of
the two indicators uses only characters chi with chi(-1)=1; their
difference uses only chi(-1)=-1. Passing to U(290)/{±1} erases exactly
the latter sector. There are 56 characters of each parity. The orientation
bit is useful because it removes the principal-character mean, not because
it supplies a universal prime bias or a new numerical formula for MRB.

This complements the earlier
[MRB contour audit](../session31_mirror_mrb_20261022/mrb/MRB_CONTOUR_AUDIT.md):
adding h_K(z)=e^(-z) product_(j=1)^K(z-j) leaves the first K samples,
the negative-axis jump and the full positive-real inverse-power/log
asymptotic expansion unchanged, but changes the alternating constant by

    (-1)^(K+1) K!/(e+1)^(K+1).

Thus even those combined observations do not determine the functional.
This is a counterexample in the stated admissible function class, not a
claim that the perturbed profile is still n^(1/n)-1. The existing CM
connection at A(8) survives; equality of the actual defect with q^4 was
already disproved, and MRB's arithmetic nature remains open.

## Reproduction

```sh
python3 -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/prime_wait.py
python3 -O -B 05-knowledge/results/jacobian_eth_prime_wait_20261102/prime_wait.py
```

Sieve counts are cross-checked by independent trial division on every
supplied row and the m=210 hostile. The infinite CRT claims have explicit
divisor certificates; five parameter values are additional checks only.
The finite census never exhausts its four-million prime sieve. The group
order histogram, Stern obstruction and all three periodic means are exact.
