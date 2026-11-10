# Primorial centers, the offset-one bridge, and the two Brocard problems

Session243, 10 November 2026. **PROVED elementary statements; FINITE-EXACT controls.**
The independent arithmetic audit accepts the stated scopes; see
`../arithmetic_audit/REPORT.md` and THM-5110.
Polignac, prime-square Brocard, factorial-square Brocard, and infinitude of
primorial primes remain OPEN. No literature-priority claim is made.

## 1. Recovered objects and provenance

The closest inherited mechanism is
`05-knowledge/results/session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md`,
Sections2–6: exact Fortunate/lesser-Fortunate definitions, the shared two-sign
wheel, the prime-offset cutoff, and the distinction from Euclid–Mullin.
The canonical hostile is `30030+19=151*199`; surviving the wheel is not
primality. The less-used sidecar here is the **deliberately omitted offset1**.
The new synthesis below proves exactly how to restore it to recover an
ordinary consecutive-prime gap. Source blobs are pinned in `source_pins.json`.

There are two unrelated named Brocard questions in the earlier work:

* **Factorial-square/Brown problem:** `n!+1=m^2`, conjecturally only
  `(n,m)=(4,5),(5,11),(7,71)`. See
  `05-knowledge/results/collatz_doubling_tower_brocard_20260927.md`, Section3.
  The historical note records the owner's request to "consider how
  {4, 5, 7} are the only factorials one less than a square, and those squares
  are {5, 11, 71}". These last three are the **square roots**, and completeness
  is conjectural, not a proved consequence of that wording.
* **Prime-square conjecture:** at least four primes strictly between the
  squares of consecutive primes `p<p'`, with `p>=3`. See
  `05-knowledge/results/procgen_brackets_20260924_pairings_transitions.md`,
  Section4. The case `p=2,p'=3` has only5 and7 and is excluded.

No literal Polignac occurrence was recovered in the tracked `.md`, `.py`, or
`.out` search at the pinned baseline. We therefore use the standard meaning
explicitly: for every positive even g, infinitely many **consecutive** primes
have difference g. Merely finding two primes separated by g is weaker.

Write `P_n=product_{i=1}^n p_i` for the ordinary primorial. The plus Euclid
number is `P_n+1`; minus primorial numbers `P_n-1` are distinguished explicitly.
Neither must be prime: `P_6+1=30031=59*509`. Euclid–Mullin instead multiplies
its own prior prime outputs and takes the least prime factor of product+1.
Fortunate offsets are different again.

Concept board: primorial center; skipped offset1; consecutive versus arbitrary
prime pair; prime-power sign packet; physical count versus arithmetic prime.

## 2. Exact restoration of the omitted offset

For n>=2 define

    F_n = min{d>1: P_n+d is prime},
    L_n = min{d>1: P_n-d is prime}.

Both exist: the positive direction has unbounded primes, while `P_n>=6`
has at least the prime2 below `P_n-1`. Let the actual nearest-prime distances be

    f_n = min{d>=1: P_n+d is prime},
    l_n = min{d>=1: P_n-d is prime}.

Then the exact identities are

    f_n = 1 if P_n+1 is prime, and F_n otherwise;
    l_n = 1 if P_n-1 is prime, and L_n otherwise.             (1)

Since P_n is composite, the nearest primes on its two sides are consecutive.
Their gap is therefore exactly `f_n+l_n`. This is an unconditional statement
about this specified gap, not a prime-gap asymptotic or an infinitude theorem.

**Corollary (equal offsets and the missing coordinate).** If `F_n=L_n=d`,
the primes `P_n-d` and `P_n+d` are consecutive if and only if **both**
`P_n-1` and `P_n+1` are composite. The full set of intervening primes is
exactly the prime members of `{P_n-1,P_n+1}`.

**Proof.** Every interior integer has offset0 or an integer offset between1
and d-1. Offset0 is P_n, composite. Minimality of F_n and L_n excludes all
offsets2,...,d-1 on both sides. The two offset1 cases are exactly those not
tested by their definitions. Equation(1) follows from the same exhaustive
partition. QED.

This is stronger than merely warning that the two endpoints might not be
consecutive: it identifies every missing intermediate-prime certificate.

### Exact hostile and positive controls

At `P_6=30030`, inherited certificates give

    F_6=L_6=17;
    P_6-17=30013 prime, P_6+17=30047 prime;
    P_6-1=30029 prime, P_6+1=30031=59*509.

Thus the pair at offsets±17 has difference34 but is **not** consecutive.
Its only intervening prime is30029. The actual crossing gap at P_6 is
`30047-30029=18`. At `P_3=30`, F=L=7, but both29 and31 intervene between23
and37; the crossing gap is2. At `P_4=210`, the two offset1 numbers209 and211
have respectively composite/prime status, and the crossing gap is12
(199 to211), not F+L=24.

The converse is nonvacuous: at `P_8=9699690`, F=L=23 while both P_8±1 are
composite, so the two offset23 endpoints are consecutive, with gap46.
The verifier checks these facts by trial division and scans the full intervals.

## 3. Fixed non-twin gaps cannot persist at primorial centers

**Proposition.** Fix an integer d>1. Let r be any prime divisor of d.
For every n with `r|P_n` and `P_n-d>r`, both `P_n-d` and `P_n+d` are composite.
Consequently there are only finitely many n for which both are prime.

**Proof.** Both integers are divisible by r. The strict lower inequality
makes r a proper divisor of both. The prime r is eventually among the first
n primes, and P_n tends to infinity, so both conditions eventually hold.
QED.

The exceptional equal-to-divisor endpoint is essential: `P_2-3=3` is prime
despite a common divisor. It has not been erased by a coprimality slogan.
For even d, parity already excludes every sufficiently large pair, since
P_n is even. For odd d>1, the first prime divisor of d supplies the eventual
obstruction.

**Consequence for Polignac.** The fixed-center template `P_n-d,P_n+d`
cannot prove infinitely many prime pairs of any fixed gap `2d>2`, even before
the additional consecutive-prime requirement is imposed. This is a no-go for
this construction, not a disproof of Polignac. The center of a general prime
pair need not be a primorial.

The same argument applies to a fixed multiplier: for any fixed positive
integer c and d>1, the pair `cP_n±d` is eventually composite. Allowing d to
grow avoids the theorem but changes the fixed-gap quantifier.

### The twin boundary is exactly offset1

At d=1 no prime divides d; every prime divisor of P_n is avoided by both
P_n±1. If both endpoints are prime, they are automatically consecutive
because the only interior integer is the composite P_n. Thus simultaneous
primality gives a genuine gap2, with no omitted-interior issue. It is not
known here to occur infinitely often; the finite wheel alone cannot prove
either endpoint prime. Examples are5/7 at P_2,29/31 at P_3, and2309/2311 at
P_5. The sixth plus endpoint30031 gives an immediate composite control.

The numbers also have to be typed before entering the tournament spectrum:
`P_2+1=7` is prime but is not a Hamiltonian-path count (inherited THM-1370).
Primality or Euclid's coprimality mechanism does not manufacture a tournament.

## 4. Brown numbers retain whole prime-power packets

For n>=4, let `e_p=v_p(n!)`, and put `M=n!/4`. This is a positive integer.
Then

    n!+1=m^2 for a positive integer m

holds if and only if the complete prime-power factors of M can be partitioned
into two products x<y with `y-x=1`. The corresponding root is `m=x+y`.
Every prime power must be assigned whole to one side; splitting a prime's
valuation between the sides is not permitted.

**Proof.** A solution has m odd, and

    x=(m-1)/2, y=(m+1)/2,
    xy=M, y-x=1, gcd(x,y)=1.

Coprimality sends each full prime-power factor of M to one side. Conversely,
if xy=M and y=x+1, then `(x+y)^2-4xy=(y-x)^2=1`, so `(x+y)^2=n!+1`.
QED.

Thus the odd prime p-packet is `p^e_p`, while the 2-packet is
`2^(e_2-2)`; the latter exponent is positive for n>=4. Equivalently every
odd prime power dividing n! has one definite sign in

    m=+1 or -1 mod p^e_p.

This is the exact sign partition behind the inherited divisor involution.
It is much stronger than retaining the radical alone, but it is not a new
solution of the factorial-square problem. The consecutive-product condition
is a global equation that residue survival alone does not impose.

The three known solutions give

    n=4: M=6=2*3,              (x,y)=(2,3),   m=5;
    n=5: M=30=2*3*5,           (x,y)=(5,6),   m=11;
    n=7: M=1260=4*9*5*7,       (x,y)=(35,36), m=71.

For any such solution, `gcd(m,n!)=1`, because a common prime would divide
both m^2 and n! and hence their difference1. Thus every prime factor of m is
greater than n. This says that m is n-rough, **not** that m is prime.

The verifier enumerates every subset of the actual prime-power packet for
each4<=n<=40, rather than searching m below an arbitrary bound. It agrees
with exact integer-square testing and finds exactly the three displayed
solutions in that finite n-range. Both methods have finite scope.

## 5. Prime-square Brocard and the prime-gap coordinate

For consecutive odd primes p<p' write `p'-p=2d`. Then `(p^2,p'^2)` is
exactly the union of d odd-square brackets

    ((p+2j)^2,(p+2j+2)^2), 0<=j<d,

plus the shared square endpoints, which are composite. Its length is
`4d(p+d)`. Therefore the number of constituent brackets is half the original
prime gap. This is an exact interval decomposition, inherited in substance
from the September24 bracket note; it gives no prime lower bound by itself.

If Legendre's conjecture is assumed, each integer-square interval contributes
at least one prime. There are p'-p=2d such intervals. Hence Legendre implies
the Brocard four-prime assertion for every **non-twin** consecutive odd-prime
pair (d>=2). Twin pairs retain only a two-prime guarantee under that argument.
Oppermann's stronger four-per-odd-bracket guarantee, also cited in the prior
note, suffices for all these Brocard intervals.

Polignac for a fixed gap2d would provide infinitely many such unions with
exactly d brackets, but says nothing itself about how many primes lie inside
their squared intervals. These are conditional relations with explicitly
different predicates. No implication from either Brocard problem to
Chan–Pak is asserted.

## 6. Genuine connection to a physical count construction

The reusable obstruction is the preservation of the exact parameter and
witness predicate, not the coincidence of small prime values. The proved
four-count family in the companion construction lane is

    6Q(T)-H(T), 6Q(T)+H(T), 6Q(T)+3H(T), 6Q(T)+5H(T).

Its literal graph proof and independent physical audit are retained in
`../constructions/PROOF.md` and `../construction_audit/AUDIT.md`. It is an arithmetic
progression with step `2H(T)`, giving a concrete map from one graph to four
graphs with prescribed **count differences**. Those output integers need not
be prime or consecutive primes. The step depends on the actual input H.

The comparison is now exact:

* Primorial expansion with d fixed eventually destroys endpoint primality.
* The physical four-count construction preserves a finite count difference
  by a graph operation; it does not expand a prime sieve or promise primes.
* Fixing H would keep the physical count gap fixed, but H alone does not fix
  Q or supply infinitely many distinct outputs. The same scalar with a
  different response is meaningful data, not automatically a missing target.
* Conversely, an arithmetic prime pair supplies no graph witness; the
  Euclid endpoint7 is the decisive control.

Connection contract: source is a primorial and two exact nearest-prime
supports; target is its actual crossing gap; map is(1); preserved are both
endpoints and consecutiveness; discarded by F,L alone is offset1 primality;
repair sidecar is two primality bits; cheapest hostile is30030. For Brown
numbers the source is the full prime-power partition, target the root, map
`m=x+y`, predicate `y-x=1`; discarding valuations or this global equation
breaks the map. These are explicit arithmetic operations, not a proposed
unproved number-to-tournament equivalence.

## 7. Reproduction and limits

    python3 -I -B verify.py > CHECK.json
    python3 -I -B -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

The script uses trial division, complete interval scans and integer square
roots, with explicit failure exceptions. It checks primorial indices2..8,
the known30030 hostile, fixed offsets2..100 at a verified obstructing center,
all prime-power partitions for4<=n<=40, and actual prime-square interval
decompositions for the 29 consecutive odd-prime pairs from (3,5) through
(113,127). General statements follow from
the proofs above, not extrapolation from these controls. No conjecture,
global coverage theorem or literature novelty is claimed.
