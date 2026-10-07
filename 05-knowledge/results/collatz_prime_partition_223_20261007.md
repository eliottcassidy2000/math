# 223, 233, 322: prime-factor tubes and source-preserving certificate sections

**Status.** PROVED: the elementary tube invariant, simultaneous join-fiber
parameterization, smaller-child transport, and all scope obstructions below.
CONDITIONAL ON AN ACCEPTED PAPER PREMISE: the Poisson–Dirichlet transfer in
section 3. FINITE-EXACT: the stated script universe. The supplied papers are
accepted as premises for this investigation; their proofs are not audited.
Universal positive-integer Collatz coverage remains OPEN.

## 1. Inheritance and numeric roles

The closest proved mechanism is the ordered carry plus exact source cylinder
in [THM-4489 — Fermat sextics and 223 returns](../../01-canon/theorems/THM-4489-fermat-sextic-223-collatz-return-budget.md).
Its cost-nine expanding return is the canonical hostile to universal modular
return drift. The corrected near miss is that all five completed 233-return
tails in [modular return debt, section 4](modular_return_debt_20261004.md)
already had an earlier ordinary descent. Added marked returns were not added
generic stopping coverage. The underused sidecar is the whole source/word/
endpoint incidence relation, not its quotient cardinality.

The concept board is **prime-predecessor factors / initial one-run / ordered
carry / common endpoint / source section / paid dependency**. The anchor is
an exact supplied-source certificate interface, the niche is a factor-law
transfer, and the wildcard is the difference between a reverse injection and
a section. The relevant method is controlled forgetting with an explicit
restoration sidecar.

These are actual numeric roles, not references to matching theorem numbers:

* The unreduced odd-word carries 223 and 233 occur for words (2,3,1,1) and
  (5,2,1), with maps (81n+223)/128 and (27n+233)/256. Interchanging the first
  two entries of the second word gives carry 149, despite identical length
  and total cost. This is inherited in
  [the triplet/CRT comparison](triplet_crt_223_233_20261004.md).
* The same note recovers binary orders 37 and 29 at primes 223 and 233,
  their sextic/octic coset observers, and the marked tournament response
  233,291,123,223. Those numerical observables retain different predicates.
* The legacy [golden numerical record](golden_oscillation_s116l.out) contains
  F13=233 and L12=322. Here their exact recurrence check is
  L12=F11+F13=89+233. Also 322=2·7·23; 323=17·19 is not prime, so 322 is not
  an example of a prime predecessor in the paper's sampling universe.
* Every nonempty positive odd-Collatz word has odd unreduced carry: its first
  carry summand is odd and every later summand is even. Therefore 322 cannot
  play that same carry role. Its direct dynamical role is an even source:

      322 -> 161 -> 121 -> 91 -> 137 -> 103 -> 155 -> 233.

  The initial arrow is one halving; the six odd-map valuations are
  (2,2,1,2,1,1). This literal connection needs neither primality nor Fibonacci
  interpretation.

## 2. What is imported from the two papers

**Accepted prime-factor premise.** Theorem 1.1 of
[The Poisson–Dirichlet law for prime predecessors](https://github.com/openai/math/blob/main/preprints/The-Poisson-Dirichlet-Law-for-Prime-Predecessors-September-24-2026/paper.pdf)
states that, for uniformly sampled primes 3<=p<=X, the ranked logarithms of
the prime factors of p-1, with multiplicity and divided by log(p-1), converge
in every fixed finite joint distribution to PD(1). This is a prime average,
not a pointwise assertion or a stated law on an arbitrarily chosen thin
progression. The prior
[seven/twenty-one note, section 4](seven_twentyone_mersenne_openai_math_20261006.md)
already separates broad prime statistics from thin Collatz clocks.

**Accepted set-theoretic premise.** Theorem 1.1 of
[The Partition Principle does not imply Choice](https://github.com/openai/math/blob/main/preprints/The-Partition-Principle-does-not-imply-Choice-September-24-2026/partition-principle-without-choice.pdf)
is a relative-consistency statement for ZF+PP+ACWO+not-AC. Its introduction
defines PP by: a surjection f:X->Y implies some injection j:Y->X. It stresses
that this injection need not satisfy f∘j=Id. Our section construction below
is explicit and does not require this set-theoretic theorem; the distinction
identifies precisely which certificate coordinate a mere injection loses.

Both primary PDFs were read through their actual theorem statements. No
proof-validity conclusion about either paper is made here.

## 3. A factor law that survives adaptive stopping inside the one-run

For a positive integer m, let K6(m) be the integer obtained by deleting all
factors 2 and 3, including their multiplicities. Let p>=3 be prime, and set

    p-1=2^a t,    a=v2(p-1),    t odd,    n0=p-2.

Then the first a-1 odd-map valuations are exactly one, and

    n_j+1 = 3^j 2^(a-j) t,       0<=j<=a-1.             (1)

Proof: before the last index, v2(n_j+1)>=2, so n_j=3 mod4 and
v2(3n_j+1)=1. The identity follows by induction. At the last index
v2(n_j+1)=1, so the next valuation is at least two. In particular,

    K6(n_j+1)=K6(p-1)                                  (2)

throughout this complete initial tube. Therefore (2) remains exact at any
selected index J(p) in {0,...,a-1}, even if the selection depends on the
entire source or the observed prefix.

Write q_i(p) for the ranked prime factors of p-1 with multiplicity, padded
by ones. Write r_i(p) for those of K6(p-1), also padded by ones. Every
factor above 3 stays in its original leading position; after these factors
are exhausted, r_i=1 and q_i<=3. Thus, for every i,

    |log q_i(p)-log r_i(p)| <= log3.                    (3)

Keeping the original denominator log(p-1), (3) tends uniformly to zero as
p grows. Uniform continuity of a bounded continuous test on [0,1]^k and
the accepted Theorem 1.1 prove:

> The first k ranked logarithmic factor coordinates of K6(n_J(p)+1),
> divided by log(p-1), have the same PD(1) limit over all sampled primes,
> for every fixed k and every such adaptive stopping rule J.

This conditional transfer neither normalizes by the new state's size nor
conditions on an arbitrary source cylinder. The retained original scale
and the stopping restriction are essential.

**Reset and drift hostiles.** For p=233, the tube is 231,347,521 and its
kernel is 29. The next odd state is 391, with K6(392)=49: the invariant does
not pass the reset. Moreover p=59 and p=233 both have K6(p-1)=29, but their
sources 57 and 231 satisfy U(57)=43<57 and U(231)=347>231. The factor
quotient has erased binary fuel. Neither its exact value nor its limiting
distribution supplies a paid source guard.

## 4. Exact simultaneous endpoint fibers

For each i in a finite nonempty bank, suppose a supplied positive odd n_i
has the actual valuation word w_i of length r_i and cost A_i, ending at the
same positive odd H. Let its exact carrier be

    F_i(n)=(3^r_i n+B_i)/2^A_i,
    R=max_i r_i.

The complete family of simultaneous positive odd endpoint lifts is

    H(t)=H+2·3^R t,
    n_i(t)=n_i+2^(A_i+1)3^(R-r_i)t,                    (4)

where t is an integer and every displayed source and endpoint is positive.
In particular t>=0 always works.

**Necessity.** Given a common odd endpoint z, source integrality requires
2^A_i z-B_i=0 mod3^r_i, hence z=H mod3^r_i. The strongest condition is
z=H mod3^R. Oddness combines it with z=H mod2 to give exactly (4).

**Sufficiency.** These equations give integer positive sources. They also
give actual words, not merely formal rational maps. Working backward from
the odd endpoint, each intermediate formal inverse is integral: its
numerator is divisible by the remaining power of 3 by reduction of the
full source congruence. Each intermediate is odd because its numerator is
2^a times an odd integer minus one, with a>=1, divided by 3. Positivity
follows forward from the positive starting source. Thus every division has
exactly the prescribed valuation. Empty words cause no exception.

For a selected row, the program solves t from the supplied source and
rejects a nonintegral, negative, or mismatched parameter. It then returns
the actual other sources. Its forward-ray API explicitly uses t>=0; the
mathematical formula also describes any permitted negative parameters.

## 5. The concrete 223/233/322 packet and a paid splice

The inherited common future is H=425:

    223 --(1,1,1,1,3)--> 425,
    233 --(2,1,1,1,2,3,1,1,2,1)--> 425.

Their odd lengths and Terras costs are respectively (5,7) and (10,15).
The pair fiber (4) is exactly

    child m=223+62208t,
    source n=233+65536t,
    common endpoint=425+118098t,       t>=0.             (5)

It provides the uniformly smaller dependency

    m=(243n+469)/256<n,

since n>=233>469/13. Given an authenticated first-ROOT word for m, strip its
five-letter join prefix and prepend n's ten-letter prefix. The result is
n's actual first-ROOT certificate. The child prefix cannot hit ROOT early:
it ends at a value at least 425. This is a conditional certificate transport
with a supplied child proof, not an assumption that every member of (5) is
already rooted.

For the full packet, adjoin the odd source 161 whose word to H is the
six-letter word from section 1 followed by the 233 word. Formula (4) gives

    161+33554432t,
    223+45349632t,
    233+47775744t,
    common endpoint=425+86093442t.                     (6)

The even source corresponding to the first row is 322+67108864t. It carries
one initial halving as a separate coordinate. The three original source
Terras clocks to 425 are 7,15,25; their odd-edge counts are 5,10,16. These
are asynchronous receipts, not an identification of clocks.

No new generic first-descent coverage is claimed: every source in (5)
already descends on its first odd step, whose valuation is two. The useful
output is the exact reusable smaller-child splice and the complete fiber,
not a purported increase in the existing descent atlas.

## 6. Sections, cardinal injections, and uniqueness

Let X be a bank of verified records (source, ROOT word), and let pi:X->Y
forget the word. A valid source-specific selector requires

    pi(s(y))=y.                                        (7)

With the two records for 223 and 233, swapping them is an injection Y->X
and both records are valid. It violates (7) at both inputs. The checker
rejects using the 233 proof as the supplied proof for child 223. This is a
finite mathematical instance of the injection/section distinction, with
actual source guards retained.

Formula (4) is a stronger positive repair: it constructs a section of the
common-endpoint relation and a computable inverse parameter, so (7) can be
checked directly. No form of Choice is needed for this explicit arithmetic
construction. More generally, finite certificate codes are naturally
well-ordered; searching for their first valid witness is a partial
algorithm on the set where a witness exists. A uniqueness convention does
not establish existence or termination on every supplied source.

## 7. Tests and unresolved interface

Run the matching script normally and optimized:

    python -B 04-computation/experiments/collatz_prime_partition_223_20261007.py
    python -B -O 04-computation/experiments/collatz_prime_partition_223_20261007.py

The explicit universe is: all 1,364 words of lengths 1..5 with exponents
1..4 for carry parity; all 60,000 positive odd endpoints below 120000 for
the pair-fiber iff; t=0..63 for the pair/triple fibers and supplied-child
ROOT splices; all 668 odd primes through 5000 for 1,313 initial-tube states
and eight factor coordinates; and ten type/source hostiles. The literal
ROOT search only supplies independent finite controls; it is not called by
the conditional transport API. Saved normal/-O outputs agree.

The live missing interface is source-specific: an independently proved
factor or partition statement would need to yield a member of the correct
guarded certificate fiber, together with its height, ordered word, and
grounded terminal. Neither accepted paper supplies that assertion. The
new tube transfer and explicit fiber construction identify what does pass
between the subjects without silently providing the missing witness.
