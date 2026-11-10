# Session241: natural-number relations inside tournament constructions

**COMPLETE, 10 November 2026. Independently audited; all twelve normal/optimized replays PASS. Chan–Pak completeness remains OPEN.**

The answer is precise: tournament path counts can faithfully encode several
relations of the natural numbers. The relation and the representation matter.
An encoding is not an inverse that realizes any requested number as a path count.

## 1. Multiplication, divisibility, gcd and lcm have an exact carrier

[THM-5103](../../../01-canon/theorems/THM-5103-faithful-fermat-path-count-arithmetic-and-sharp-coprime-bank-cost.md)
uses inherited near-transitive gadgets with Fermat path counts
F_i=2^(2^i)+1. Assign F_i to the i-th ordinary prime and combine prime powers
by ordinal sums. Pairwise coprimality, not Fermat primality, proves that the
resulting map preserves multiplication, numerical gcd/lcm, divisibility and
coprimality for **every positive integer**.

The new synthesis includes a complete criterion: a scalar multiplicative bank
reflects divisibility exactly when each generator has a private prime divisor.
Injectivity alone fails, even for actual path-count weights 15, 51 and 85.
Within the inherited endpoint-reversal family, the first k Fermat gadgets
uniquely minimize total order among coprime banks: **2^k−1+2k**, for k>=1.
This is a family-scoped optimum, not a global tournament lower bound.

[Proof, construction and controls](constructions/PROOF.md).

## 2. A strong-tournament divisor lattice recovers the totient predicate

[THM-5104](../../../01-canon/theorems/THM-5104-strong-fan-divisor-lattice-chart-and-totient-predicate.md)
uses the inherited Fibonacci fans to construct a strong Q_n with
H(Q_n)=4^n−1. It needs 2n+2 vertices for n>=2, and the triangle for n=1.
The classical repunit gcd identity then gives

    gcd(H(Q_a),H(Q_b)) = H(Q_gcd(a,b)),
    a divides b iff H(Q_a) divides H(Q_b),
    phi(n) = #{1<=a<=n : gcd(H(Q_a),H(Q_n))=3}.

This is the entire divisor lattice in an explicit linear-size chart, with
no factorization of n needed for construction. Its internal lcm is not the
ordinary numerical lcm. The normalized scalar (4^n−1)/3 connects to the prior
Collatz trunk, but at n=3 equals the forbidden path count 21. Normalization
is not a lawful tournament operation.

The inherited constructibility theorem additionally implies phi(H(Q_n)) is
a power of two exactly for n in {1,2,4,8,16}. This is a different totient from
phi(n), and the finite list is unconditional because F_5 is composite.

[Proof, exact map and losses](repunit_chart/PROOF.md).

## 3. A rooted arithmetic family now has an exact inverse and coverage law

[THM-5105](../../../01-canon/theorems/THM-5105-two-burst-shifted-product-inverse-and-exact-dyadic-coverage.md)
uses a literal strong one-root tournament with spine word O O^a Z O^b Z,
a,b>=1. The initial O is essential. Its actual path count satisfies

    2H−3=(2^(a+2)−1)(2^(b+1)+1).

Binary runs uniquely recover both parameters, without factoring. For every
k>=4, the exact number of resulting counts at most 2^k is

    (k−4)(k−3)/2 + floor(k/2)−2.

Hence there are only (1/2)(log_2 X)^2+O(log X) targets up to X. This is both
an all-order constructor and an all-order sparse-image boundary, stronger
than a finite modular observation. Multiplicative order, Carmichael periods
and LTE retain precise prime-power information about the **shifted** count,
not primality of H itself.

[Proof and inverse](coverage/PROOF.md) · [Independent audit](coverage_audit/REVIEW.md).

## 4. Prior research was merged, not relabelled as new

The [provenance synthesis](provenance/README.md) and [16 source pins](provenance/source_pins.json)
recover the actual chains: Moon/near-transitive counts; prime circulant and
Wieferich automorphisms; strong-labelled Camion/Korselt congruences; repunit
multiplicity degrees; Sessions218–219's totient fibers and Fermat-factor
escape; and prior spectrum-monoid factorization/elasticity.

In particular, aggregate strong-labelled counts s(n), a fixed tournament's
H(T), totient inverse fibers, automorphism groups and physical recipe states
remain distinct. Prior smooth atoms and elasticity 5/2 are explicitly credited,
not claimed again. Equal-order Fermat factors 641 and 6700417 have different
richer continuations; a prime index may map to a composite count.

## 5. What advances, and what remains open

The concrete advances are two faithful arithmetic representations with exact
loss ledgers, a sharp gadget-bank cost, and a factorization-free physical
inverse with exact all-order image size. The proof target remains unchanged:
THM-4094's prime and seven-times-prime lanes still need arbitrary target
placement. These encodings preserve specified relations but do not achieve
Chan–Pak surjectivity. No new unrestricted finite coverage cutoff is claimed.

Useful next questions are structural rather than numerical coincidence:
can the shifted addition law h_(a+b)=h_a h_b+h_a+h_b be implemented by a
uniform physical gadget with controlled ports; can a richer burst grammar
retain a tractable inverse while increasing image size; and can a different
coprime bank beat the endpoint-family cost? These are questions, not claims.

## Reproduction and audit

From this directory run `python3 verify_all.py`. The portable standard-library
packet replays six programs, including independently written adjacency path
counts and inverses. Producer checks include 500,000 odd inverse targets,
40,000 chart pairs, 10,000 burst parameters and 72,000 prime-order predicates.
The independent audits use different algorithms and explicit finite universes.

[Embedding and chart audit](embedding_audit/AUDIT.md) · [Coverage audit](coverage_audit/REVIEW.md).
Finite controls supplement all-order proofs; no full tournament census is implied.

The inherited base is commit `8a6b203c9b33c87a0a75cfc1fd6038308ad92646`.
The dirty main checkout is preserved; this session uses an isolated worktree.

## Close-out

All six live concepts were revisited: valuation supports give the faithful map,
SCC order explains the sorting loss, rooted words give the inverse, prime
orders govern shifted divisibility, and aggregate Camion/totient fibers remain
distinct targets. The existing research-method cards “Separate local support
from bounded-height coverage” and “Keep observer, recurrence class, and finite
head” were used. No new method card is promoted from a single session.

The first two verified main checkpoints are
`3d42677c554d8c026229206c3497cd098cefdef8` and
`cc8fbf1f70e4a88522dc5443905af27d2f1fe18b`. This final checkpoint adds
maintained navigation, frozen checksums and close-out. No theorem completion
is inferred from finite testing. **Session241 has ended.**
