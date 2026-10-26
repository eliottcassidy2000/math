# Multiperfect numbers: independent scoped audit

**PASS: elementary exact identities and independent finite controls.**
No odd multiperfect nonexistence claim, completeness classification, or new-priority claim is made.

## 1. Exact family theorem

For integer `r>=2`, let `M=2^r-1` and `D=binom(2^r,2)=2^(r-1)M`. Since `M` is odd, divisors of the two coprime factors pair uniquely, and

```
sigma(D) = (1+2+...+2^(r-1))*sigma(M) = M*sigma(M),
sigma(D)/D = sigma(M)/2^(r-1).
```

The denominator after the cancellation is solely a power of2. Therefore
`D | sigma(D)` iff `v2(sigma(M)) >= r-1`—both directions, for this exact family. Since `M>1`, its divisors1 andM already sum to`2^r`; the integer multiplier is at least2. It equals2 iff these are the only divisors, equivalently iffM is prime.

The report properly preserves the distinction between perfect and multiply perfect. `r=4` gives `D=120`, `sigma(D)=360`, multiplier3. The recorded `r=10` example gives `M=1023=3*11*31`, `sigma(M)=4*12*32=1536`, `D=523776` and `sigma(D)=1571328=3D`. Thus even the tested range has more than one composite-apex multiperfect case. The root added this second example to the prose without changing the computation.

The transport retains an arithmetic predicate on a specified scalar arc count and discards every tournament orientation. It is not a structural theorem about tournament paths, LRC or topology. Every`D` here is even, so no odd multiperfect witness is represented by this family.

## 2. Odd valuation ledger and exact failure of sufficiency

For odd`n=product p^a`, multiplicativity gives

```
v2(sigma(n)) = sum over odd a of (v2(p+1)+v2(a+1)-1).
```

For even`a`, the geometric sum has an odd number of odd terms, so valuation0. For odd`a`, write `a+1=2^s*u` with`u` odd and`s>=1`. Factor `p^(2^s*u)-1` into `p^u-1`, `p^u+1`, and`s-1` factors of the form `p^(2^j*u)+1`. The first factor has valuation`v2(p-1)`; the second has valuation`v2(p+1)`; each later factor has valuation1 because an odd square is1 modulo8. Dividing by`p-1` proves the displayed formula. This elementary argument independently establishes the ledger, rather than assuming a paper's full claims.

If `sigma(n)=k*n`, oddness of`n` forces this ledger to equal`v2(k)`. Each odd exponent contributes at least1. Consequently odd`k` forces every exponent even, so`n` is a square;`k=2` forces exactly one odd exponent with its prime and exponent both1 modulo4. The latter is Euler's necessary form, not an existence theorem. The general odd-exponent count is at most`v2(k)`.

The source citation is accurate: Chen–Luo section2, Theorem2.1(1), contains the prime-power identity, and Theorem2.2 develops the odd structure. The report does not import a current nonexistence result from this2011 preprint or misread restricted nonexistence examples as the global conjecture. [Chen–Luo primary paper, section2](https://arxiv.org/html/1102.4396v1)

Three decisive scope hostiles:

* `21`: `sigma(21)=32`, but`21` does not divide32. Large binary valuation does not supply the odd divisibility outside the arc family.
* `5`: correct Euler-form congruences and`v2(sigma(5))=1`, yet`sigma(5)=6`, not10. Even the exact binary budget for`k=2` is insufficient.
* `9`: a square and`sigma(9)=13` odd, yet not multiperfect. Squareness is not sufficient for odd-multiplier existence.

## 3. Independent controls

The new checker never runs or imports the producer's checker. It uses a divisor-incidence sieve through100,000, independently matching paired-divisor enumeration on1..10,000. It checks the valuation ledger and odd-sigma/square equivalence for all50,000 odd inputs through99,999, plus the Euler budget-one condition and odd-exponent-count bound.

For each of the23 exponents2..24, it independently factorizes and enumerates the full divisor sets of bothM andD. Every abundancy, valuation and family decision matches the producer's stored rows. Normal and`-O` outputs are byte-identical; checks use exceptions, not removable assertions.

The independent multiperfect census through100,000 is

```
(n,k) = (6,2), (28,2), (120,3), (496,2), (672,3),
        (8128,2), (30240,4), (32760,4).
```

There is no odd witness in this finite interval. That is a bounded exact observation only, not evidence of exhaustive odd nonexistence. Likewise, the family test is complete for each suppliedr, but testing2..24 is not a classification for allr.

## 4. Inheritance and search semantics

The inherited HYP-3546 snapshot byte-matches the pinned repository source. Part(A)'s Euclid–Euler arithmetic is valid. Its additional parity/Ky Fan, Ham-Sandwich, LRC and topology prose is retained only as historical mixed-status material, not used in the proof above. Scalar coincidences do not recover the erased tournament data.

For any specified positive integer, divisor enumeration terminates. Scanning odd`n>1` and halting on `n|sigma(n)` is therefore a correct semidecision description for existence of an odd multiperfect number; because`sigma(n)>n`, any found integer multiplier is at least2. It does **not** verify a particular compact machine table. Odd-perfect existence implies odd-multiperfect existence by`k=2`; the contrapositive nonexistence implication is sound. No converse, mutual reduction or finite-graph isomorphism was established.
