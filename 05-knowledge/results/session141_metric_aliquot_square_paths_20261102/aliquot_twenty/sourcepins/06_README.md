# The unique perfect target of a prime square, and certified small aliquot fibres

2026-10-26. **PROVED elementary deduction from the CITED Euclid–Euler theorem;
PROVED inherited inverse bounds; FINITE-EXACT recertification through60.**
Independent audit [PASS](../chains/AUDIT_ALIQUOT.md), including a separate
additive divisor-sieve verification. No novelty claim for the inverse bounds or the
small fibres, and no resolution of the odd-untouchable or odd-perfect-number
problems.

Let s(n)=sigma(n)−n be the sum of proper positive divisors, with domain
n≥1. Put F(m)={n≥1:s(n)=m}. Zero is an output; iterating beyond that terminal
output would require a separate convention for s(0).

## 1. Inheritance and the precise next result

[Session66, aliquot report §§1–4](../../session66_root27_aliquot_graph_interfaces_20261026/session66-aliquot/REPORT.md)
already proves complete inverse fibres through13, the sharp finite bound,
its stronger non-prime-square form, and a complete guarded inversion
algorithm. Its frozen computation supplies fibres through500. The present
checker independently recertifies its prefix through60. The
[S57 report](../../s57_untouchables_root27_20261026/README.md) explicitly
inherits those statements; its MISTAKE-672 correction preserves the small
single-source exceptions3 and7. We do not treat those statements as new.

The new local deduction recorded here answers a different question: why
the fixed point6 also receives25, and whether every perfect number receives
a prime square. The answer is that **6 is the only perfect number that does**.
This is a short consequence of classical structure, not a claim of
historical priority.

Board: exact inverse fibre; prime-square branch; perfect fixed point;
binary divisor chain; missing-target certificate. Anchor: classify perfect
targets of prime squares. Niche: certify the supplied small fibres.
Wildcard: keep the prime weights when passing from divisor chains to
arithmetic. Canonical hostile:28 is perfect but27 is composite, and
s(27²)≠28. Corrected near miss:9 belongs to F(4), not F(3). Least-used
sidecar: the residue modulo3 forced by the Euclid–Euler exponent.

## 2. Exactly one perfect target receives a prime square

**Theorem A.** For a prime p and positive integer N,

    s(p²)=N and s(N)=N   if and only if   (p,N)=(5,6).

**Proof.** The proper divisors of p² are1,p, so N=p+1. If p=2 then N=3,
which is not perfect. Otherwise p is odd, and N is even. By the classical
Euclid–Euler theorem,

    N=2^(r−1)(2^r−1),

where 2^r−1 is prime and consequently r is prime. If r=2 we obtain N=6,
p=5. For any remaining r, the prime r is odd and r≥3. Modulo3,

    2^(r−1)=1,   2^r−1=1,   hence N−1=0.

But N≥28, so p=N−1 is divisible by3 and greater than3, contradicting
primality. Conversely s(25)=1+5=6 and s(6)=1+2+3=6. QED.

The cited input is [Keith Conrad, *Perfect numbers and Mersenne primes*,
Theorem1.1, pages1–2](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/perfect-numbers-mersenne.pdf).
Its forward and converse proofs and the elementary fact that r must be
prime were read. The modulo3 deduction is given above in full.

The theorem quantifies over **all** perfect N, without assuming that odd
perfect numbers do not exist: the prime-square incoming edge itself
forces N to be even apart from the explicitly rejected N=3. It concerns
an incoming square of a **prime**, not every square or every prime power.
For example, replacing p by an arbitrary integer invalidates the identity
s(p²)=p+1; already s(16)=15≠5.

The exact local graph is therefore

    25 -> 6 -> 6.

The loop is perfectness, while the additional incoming arrow is a different
property. The theorem explains their coincidence at6 without identifying
the two properties generally.

## 3. Complete inverse fibres, rather than a source cutoff

We repeat the inherited bound to make the finite certificate self-contained.
If m>1 and s(n)=m, then n is composite. Let p be its least prime factor and
q=n/p. Since p≤q and1,q are distinct proper divisors,

    q≤m−1,   n=pq≤q²≤(m−1)².                         (1)

Equality holds exactly for n=p² with p=m−1 prime. If n is not a prime
square, then1,p,q are distinct proper divisors, giving

    4n=4pq≤(p+q)²≤(m−1)².                           (2)

These bounds prove finiteness of every F(m) for m>1 and certify that a
search has no missing sources. They do not merely support an empirical
choice of cutoff. For m=0, F(0)={1}, since every n>1 has divisor1. For
m=1, F(1) is the infinite set of primes: a composite also has a proper
prime divisor p and thus s(n)≥1+p>1.

| Target m | Complete fibre F(m) |
|---:|:---|
|0|{1}|
|1|all primes|
|2|empty|
|3|{4}|
|4|{9}|
|5|empty|
|6|{6,25}|
|7|{8}|
|8|{10,49}|
|9|{15}|
|10|{14}|
|11|{21}|
|12|{121}|
|13|{27,35}|

The full certified table for2≤m≤60 is in [results.normal.json](results.normal.json).
The only untouchables in that interval are2,5,52; the only **prime**
untouchables in that interval are2,5. These are explicitly finite claims.

Two empty fibres can also be settled without a table. For m=2, (1) gives
n≤1, inconsistent with n composite. For m=5, a prime square would require
the nonprime base4. Every other possibility satisfies n≤4 by (2), and
the only composite there is4, itself a prime square. Thus F(5) is empty.

To expose an inadequate search, F(60) contains3481=59², although3481>60.
Scanning sources only through the target bound60 would miss it. A target
bound must be transported by (1), or by a proved sharper enumeration.

## 4. The binary sum and the open prime-target boundary

The divisors of2^r form the chain1,2,4,...,2^r. For r≥1,

    s(2^r)=sum_(j=0)^(r−1) 2^j=2^r−1.

Thus1+2+4=7 is the aliquot image of8. It gives the terminal chain
8→7→1→0; it does not make7 or8 a perfect fixed point. More generally,
when2^r−1 is a Mersenne prime, the same chain has three arrows.

For **distinct** primes a,b, one has s(ab)=1+a+b. Consequently the
conjecture that every even integer at least8 is a sum of two distinct
primes would imply that every odd target at least9 is touchable. The small
odd targets1,3,7 already have sources2,4,8, leaving5 as the only odd
untouchable under that conjecture. Equal primes cannot be substituted:
s(3²)=4, although1+3+3=7. Failure to find a distinct Goldbach pair would
also not exclude every other source shape.

This is the precise classical sufficient implication stated in
[Pomerance–Yang, *On untouchable numbers and related problems*, §1](https://math.dartmouth.edu/~carlp/uupaper3.pdf)
and [Pomerance's2020 lecture, slide36](https://math.dartmouth.edu/~carlp/upintconf.pdf).
Both primary texts were read. The uniqueness of5 among odd untouchables
remains OPEN. The claim that2 and5 are the only prime untouchables is weaker
and is not established here either. No equivalence with Goldbach is asserted.

## 5. Verification, transport and remaining scope

Run from this directory:

```sh
python3 -B check_aliquot.py --output results.normal.json
python3 -B -O check_aliquot.py --output results.optimized.json
cmp results.normal.json results.optimized.json
```

The first exact path computes sigma multiplicatively by prime factors for
every1≤n≤3481=(60−1)². The second computes proper divisors directly and
also constructs each fibre separately using (2) plus the prime-square
exception. Both agree with the frozen inherited prefix in
[INHERITED_FIBRES.json](INHERITED_FIBRES.json). No inherited producer code
is imported. Normal and optimized receipts are byte-identical; the
manifest records their gate count and hashes.

Controls include the bound's equality case, the missed source3481,
composite-base and repeated-prime failures, exact binary chains, the first
eight Euclid–Euler numbers, and all prime bases at most1000 in Theorem A.
The finite perfect controls do not replace its quantified proof.

Source→target: a positive integer with its prime-labelled divisor diagram
maps to the sum of its proper divisor weights, then to its exact inverse
fibre or its fixed-point test. Preserved: s(n), the incoming source, and
primality guards. Lost by forgetting weights: even the value of s(p²),
as the equal-shaped chains for4 and9 give3 and4. Lost by retaining only
the output6: the distinct source diagrams of6 and25. Necessary sidecars:
prime weights, source/target direction, and the finite inverse bound.
The [inherited weighted-divisor/F1 bridge](../../session66_root27_aliquot_graph_interfaces_20261026/session66-root/DIVISOR_F1_BRIDGE.md)
develops this interface; no new F1 model is claimed here.

Remaining open questions include the complete set of prime untouchables,
odd-untouchable uniqueness and global aliquot dynamics. The prime-square
theorem does not settle the different problem of arbitrary prime-power
predecessors of perfect numbers.
