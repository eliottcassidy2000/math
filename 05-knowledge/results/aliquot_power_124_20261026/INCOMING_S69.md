# Late integration: twin-prime collisions and the binary abc map

**PROVED elementary identities below; CITED independently audited incoming
classifiers.** The [S69 packet](../session69_prime_ports_power_aliquot_20261026/README.md)
arrived at commit2f483a9201 during final integration. It independently recovers
the complete small aliquot fibres, operation models, and open-boundary guards.
Its introductory synthesis was read. The following short proofs make the
specific additional connections imported here independently reviewable.

## Three meanings of a later occurrence

Our delayed-row theorem asks whether the **minimum** doubling-row root of
n² exceeds n. It includes odd n, and first occurs at3. Restricting to even n
gives first6; every representation then uses a positive number of doublings.
S69 uses precisely this even/positive-doubling predicate on odd roots.
Its natural density among all integers is1/2, consistent with our relative
density1 among even integers. Including odd bases in our broader predicate
raises its density to1.

If one merely asks whether n² appears in **some** later doubling row, while
allowing even row roots, the first positive-doubling example is instead
`4²=2·8`, with8>4. It fails the minimum-root predicate, since16 first appears
in row1. Thus the apparent competing first examples3,4,6 correspond to three
explicit predicates, not a contradiction or an ordering of infinite rows.

## Twin primes explain the first double finite fibres

For odd twin primes p,p+2, the proper-divisor formulas give

    s(2p)=1+2+p=p+3=s((p+2)²).

This constructs two different preimages of the same target. At p3 it gives
the exact fibre `{6,25}` of target6; at p5 it gives `{10,49}` of target8.
At p11 it gives preimages22,169 of target14, without claiming these exhaust
that fibre. The complete inverse bounds in the aliquot note are still
required to certify exhaustion. These witnesses assume no infinitude of
twin primes.

The semiprime source itself is a fixed point exactly when `2p=p+3`, namely
p3. Our unique6 theorem is stronger than this observation: it excludes
**every** perfect target of a prime square outside6, rather than only
perfect targets in this twin-prime family. Complete fibres before6 show
that6 is also the first target with a finite fibre of size at least two;
F(1) was already infinite and is excluded from that statement.

## A third abc map, and an exact overlap at1+8=9

For a delayed even base `n=2^a u`, with a>=1, u odd and u>2^a, retain the
positive difference and form

    (A,B,C)=(u-2^a,2^a,u).

It is primitive because an odd number is coprime to a power of two, and
A+B=C. Conversely any primitive positive abc triple whose distinguished
second term is `B=2^a`, a>=1, has odd C>B and reconstructs the unique
delayed even base `n=BC`. This proves a bijection with that **tagged binary
subfamily** of abc triples, not with all primitive triples.

The base6 yields1+2=3. The base72=8·9 yields1+8=9. Independently, our
prime-power aliquot identity at `p=3,k=2` gives the very same triple:

    1+(3-1)s(9)=1+8=9.

This is an exact overlap of two maps, with distinct source integers72 and9.
The ordinary normalized aliquot map at9 instead gives9+4=13. A proposed
unification therefore needs the **map identity** as a sidecar in addition
to any removed gcd. Equality of output triples alone cannot reconstruct
which arithmetic construction produced them.

These identities prove neither abc nor a uniform bound on triple quality.
The incoming fixed-core estimates keep core-dependent constants, while our
prime-power radical estimate retains the abc conjecture as an antecedent.
The two bounds are not interchangeable. No additional global theorem is
imported from S69 without its hypotheses.

The late note was proof-checked by the integrating root against literal
divisor sums and the inverse formulas. The earlier independent integration
audit remains pinned to the pre-appendix files; this addition does not
silently change the files or claims that audit reviewed.
