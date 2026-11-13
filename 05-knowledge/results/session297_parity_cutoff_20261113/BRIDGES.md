# Parity carriers: prime pairs, gap coordinates, and curves

## Prime wheels: an exact equality of shadows

Let M be odd, and require x,y units modulo M. The local fibers of x+y=t
and x+2y=t are bijective by multiplying the second coordinate by two.
For p^a dividing M, their common cardinality is
p^(a-1)(p-1) if p divides t and p^(a-1)(p-2) otherwise: choose y avoiding
0 and the unique second forbidden residue modulo p (these coincide
when p divides t). Multiply these counts by CRT. Modulo 2^a (a>=1),
with x,y odd, each target of the appropriate parity has 2^(a-1) pairs;
the other parity has none. Thus after keeping the odd-component target
and swapping the parity channel, every finite wheel has the same
Goldbach/Lemoine marginal counts. The joint map (E,O)=(x+y,x+2y) is
invertible over any ring: y=O-E and x=2E-O. Inherited HYP-2218/2219
keeps the joint prime-pair coordinates. Doubling a unit modulo odd M
is a permutation, but doubling an odd prime in Z is not prime. The
wheel does not lift the bijection to prime pairs.

Dirichlet and CRT imply every parity-compatible class modulo any
modulus contains arbitrarily large actual sums of two odd primes, and
likewise values p+2q with both primes odd. This does not prove binary
Goldbach or Lemoine. Actual odd binary sums n>=5 are much sparser:
n has such a representation iff n-2 is prime. For even n, reversal of
an ordered pair (p,q) gives the standard unordered count
(R_2(n)+1_{n/2 is prime})/2; fixed points are exactly the diagonal
solutions. This is a reflection-orbit calculation, not a primality test.
Odd ternary Goldbach also has three-odd-prime options and is not a
translation of even binary Goldbach.

## Where the two sevens live

Session210 and THM-4940 use q=(H-1)/2, so missing H=7,21 gives exactly
q=3,10. Session211's square-sum counts at orders 23,25 are 3,10;
its no-transport result already rules out the proposed global
H=2p+1 converter. The repo's Collatz mod-6 sign-specific probes give
a literal {4,11} cycle of the least-exponent greedy inverse for 3x-1:
4 -> (8*4+1)/3 =11; 11 -> (11+1)/3 =4, taking the least nonnegative
exponent that yields an integer not divisible by three. For 4, the
smaller integer candidate at exponent 1 is 3 and is disallowed.
The {3,10} and {4,11} pairs also occur among the CRT runner pairs of
the n=14 construction (lrc-n14-metagraph-crt-proof-s524).
Both have distance seven, so both collapse modulo seven. As
reflection pairs their sums are 13 and 15: under x->x+1 the reflection
x->13-x conjugates to x->15-x on Z. Thus there IS a shift of the
reflection axes on the common integer line, but not a shift identifying
the two residue wheels Z/13 and Z/15 or transporting either physical
construction. The companion pair is not a second proved tournament
gap spectrum.

## A curve needs a sign or a contraction, not only a smooth picture

Session127's volume theory gives V(n)=pi^(n/2)/Gamma(n/2+1) and
V(n+2)/V(n)=2pi/(n+2). The ratio is strictly decreasing and crosses
one once per parity subsequence: peaks are 6 on the even subsequence,
5 on the odd, and the integer global peak is 5. This is a *proof*
of a finite cutoff in growth because the recurrence controls ratios
at every later n. It is not evidence that a histogram peak of missing
tournament values supplies such a recurrence. A more literal
tournament microcosm is THM-4900's marked gate resolvent: its
golden Möbius observer converges and eventually forgets the core at
that scale, while the full tournament still remembers it. THM-4527's
Syracuse label microcosm instead has the exact replication
F(4n-1)=F(n). An exact replication and a convergent quotient cannot
be exchanged without proving an appropriate lift.

Session31's entire MRB deformation gives an audited stress test:
the CM-adjacent approximation's small residual changes sign under a
roughly one-part-in-a-million kernel perturbation, and the independent
audit isolates a nearby unique zero. Near equality is therefore not a
transport principle. The complementary arithmetic fact is sharper:
if a discrete count is rigorously approximated by an integer predictor
to error strictly less than 1, it equals the predictor. Here one needs
the discrete lattice and a certified uniform tail bound, not the
shape of the analytic curve alone. Compare the mod-m lift condition in
[CUTOFF](CUTOFF.md): the lattice spacing becomes m.

Sources: `microcosm_gate_doubling_resolvent_20261106.md`;
`session127_gsma_norm_collapse_20261031/volume/THEORY.md`;
`session31_mirror_mrb_20261022/ENTIRE_MRB_DEFORMATION.md` and its
`analytic-audit/AUDIT.md`; THM-4527; HYP-2218/2219.

## A moving scale gives a rigorously located transition

For radius R>0 let W_n(R)=R^n V(n) including dimension 0. The unique
crossover radius between n and n+1 is

    R_n = Gamma((n+3)/2)/(sqrt(pi) Gamma((n+2)/2)).

These radii strictly increase by strict log-convexity of Gamma:
x -> log Gamma(x+1/2) - log Gamma(x) is strictly increasing on x>0.
Stirling gives divergence and the stronger expansion
2*pi*R_n^2 = n+3/2 + O(1/n). Thus, at radius R between R_(n-1)
and R_n, dimension n is the unique maximizer; only adjacent dimensions
tie at a crossover. Dimension 0 is maximal below R_0. In particular,
the maximizing dimension is 2*pi*R^2+O(1) as R tends to infinity.
Here phase transitions persist unboundedly, but each *individual*
dimension's tenure is finite. There is no fixed period in dimension
hidden in this construction; a monotone ordered family of thresholds
is enough to organize the infinite sequence.
