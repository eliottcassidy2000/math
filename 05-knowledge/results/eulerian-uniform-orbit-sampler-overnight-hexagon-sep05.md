# Direct uniform Eulerian orbit sampling from the fixed-pair lift

**PROVED ELEMENTARY / INDEPENDENTLY AUDITED.** The measure identity is
the independently audited Section5 of the
[native Boolean flow note](overnight_hexagon_sep05_boolean_flow.md).
This continuation implements that identity, proves its resource boundary,
and supplies exhaustive distribution controls. It claims no literature
priority and no polynomial preprocessing or Markov mixing theorem.

## Statement and algorithm

For each n>=1, there is an exact sampler of a labelled Eulerian graph on
n vertices whose **isomorphism class is uniform**. It requires
`exp(O(sqrt(n)))` bit operations and storage for preprocessing, followed
by expected polynomially many bit operations and random bits per sample.
It uses neither graph canonicalization nor automorphism computation.
The word Eulerian here means all degrees even, without connectedness.

For a partition lambda=(l_1,...,l_r) of n, set

```text
b(lambda)=sum_i floor(l_i/2)+sum_(i<j) gcd(l_i,l_j),
f(lambda)=b(lambda)-r+1_{some l_i odd},
z(lambda)=product_l l^(m_l) m_l!,
w(lambda)=(n!/z(lambda)) 2^f(lambda).
```

The inherited closed-rank proof shows that f is the dimension over F_2
of the Eulerian graphs fixed by a permutation of this cycle type. Store
all positive integer weights w in a cumulative table. For each sample:

1. Choose lambda with probability w(lambda)/sum_mu w(mu), using a
   uniform integer below the total and binary search in the table.
2. Form a canonical permutation g_lambda with these consecutive cycles.
   Build its orbits on unordered vertex pairs and their degree-parity
   matrix. Row-reduce over F_2, choose every free coordinate by an
   independent unbiased bit, and expand the resulting kernel vector
   into an edge set F fixed by g_lambda.
3. Choose a uniform vertex permutation sigma and output sigma(F).
   Its class is the required sample. The implementation also returns
   sigma g_lambda sigma^(-1) as a checkable fixed-permutation witness.

There is no rejection based on an unknown stabilizer. Random integer
generation uses ordinary rejection of an out-of-range bit string, with
expected fewer than two trials. This is why the random-bit complexity
is stated in expectation, not as a deterministic worst-case bound.

## Exact measure and resource proof

For fixed lambda, every conjugate g has exactly z(lambda) preimages
under sigma -> sigma g_lambda sigma^(-1). Relabelling maps the kernel
bijection onto Fix(g), so each pair (g,F) of this type has probability

```text
[w(lambda)/W] [z(lambda)/n!] [1/2^f(lambda)] = 1/W.
```

Thus all fixed pairs are uniform. Each isomorphism class O contributes
`|O| |Aut(F)|=n!` such pairs. Its output mass is `n!/W=1/q`, where q
is the number of classes. Individual labelled graphs instead have
probability `1/[q |O|]`. These are different distributions: at n=4
the three class sizes are1,3,4, so uniform labels would give masses
1/8,3/8,1/2, not1/3 each.

There are p(n) partitions. For t>0, the elementary partition product
gives

```text
p(n) exp(-nt) <= product_(j>=1)(1-exp(-jt))^(-1),
log(product)=sum_(k>=1) 1/[k(exp(kt)-1)]
            <= (1/t) sum_(k>=1)1/k^2 <=2/t.
```

Taking t=1/sqrt(n) proves p(n)<=exp(3sqrt(n)); no asymptotic partition
formula is a dependency. The table entries and total have O(n^2+n log n)
bits, since there are at most `n! 2^binom(n,2)` fixed pairs. Enumerating
partitions and computing each weight costs polynomial time in n per
partition, proving the stated preprocessing bound. Each sample searches
that table, constructs at most binom(n,2) edge orbits, eliminates an
n-by-b binary matrix, and permutes at most binom(n,2) edges. All of
these have polynomial bit complexity. Kernel caching is optional and
unneeded for the theorem.

## Controls and exact limits of the connection

The standalone checker exhausts all labelled Eulerian graphs and all
vertex permutations for n=1,...,5. It separately compares every
canonical fixed kernel with literal fixed graphs, checks the free-bit
bijection, proves equal class mass by summing exact type weights, and
checks each labelled stabilizer against the inverse-class-size law.
This is a complete small-universe distribution check, not a histogram.

It also compares the closed rank to freshly built literal matrices for
all271 partitions through n=12. The default n=30 run preprocesses5604
types and draws64 seeded samples, checking degree parity and the returned
fixed permutation. The exact number of classes obtained is

```text
623045225369818212764392744508736687716211191500581246494929916017252505791800234321992704.
```

Those large samples are implementation controls, not statistical evidence
for uniformity; the preceding fixed-pair calculation proves the law.

The source-to-target map now supplies exactly the endpoint measure used
by the native fractional-flow formula. It does not bound that flow's
congestion, and fast independent sampling does not imply that local
triangle moves mix rapidly. Nor does the output provide a canonical
encoding of its graph class. These are the destroyed coordinates and
remaining obligations, not implicit consequences of the sampler.

```bash
python3 -B 04-computation/eulerian_uniform_orbit_sampler_overnight_hexagon_sep05.py
python3 -B -O 04-computation/eulerian_uniform_orbit_sampler_overnight_hexagon_sep05.py
```

The default producer has545 optimization-live gates; normal and optimized
runs agree. A separate agent audited the full measure/resource proof and
literal implementation and independently replayed all545 gates. Raw LF
source SHA256 is
`10f909d32eb4f2b0cc1d636a568a951f14e8169e1627ffadd7b5351fad45dd51`;
[output](eulerian_uniform_orbit_sampler_overnight_hexagon_sep05.out) SHA256 is
`b05cc16c33025c00fecba832bee2a7cbba9b552bf190680000cb051476f5f857`.
The shared semantic digest is
`bd53537ed6ab901a5bffd55f2afcfda27a8c6e1d9e64c0b0e17a7a7604d358f2`.
