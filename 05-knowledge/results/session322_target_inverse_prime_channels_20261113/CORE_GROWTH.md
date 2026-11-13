# A quantitative large-core obligation for absolute selection

This concerns the two fixed [selective quotient charts](DYADIC_INVERSE.md),
but allows **arbitrary** tournaments in both modules. It is independent
of the [one-apex word obstruction](ONE_APEX_CAPACITY.md). In particular
it does not assume that an SCC is an apex over a chain.

## Size costs and what the output actually sees

For either chart, with h=H(U), the coefficient A satisfies
A>=9h+2Q_U(3) since Q_U(2)>=2h. Since Q_V(2)>=2H(V), the resulting count
obeys
```
N = A Q_V(2) − h H(V) >= 2 Q_U(3) Q_V(2)
  >= 2 · 3^|U| · 2^|V|.
```
The second inequality follows by assigning each vertex to one of the
optional slots and using at least one Hamilton path on each nonempty
assigned subset. Thus N<=X bounds both module orders by O(log X).

All the chart can see of U is the product of the triples
(H,Q2,Q3) of its strong components. For V it sees only the product of
the pairs (H,Q2). This follows from ordinal multiplicativity, and
the component order does not affect these products. The output
graph may be strong nevertheless; the restriction is on the two
inserted modules, not on that output.

## Core-size thresholds that still give only X^(o(1)) outputs

Let b_U,b_V bound the largest SCC allowed in the two modules.
A j-vertex tournament has
```
1 <= H <= j!,    Q_t <= j! binom(j+t−1,t−1).
```
For the upper bound, concatenate the ordered paths to one permutation
and record all slot lengths, including zeros; this coding is injective.
The number of possible triples of order <=b is consequently at most
exp((3+o(1)) b log b), and the number of pairs at most
exp((2+o(1)) b log b). This does not require knowing which tuples
are physically realizable. There are O(log X) possible multiplicities
per such tuple; ignoring their total size constraint can only
increase the estimate.

Fix any eps>0 smaller than 1/3. If, as X grows, the permitted sizes satisfy
```
b_U <= (1/3−eps) log log X / log log log X,
b_V <= (1/2−eps) log log X / log log log X,
```
then the number of available tuple types on each side is at most
(log X)^(1−c+o(1)) for some c>0. Raising O(log X) to the sum of
these type counts yields exp(o(log X))=X^(o(1)) possible target
integers <=X. A fixed finite bank is an even smaller special case.

Thus even allowing the SCC size cutoff to grow slowly, no selector
restricted to those cutoffs can cover all primes (or all 7p targets).
To cover a positive fraction of either class up to X, most witnessed
targets must escape at least one of these two limits. This is a
necessary growth requirement, not a claim that these thresholds
are sharp or that larger cores suffice. It also does not constrain
other physical quotient charts with additional ordered contacts.
