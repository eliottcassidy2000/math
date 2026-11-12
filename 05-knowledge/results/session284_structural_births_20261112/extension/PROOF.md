# Signature-run grammar: batch growth, family births, and logarithmic fracture

**Status: PROVED below; independently reviewed.** Exact finite controls are
in `verify.py` and `receipt.json`. This is a complete *marked-extension*
grammar and a rational enumeration of each fixed lineage family. It is not a
finite list of all unmarked tournament cores, not a fast A000568 algorithm,
and not a Chan–Pak proof.

## 0. Inheritance and precise new scope

The canonical chain decomposition and the one-vertex alternating grammar are
inherited from [weighted-core transport](../../weighted_core_transport_20261020/README.md),
using [THM-4463](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).
The old packet already proves the exact labelled one-vertex rank kernel.
Session283 treats arc edits, not the batch extension considered here.

This note makes the all-at-once interface explicit: signature words, their
positive run lengths, the newborn tournament, and their symmetry action form
a lossless family language. It derives the rational marked-family series and
a closed all-batch survival law. The stochastic law is also a consequence of
iterating the inherited insertion kernel; no priority claim is made for its
probability mechanism. The connection to shattering adapts the repository's
[THM-1450](../../../../01-canon/theorems/THM-1450-linearity-trivialises-vc-and-the-two-vc-one-tournaments.md),
not a claim to have discovered VC dimension.

## 1. The complete batch normal form

A homogeneous pair `{u,v}` has identical orientation to every third vertex.
Its graph `J(T)` is a disjoint union of directed paths. Their vertex sets are
the canonical transitive modules. Hence every nonempty tournament has the
unique form

    T = Q[TT_(l_1), ..., TT_(l_b)],  l_i >= 1,

where `Q` has no homogeneous pair. Write `b(T)=b`. These are not necessarily
prime modules in the modular-decomposition sense.

Adjoin a marked set `B` of `k>=1` newborn vertices. Retain the tournament `R`
on `B`. For every old vertex `v`, write its signature

    sigma(v)=(1_[x->v] : x in B) in {0,1}^k.

Read signatures along the intrinsic order of each old chain. Compress each
maximal constant interval into a run. Let the resulting proper words be
`w_1,...,w_b`; consecutive letters in each word differ. Let `r_i=|w_i|` and
`r=sum r_i`. Retain one positive length `a_ij` for every run.

**Theorem 1 (lossless signature-run grammar).** The data

    (Q, R, (w_i), (a_ij))                                      (1)

classify every batch extension up to isomorphism preserving the old vertex
set. Their reconstruction is as follows. Form a skeleton `S` on `r+k`
vertices, one per old run and one per newborn. Runs from the same old chain
are ordered transitively; runs from different chains inherit `Q`; newborns
inherit `R`; cross arcs are prescribed by the signature letters. Replace
run vertices by transitive chains of lengths `a_ij`, and leave newborns
singletons. This recovers the exact child tournament.

Conversely every `Q` without homogeneous pairs, every `R`, every proper
signature-word tuple, and every positive run-length tuple realizes such a
marked extension. The parent weights are `l_i=sum_j a_ij`.

Two such descriptions are equivalent exactly when a core isomorphism on `Q`
and a newborn isomorphism on `R` send each ordered signature word to the
corresponding word and carry its run lengths along with it.

**Proof.** An old constant-signature run remains a transitive module: old
vertices outside it see it uniformly by the parent chain decomposition, and
new vertices see it uniformly by its signature. The reconstruction is thus
literal substitution. The parent and its canonical chains are intrinsically
recoverable from the child together with the old subset. Every isomorphism
preserving that subset induces the stated core permutation, preserves each
chain's transitive order, and permutes the newborn coordinates by an
isomorphism of `R`. It therefore preserves exactly the run words and their
lengths. This also proves the converse. QED.

The number of descriptions is not the number of *unmarked* children. Several
different old subsets in the same child can give different lineage types.
Canonical deduplication after forgetting the mark remains necessary.

### What merging can and cannot do

The skeleton's canonical chain quotient is independent of the positive run
lengths, by uniqueness of the inherited transitive-chain quotient. Therefore

    b(child) = b(S),                                          (2)

even when the original tournament has millions of vertices. Its exact core
weights are the sums of the run lengths and newborn singleton weights along
the corresponding skeleton chains.

A child chain cannot contain vertices from two distinct parent chains: its
old restriction would be a transitive module in the parent, contradicting
maximality there. But it **can** contain multiple signature runs from one old
chain, joined through newborn vertices. The smallest witness is an old arc
`a->b` with a new vertex inserted as `a->x->b`. The old signatures differ,
while the whole child remains transitive. Thus the tempting formula
`b(child)=r+number_of_new_only_chains` is false. Contracting the actual
skeleton, rather than counting runs alone, retains the missing bridges.

## 2. A rational series for each whole lineage family

Fix `Q,R,(w_i)` and vary all positive run lengths. Let `G` consist of all pairs
`(alpha,beta)` of automorphisms of `Q` and `R` preserving the signature-word
tuple. It acts on the `r` run positions. The action need not be faithful:
newborn-only symmetries may lie in its kernel.

**Theorem 2 (marked family count).** The ordinary generating function for
isomorphism classes of marked extensions in this fixed family, by total
vertex count, is

    F(z) = z^k / |G| * sum_(g in G)
                    product_(C cycle of g on runs) z^|C|/(1-z^|C|).  (3)

Its first nonzero coefficient occurs exactly at `r+k`. Every later
coefficient is computable from the same finite rational function, without
a further tournament census.

**Proof.** Theorem 1 identifies classes with `G`-orbits of positive integer
run vectors. A vector fixed by `g` is constant on each run cycle `C`. Its
contribution to total old size is `|C|` times a positive integer, giving the
geometric factor in (3). Burnside's lemma yields the expression. All run
lengths one realize order `r+k`, and fewer vertices are impossible. QED.

For example `Q=C3`, one newborn sink, and a single signature run in each
chain gives

    F(z) = z/3 * ((z/(1-z))^3 + 2 z^3/(1-z^3)).                (4)

Its coefficients at total orders4,5,6,7 are1,1,2,4. A one-chain signature
word `010` with one newborn instead has the rigid family series
`z^4/(1-z)^3`, with coefficients1,3,6,10 at those orders. These are exact
marked-family counts, not the complete A000568 values at those orders.

The grammar is concise in its *rules* but has infinitely many templates:
`Q`, `R`, and the proper words can all grow. This is the distinction between
an exhaustive recursive description and a finite closed classification.
Nothing forces an unbounded sequence of future layers to be predicted from
a fixed bounded core catalog alone.

### Root forgetting with its exact symmetry correction

For a fixed parent `T` of order `n`, let `D(T,k)` be all labelled choices of
`R` on `k` named newborns and all `nk` cross arcs; write `U_d` for the child.
Then the following is a generalization of the inherited one-vertex identity:

    A000568(n+k) = 1/binom(n+k,k)
       * sum_[T], |T|=n sum_(d in D(T,k))
           |Aut(U_d)| / (|Aut(T)| k!).                       (5)

Indeed `Aut(T) x S_k` acts on the batch data. The stabilizer of `d` is
`Aut(U_d; old subset)`. One orbit contributes the size of the corresponding
orbit of old subsets under `Aut(U_d)`. Summing those orbit sizes for an
unmarked child gives exactly `binom(n+k,k)`. The formula is exact, but needs
the actual child symmetries; simply dividing marked counts by that binomial
coefficient is generally wrong.

## 3. Hamiltonian semantics: do not discard the weights

The skeleton and its positive weights retain the entire child, not merely
its core rank. They can be evaluated through the inherited
[THM-3121 walk-content substitution kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
For skeleton adjacency matrix `A_S`, set

    W_S(x)=1^T X(I-A_S X)^(-1)1,  X=diag(x_i).

If `a_i` are all skeleton weights, including `a_i=1` for newborns, then

    H(S[TT_(a_i)]) = sum_(1<=c_i<=a_i) [x^c] W_S(x)
                                  product_i c_i! {a_i brace c_i}.  (6)

Here the braces are Stirling numbers of the second kind: a spanning cover
of a transitive tournament by `c` unordered directed paths is exactly a
partition of its vertices into `c` nonempty sets, each with its unique
increasing path. This specializes the already proved substitution law.
Endpoint-refined cover data must be retained for endpoint-sensitive
operations; scalar `H` or core rank alone cannot replace this information.
Formula (6) is an exact semantics, not an efficient all-target inverse.
The present grammar supplies no proof that all desired odd outputs occur.

## 4. What random growth does to an existing large structure

This section concerns *labelled uniform independent arcs*, not uniform
isomorphism classes or uniform choices of distinct child classes. Fix any
parent `T` on `n>=1` vertices with `b=b(T)`. Add `k>=1` new vertices with all
cross arcs and newborn arcs independently uniform. Put `N=n+k` and let `Y`
be the number of old homogeneous edges that survive in `J(child)`.

**Theorem 3 (exact batch survival and controlled rank law).**

    Y ~ Binomial(n-b, 2^(-k)).                               (7)

Writing

    epsilon=(nk+binom(k,2)) / 2^(N-2),

one has

    E b(child) = N - (n-b)/2^k - epsilon,                     (8)
    d_TV(Law(b(child)), Law(N-Y)) <= min(1,epsilon).          (9)

In particular

    (1-2^(-k))^(n-b)-epsilon <= P(b(child)=N)
                            <= (1-2^(-k))^(n-b).            (10)

**Proof.** No old non-homogeneous pair can become homogeneous: its old
separator remains. An old homogeneous edge survives exactly when the `k`
newborn vertices give its endpoints equal bits. For each newborn, the
adjacency bits at the roots of the parent chain forest and their differences
across its `n-b` edges are independent uniform bits: this is an invertible
linear change of variables over `F_2`. Consequently survival indicators
for all old homogeneous edges are independent, each with probability
`2^-k`, proving (7).

Every pair with at least one newborn vertex is homogeneous with probability
`2^(-(N-2))`. For an old/new pair the `n-1` old exterior tests constrain
independent new cross arcs; the other `k-1` tests compare fresh newborn and
cross arcs. For a new/new pair all `N-2` tests are independent equal-bit
comparisons. There are `nk+binom(k,2)` such pairs. Let `Z` count their
homogeneous edges. Since `J(child)` is a path forest,
`b(child)=N-Y-Z`, and `E Z=epsilon`, proving (8). Coupling with `N-Y` gives
(9) from `P(Z>0)<=E Z`. Pair-freeness requires `Y=0`; failure of sufficiency
has probability at most `P(Z>0)`, proving (10). QED.

For any sequence of parents with `m=n-b -> infinity`, choose integers `k`
so `m/2^k -> lambda in (0,infinity)`. Then `Y` tends to
`Poisson(lambda)`, `epsilon ->0`, and

    P(child is pair-free) -> exp(-lambda).                  (11)

This gives a precise logarithmic growth scale, approximately `log_2(n-b)`,
on which pre-existing transitive-chain compression is erased. It does not
say that pair-free children have been fully classified: almost all remaining
structure is then in the core itself. If `n-b=0`, there is no old homogeneous
edge to erase, and (8)--(10), not this Poisson scaling, are the appropriate
statements.

## 5. An explicit endless schedule of interface births

Define the **external interface width** of `T` as the largest `k` for which
there is a boundary set `B` of `k` vertices such that vertices outside `B`
realize every one of the `2^k` orientation signatures toward `B`.

**Theorem 4 (sharp external-interface onset).** The smallest tournament
order realizing width at least `k>=1` is exactly

    n_k = k + 2^k.                                          (12)

**Proof.** Different signatures require distinct witnesses, and every witness
is outside the `k` boundary vertices. Conversely take any tournament on `B`,
a transitive chain on `2^k` old witnesses, and assign each of the `2^k`
signatures once. Every cross arc is freely specified, so this is a tournament
of the required order. QED.

Thus successive exact onsets include

    k:      1  2   3   4   5   6    7    8     9
    n_k:    3  6  11  20  37  70  135  264   521.

This is a meaningful unbounded family of new relational patterns beyond
orders50 and200. It is **not** a claim that no other new patterns appear
between these thresholds, nor that every tournament at a threshold has the
property. The property is intrinsic and existential over `B`, so it is not
just a novelty of how one labels a construction history.

Ordinary out-neighborhood VC dimension allows the boundary vertices
themselves to serve as witnesses. Its sharp onset is instead `2^k`, already
proved in THM-1450. The extra `k` in (12) pays for a *disjoint* interface and
witness set. The distinction is essential.

A Gray-code ordering of the witnesses gives a particularly transparent
construction: successive old signatures differ in exactly one newborn
coordinate. Their total boundary complexity grows while each adjacent step
uses the simplest possible local change. This provides a concrete answer to
how higher-order relational novelty can coexist with a fixed elementary
growth rule. It does not show that finite low-order data uniquely determine
all possible future systems without the rule.

## 6. Reproduction and limits

Run `python3 verify.py` and `python3 -O verify.py` from this directory. The
script uses explicit failures rather than removable `assert`s. It verifies
literal skeleton expansion, canonical-rank preservation, exact batch means
and binomial survival laws on the finite universes listed in the receipt;
three marked-family orbit coefficient examples; independent vertex-DP
Hamiltonian counts versus (6); and Gray-code witnesses for `k=1..8`.

The all-order statements follow from the proofs, not extrapolation of those
controls. The rules classify every marked batch extension when the necessary
core/word data are retained. They do not yield a finite catalog of all cores,
a scalar recurrence depending only on old A000568 values, or Chan–Pak
surjectivity.
