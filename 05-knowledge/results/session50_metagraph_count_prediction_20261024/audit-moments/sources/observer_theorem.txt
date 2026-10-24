# Observer-weighted score moments and a nonreal-root boundary

Status: **PROVED** all-size identities and explicit counterexamples;
**FINITE-EXACT** independent labelled, full-record, and original-metagraph
controls. No literature-priority or new isomorphism-enumeration claim is made.

## 1. Inheritance, objects, and retained information

The frozen `sources/observer_completions_theorem.txt` supplies the exact
free-host completion bijection. Its fixed-host predecessor is retained in
`sources/fixed_host_global_observer_theorem.txt`. The present extension does
not replace the fixed-host transition test by a free-host argument.

The source objects are labelled tournament/observer records `(T,P)`, where
`P` is an actual directed Hamiltonian vertex order. Write `H(T)` for the
number of these orders, `d_v` for the outdegree of owner `v`, and

    X_v(T) = 2 d_v - n + 1,          Q(T) = sum_v X_v(T)^2.

Let `c3(T)` count three-element owner sets inducing a directed cycle. Fix
`N=binom(n,2)` and, for the unconditioned model, `f0=N-n+1`.

There are two distinct probability measures:

- `E_U` averages uniformly over the `2^N` labelled tournaments.
- `E_H` averages uniformly over all marked records `(T,P)`. Equivalently,
  its host probability is `H(T)/sum_S H(S)`.

The map `(T,P) -> (P, nonpath arc bits)` is bijective. It preserves the
entire named source and observer, not just their statistics. Projecting to
`T` loses the marked observer and has fibre size `H(T)`.

There is also an exact original-operation interpretation. Define `K_P(T)`
by retaining path arcs and reversing all other arcs. The pair of labelled
tournaments `(T,K_P(T))` determines `P`: their agreeing oriented arcs are
exactly its directed Hamiltonian path. Thus `(T,P) -> (T,K_P(T))` is injective
and counts compatible labelled source/partner pairs as well as marked records.
This recovery requires BOTH tournaments with their common owner alignment.
It need not survive independent passage to isomorphism classes; retaining only
the first tournament loses `H(T)` again. This is an agreement-path map, not a
claim about arbitrary graph complementation or arbitrary pairs of classes.

Live concept board: marked observer weighting; score/triangle identities;
free-arc conditional variance; fixed-path class fibres; scalar polynomial
root geometry. Their common mechanism is the explicit path-plus-free-bits map,
not native harmonic transport.

## 2. Exact unconditioned observer moments

**Theorem 1.** For every `n>=2`,

    E_H Q = n^2 - 3n + 4,
    E_H c3 = (binom(n,3) + n - 2)/4,
    E_U H = n!/2^(n-1),
    Cov_U(H,c3) = E_U H * (n-2)/4.                         (1)

For `n=1`, `H=1`, `Q=c3=0`; the first two displayed formulas are not used.
The `n=2` case gives `Q=2`, `c3=0`, and zero covariance.

**Proof.** Fix an observer order `P`. Its `n-1` directed path arcs are
deterministic. All `f0` other arcs are independent fair orientations. An arc
oriented `u->v` contributes `e_u-e_v` to the imbalance vector. The sum of the
path contributions telescopes to `e_first-e_last`, of squared norm `2`.
Each free arc contributes an independent mean-zero vector of squared norm
`2`. All cross expectations vanish. Therefore

    E[Q | P] = 2 + 2f0 = n^2 - 3n + 4.

This conditional mean is independent of `P`.

For every individual tournament, a transitive triple has exactly one owner
pointing to its other two owners, while a cyclic triple has none. Hence

    c3 = binom(n,3) - sum_v binom(d_v,2)
       = n(n^2-1)/24 - Q/8.                               (2)

The second equality uses `sum_v d_v=N` and
`Q=4 sum_v d_v^2 - n(n-1)^2`. Substitution of the first moment gives the
claimed cyclic-triangle mean.

There are `n! 2^f0` marked records, proving the formula for `E_U H`.
Under the uniform labelled-host measure, each triple is cyclic in two of
its eight orientations, so `E_U c3=binom(n,3)/4`. By the defining weighting,
`E_U[H c3]=(E_U H)(E_H c3)`. Subtraction proves the covariance. QED.

Positive covariance for `n>2` is an aggregate statement. It is not a
pointwise monotonicity theorem relating `H(T)` and `c3(T)`, not an extremal
characterization, and not an assertion that the two probability measures
are equal. Already at `n=3` their cyclic-triangle means are `1/2` and `1/4`.

## 3. Labelled checksums, fixed-observer class fibres, and the metagraph

The exact labelled audit sums for `n>=2` are

    sum_T H(T)       = n! 2^f0,
    sum_T H(T) Q(T)  = n! 2^f0 (n^2-3n+4),
    sum_T H(T) c3(T) = n! 2^f0 (binom(n,3)+n-2)/4.          (3)

Let `[T]` range over tournament isomorphism classes, and let `a(T)` be the
automorphism-group order. Orbit-stabilizer gives exactly `n!/a(T)` labelled
hosts in that class. Dividing (3) by `n!` yields

    sum_[T] H(T)/a(T)       = 2^f0,
    sum_[T] H(T) Q(T)/a(T)  = 2^f0 (n^2-3n+4),
    sum_[T] H(T) c3(T)/a(T) = 2^f0 (binom(n,3)+n-2)/4.      (4)

These are weighted checksums, not a formula counting unmarked isomorphism
classes and not a claim to enumerate A000568.

**Fixed-observer fibre theorem.** For ANY fixed labelled vertex order `P`,
the number of labelled realizations of class `[T]` containing `P` as a
directed Hamiltonian path is exactly `H(T)/a(T)`.

**Proof.** Choose one representative of `[T]`. Each of its `H(T)` Hamiltonian
orders `Q` has a unique owner bijection sending the ordered list `Q` to `P`.
Two such bijections yield the same labelled tournament exactly when they
differ by an automorphism. An automorphism fixing an ordered Hamiltonian
list fixes every owner, so the action on the lists is free. Each resulting
labelled realization is therefore counted `a(T)` times. QED.

In particular `a(T)` divides `H(T)`. Uniform fixed-path nonpath-bit words
have the same isomorphism-invariant law as uniform marked records. The same
law holds if one restricts to any nonempty subset of observer orders, or
gives those orders arbitrary nonnegative weights of positive total mass,
provided every order still receives all its free nonpath bit assignments
uniformly. Mixing identical class laws cannot change them.

The qualification on free bits is necessary. At `n=3`, fix `P=(0,1,2)` and
also fix the nonpath arc `0->2`. The only host is transitive, with `c3=0` and
`Q=8`, instead of the unrestricted fixed-path means `1/2` and `4`.

For a converse-merged class node, define its **word-fibre multiplicity**

    m([T, T^op]) = H(T)/a(T)       if T is self-converse,
                 2 H(T)/a(T)      otherwise.

Reversing Hamiltonian orders preserves `H`, and converse preserves `a`,
`Q`, and `c3`. Thus (4) is equivalently the all-size merged-node identity

    sum_nodes m Q  = 2^f0 (n^2-3n+4),
    sum_nodes m c3 = 2^f0 (binom(n,3)+n-2)/4.                (5)

The frozen ORIGINAL order-eight `merged_nodes.tsv` is copied, unchanged,
as `sources/merged_nodes.tsv.txt`. Its pinned source commit and path appear
in `sources/SOURCE_PINS.json`. The exact checker verifies all `3,528` rows:
it computes `H` independently by Hamiltonian-subset dynamic programming,
counts automorphisms by adjacency-preserving bijections, verifies the literal
converse-isomorphism pairing, and checks EACH stored `fibre` equals the
appropriate `H/a` or `2H/a`. It then obtains

    sum fibre           = 2,097,152 = 2^21,
    sum fibre Q         = 92,274,688,
    sum fibre c3        = 32,505,856.                       (6)

This is an actual original-metagraph enumeration constraint. The `m` here
is the word-fibre count, not an unexamined Laplacian diagonal or unit-resource
capacity. No stationary-walk claim is made. Any future port-degree comparison
must explicitly account for loop incidences twice; the ordinary Laplacian
diagonal drops loop contributions. Nor does this theorem identify individual
labelled observer lifts with additional disjoint native resources.

## 4. General coherent free-host completion moments

Use the frozen completion theorem's hypotheses: a nonempty named owner set,
coherent entire local source tournaments, and actual local Hamiltonian orders
whose entire local partners must be recovered. Let `e` be the number of
uncovered pairs. The union of local consecutive arcs must be a directed
linear forest; let its components be the `p` ordered path blocks.

An allowed block order `O` respects every declared GLOBAL clock pin. Such
pins concern cumulative block lengths and are not silently substituted for
local cell clocks. For this `O`, force the `p-1` joining endpoint arcs.
They are distinct uncovered pairs. All other `f=e-p+1` uncovered pairs are
independent fair choices. Write `b_O` for the divergence vector contributed
by ALL covered arcs together with these joining arcs, each counted once.

**Theorem 2.** If the allowed-order set `A` is nonempty, the uniform marked
completion measure satisfies

    E Q = 2f + (1/|A|) sum_(O in A) ||b_O||^2,
    E c3 = n(n^2-1)/24 - f/4
           - (1/(8|A|)) sum_(O in A) ||b_O||^2.             (7)

There are `|A| 2^f` marked completions. Each order has the same completion
count, so its marginal distribution is genuinely uniform in this marked
measure. Given `O`, the imbalance vector is `b_O` plus `f` independent
mean-zero arc vectors. The same second-moment expansion proves the first
formula; (2) proves the second. This includes `f=0`. QED.

If the forest fails, or global pins leave `A` empty, there are no completions.
The corresponding counting sums and polynomial are zero, but no uniform
probability measure or conditional mean is assigned. Forgetting the observer
would reweight source hosts by their number of accepted observers, so a
uniform unmarked-host assertion cannot be substituted into (7).

The forest gate is justified by the retained FULL local records: in any
completion every local path edge is globally consecutive; a coherent linear
forest forces each whole local cell path to be a contiguous block segment.
Consequently its local nonconsecutive flags are preserved too. Cross-block
pairs are uncovered, so arbitrary block permutations are legal before global
pins. The general claim concerns this free-host problem, not fixed hosts,
electrical guards, harmonic lifts, or routing between source problems.

## 5. A complete full-record nonreal-root counterexample

Use owners `0,...,7` in their natural reference order. Give four disjoint
two-owner cells the directed local source/observer paths

    B0=(0,3), B1=(6,1), B2=(2,5), B3=(4,7).

These coherent records form four blocks, so `p=4`, `e=28-4=24`, and `f=21`.
For a block order `O`, let `a(O)` count its three joining arcs oriented
forward relative to the owner reference order. Direct enumeration of all
`24` orders yields

    B(z) = sum_O z^a(O) = 4 + 10z + 10z^2.

For a hand-auditable certificate, block-index orders with `a=0` are
`0231,1023,2310,3102`; those with `a=1` are
`1032,1230,1302,1320,2031,2301,3021,3120,3201,3210`; those with `a=2` are
`0123,0132,0213,0312,0321,1203,2013,2103,2130,3012`. None has `a=3`.

The UNCOVERED-forward-arc completion polynomial is therefore

    G(z) = (1+z)^21 (4+10z+10z^2).                          (8)

The quadratic discriminant is `100-160=-60`. Its roots are
`(-5 +/- i sqrt(15))/10`, so (8) has nonreal roots despite fully coherent
local source/observer/partner records. It counts `24*2^21=50,331,648`
marked completions. No enumeration of `2^28` hosts is needed or claimed.

Now take the disjoint union of this completion problem and its LITERAL
converse: reverse all local sources and all local observer orders. There
are no global pins here. Converse gives

    G_op(z) = z^24 G(1/z)
            = (1+z)^21 (10z+10z^2+4z^3).

The original input has THREE covered forward arcs, and the converse has
ONE. To weight ALL forward source arcs, these fixed monomials must be kept.
The disjoint-union polynomial is consequently

    F_union(z) = z^3 G(z) + z G_op(z)
               = 2 z^2 (1+z)^22 (5+2z+5z^2).               (9)

It counts `100,663,296` marked records. Its second quadratic has
discriminant `4-100=-96` and roots `(-1 +/- 2i sqrt(6))/5`, both nonreal
and on the unit circle. Thus converse symmetry does not restore real roots.
More precisely (9) is `28`-reciprocal: `z^28 F_union(1/z)=F_union(z)`.
Its actual degree is `26`, with support `2,...,26`; do not call it reciprocal
relative to degree `26`. After removing `z^2`, the degree-`24` factor is
ordinarily reciprocal. The exact length-`29` coefficient vector is

    [0,0,10,224,2408,16544,81620,308000,924616,2267232,
     4626006,7958720,11653840,14615104,15754648,14615104,
     11653840,7958720,4626006,2267232,924616,308000,81620,
     16544,2408,224,10,0,0].

Every deterministic arc in each of these two inputs is a path arc. For
each block order the deterministic arcs form one complete Hamiltonian path,
so `||b_O||^2=2`. Both conditioned inputs retain the universal class law of
Section 3, in particular `E Q=44` and `E c3=31/2`, while their label-dependent
forward-arc polynomial loses real-rootedness. This cleanly separates retained
isomorphism-invariant statistics from changed reference-order root geometry.
There is no asserted harmonic preservation or native routing between the two
source problems; the union in (9) is a union of counting domains.

## 6. The real-rooted unconditioned survivor, proved rather than assumed

With no local cells, the observer polynomial is

    (1+z)^f0 A_n(z),

where `A_n(z)` counts ascents of permutations of `n` named owners. For
completeness its real-root property follows directly from insertion:

    A_1=1,
    A_(n+1)=(1+n z) A_n + z(1-z) A_n'.                     (10)

Indeed a permutation with `k` ascents has `k+1` insertion positions for
the new largest owner preserving `k` ascents, and `n-k` positions increasing
it to `k+1`. This proves (10).

Starting with `A_2=1+z`, suppose `A_n` has its `n-1` distinct negative
roots `r_1<...<r_(n-1)`. At each root,
`A_(n+1)(r_i)=r_i(1-r_i) A_n'(r_i)`, whose sign is the negative of that
derivative's sign. The signs therefore alternate at successive roots.
The monic degree-`n` polynomial `A_(n+1)` has the opposite sign at negative
infinity to its sign at `r_1`, while its value at `0` is `1`, opposite to
its sign at `r_(n-1)`. The intermediate value theorem produces one negative
root in each of these `n` disjoint intervals, exhausting the degree. All
are simple. Induction proves all roots of `A_n` real and negative for
`n>=2`; adjoining the `(1+z)^f0` factor preserves real-rootedness.

The first failed implication in an attempted constrained extension would
be replacing the Eulerian polynomial by an arbitrary endpoint-comparison
block polynomial while assuming its real-rootedness. The four-block witness
in (8) refutes exactly that step. Nonnegative coefficients and even the
converse symmetry of (9) do not repair it.

All roots in this section and the preceding one are zeros of scalar
counting polynomials. They are NOT eigenvalues of a self-adjoint operator.
No spectrum, Laplacian positivity, or harmonic conclusion is inferred.

## 7. Exact closed controls and remaining scope

`check.py` uses only the Python standard library and explicit exception-based
checks that remain active with `python -O`. It independently verifies:

1. Every labelled host and every actual Hamiltonian observer for `n<=5`,
   by both host-first and observer-first enumerations; direct triangle
   counts; all moments, covariance, automorphism-weighted sums, fixed-path
   class fibres, converse-merged checksums, and labelled agreement-path
   injectivity.
2. All coherent zero/one/two-distinct-cell families of arbitrary nonempty cell sizes
   on `n<=4`, every covered source orientation, and every local Hamiltonian
   record. The comparison checks entire local `K` tournaments, not merely
   observer subsequences. These yield `4,854` full-record instances and
   `95,981` pin-conditioned cases, including every single owner/global-clock
   pin and impossible double/out-of-range pins. Branch and cycle controls
   are included. `80,394` pin-conditioned cases have no completion.
3. The `24` hostile block orders and literal-converse orders, exact binary
   multiplicities of all `21` free arcs, both factorizations, full coefficient
   vectors, quadratic discriminants, and reciprocal exponent bookkeeping.
4. Eulerian insertion recurrence against actual permutations through `n=8`.
5. Every row and both weighted checksums of the original order-eight atlas
   using independent `H`, automorphism, converse-isomorphism, degree, and
   directed-triple calculations.

These finite controls support, but do not replace, the all-size proofs.
They do not enumerate large labelled-host spaces, solve LRC(14), establish
new native packet wiring, or turn marked counts into unmarked counts.

`MANIFEST.json` pins the complete closed packet. `verify_packet.py` copies
ONLY its listed files to a fresh temporary directory and reruns the checker
normally and under optimization, comparing both byte-for-byte with the
frozen reports. No parent-directory file, repository checkout, Git state,
network service, third-party package, or import from an inherited checker
is needed for reproduction.
