# A constructive class bijection with odd self-converse tournaments

**Status:** PROVED by the explicit construction and group argument below;
FINITE-EXACT checks are separately scoped. This fills a specific structural
omission in the repository's Session22 count-channel packet. It is **not**
a claim of external mathematical novelty. Tournament-isomorphism complexity
equivalence with self-converse tournaments is already in the literature;
the contribution here is the explicit faithful codec and retained symmetry
sidecar connecting the present paired enumeration to that carrier.

## 1. What was already established

The exact inherited source is

`05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/README.md`,

“Verified self-converse reduction, with source correction.” It proves by
twisted Burnside that

```
SC(2n)=F_n(2,0),       SC(2n+1)=F_n(2,1).
```

Its sibling

`05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-review/COROLLARY_REVIEW.md`

proves an equivariant bijection between tournament tuples and one tournament
with disagreement graphs, but explicitly does **not** establish a bijection
with unmarked self-converse tournaments of twice the order. The present
construction supplies such an isoclass bijection, with a vertex mask for
odd order. The prior numerical identity is therefore inherited, not new.

The source **THM-283 + `01-canon/theorems/THM-283-sc-burnside-formula.md`**
contains the anti-cycle-type formula. Use its corrected endpoint-swap
criterion: an antipodal unordered pair can have odd orbit length and
return with its endpoints swapped. It is false that every allowable
unordered-pair orbit must have even length.

## 2. Input and explicit lift

Use the notation and conventions of `ENCODING_AND_ACTIONS.md`. The input is
an aligned rooted consecutive pair, equivalently

```
(T,S,d) on V=[n],       R = disagreement graph of T and S.
```

Here `d_i=1` means the larger input tournament's root points to vertex `i`.
At an unordered pair `{i<j}`, the matrix coordinates are `x_ij=M_ij`,
`y_ij=M_ji`, and `R_ij=x_ij XOR y_ij`.

Construct a tournament `W=L(T,R,d)` on

```
{0} union {p_i : i in V} union {q_i : i in V}.
```

Its arcs are exactly:

1. `p_i -> q_i` for every `i`.
2. `p_i -> p_j` iff `i -> j` in `T`.
3. `q_i -> q_j` iff `j -> i` in `T`.
4. For `i!=j`, `p_i -> q_j` iff `R_ij=1`.
   Since `R` is symmetric, `p_i -> q_j` and `p_j -> q_i` have the same bit.
5. `0 -> p_i` iff `d_i=1`, and `q_i -> 0` iff `d_i=1`.

As usual, the opposite arc is present when an “iff” bit is zero. These rules
cover every unordered pair exactly once and are compatible where written
symmetrically.

Let `tau` fix `0` and exchange `p_i,q_i` for each `i`. It is an
**anti-automorphism**:

```
a -> b in W  iff  tau(b) -> tau(a) in W.
```

This is not an ordinary automorphism. Rule 2 is paired with Rule 3; the two
cross arcs in Rule 4 are paired; the two root arcs in Rule 5 are paired;
and each Rule-1 arc returns to itself with reversed endpoints. Also
`tau^2=id`. Thus every lift is self-converse.

The lift is equivariant for `S_n`: relabeling input vertex `i` as `pi(i)`
relabels both `p_i,q_i` accordingly and fixes `0`. This uses the ordinary
graph action on `R`, not a tournament action on those cross bits.

The number of free bits after normalizing `p_i -> q_i` is

```
binom(n,2) [T] + binom(n,2) [R] + n [d] = n^2.
```

So this is an explicit, lossless lift of the fixed-frame matrix universe,
with its **natural paired action**, to marked self-converse tournaments.
It is not an equivariant reduction of ordinary matrix conjugacy.

## 3. Every odd self-converse tournament has the required marking

Here an anti-involution means an arc-reversing permutation whose square is
identity; this allows the identity permutation on the singleton.

For a tournament `W`, an automorphism preserves every arc. Its automorphism
group `A` has odd order: otherwise an involution in `A` exchanges some two
vertices and reverses their mutual arc, a contradiction.

For self-converse `W`, define the **signed automorphism group**

```
A^± = {(g,epsilon): g preserves arcs if epsilon=0,
                     g reverses arcs if epsilon=1}.
```

The sign map has kernel `A`, and its image is `C_2`; hence `|A^±|=2|A|`.
The signed formulation also covers the singleton: identity with negative
sign reverses its vacuous edge relation. For order at least two, the signs
of automorphisms and anti-automorphisms are already distinguishable from
their underlying vertex permutations.

Because `|A|` is odd, every Sylow-2 subgroup of `A^±` has order two. Its
nonidentity element has negative sign and gives an anti-involution `tau`.
An anti-involution fixes at most one vertex, because two fixed vertices
would have their mutual edge reversed. On odd order it has at least one
fixed vertex. Therefore it has **exactly one fixed vertex**, with all other
vertices paired into 2-orbits.

This argument requires `W` to be self-converse. It does not manufacture an
anti-map for an arbitrary tournament, and it is not a polynomial-time
procedure for finding one in an unmarked input.

## 4. Decoder with a supplied anti-involution

Given `(W,tau)` of order `2n+1`:

1. Name its unique fixed vertex `0`.
2. Enumerate the two-element `tau`-orbits arbitrarily as pairs `i=1,...,n`.
3. In each pair there is a unique arc; name its tail `p_i` and its head `q_i`.
   Then `p_i -> q_i` is automatic. No independent pair-flip gauge remains.
4. Read `T` on the `p` vertices.
5. Set `R_ij=1` iff `p_i -> q_j`. Anti-invariance gives
   `p_i -> q_j iff p_j -> q_i`, so `R` is an ordinary symmetric graph.
6. Set `d_i=1` iff `0 -> p_i`.
7. Recover `S` by reversing `T` on exactly the edges of `R`.

All remaining arcs of `W` are forced by anti-invariance and agree with the
lift rules. This is the inverse of the lift up to the arbitrary enumeration
of the pairs, which is exactly simultaneous input relabeling.

Consequently natural aligned-pair classes are in bijection with **marked**
isomorphism classes `(W,tau)`, where isomorphisms must intertwine the marks.

## 5. Forgetting the anti-involution preserves isomorphism classes

**Theorem.** The construction induces a well-defined bijection

```
{aligned rooted consecutive tournament pairs of sizes n,n+1}/isomorphism
     <->
{self-converse tournaments on 2n+1 vertices}/isomorphism.             (1)
```

**Proof.** Surjectivity follows from Sections 3 and 4. For injectivity,
we first show that any two anti-involutions `tau_1,tau_2` of a fixed `W`
are conjugate under its ordinary automorphism group `A`.

Their signed order-two subgroups are Sylow-2 subgroups of `A^±`, so Sylow
conjugacy gives `g tau_1 g^{-1}=tau_2` for some signed symmetry `g`. If
`g` is positive, it lies in `A`. If it is negative, replace it with
`g tau_1`, which is positive and induces the same conjugation of `tau_1`.
Thus a conjugator can always be chosen in `A`.

Now let `phi:W -> W'` be an unmarked isomorphism between two lifted
objects. The transported mark `phi tau phi^{-1}` is an anti-involution of
`W'`. By the preceding paragraph, an automorphism `a` of `W'` conjugates
it to the chosen mark `tau'`. Therefore `a phi` is an isomorphism of the
marked objects. It sends the unique root to the root and each pair's arc
tail to the corresponding arc tail. It hence induces a permutation of
pair labels preserving `T,R,d`, which is precisely input isomorphism.
QED.

The bijection on isomorphism classes is independent of all decoder choices.
It is not a claim of a distinguished anti-involution or distinguished
vertex labels inside each unmarked tournament.

The same argument, omitting `0` and `d`, gives the even-order companion
between two aligned tournaments on `n` vertices and self-converse tournaments
on `2n` vertices. That corollary uses a fixed-point-free anti-involution.
The current direct executable lift checks below target the odd construction;
the even statement follows from the same proof, not a separate reported census.

## 6. Essential symmetry sidecar: a centralizer, not the whole group

Let `G=Aut(T,S,d)`. Extension of a pair permutation to the two lifted copies
identifies

```
G  ≅  C_{Aut(W)}(tau).                                            (2)
```

Every marked automorphism fixes the root, permutes the `tau`-pairs, and
preserves their uniquely oriented tails; the lift/decoder shows `(2)` in
both directions. But an unmarked automorphism may move `tau` to another
anti-involution. Therefore one must **not** infer `G≅Aut(W)`.

Since all anti-involutions form one conjugacy orbit,

```
# {anti-involutions of W} = |Aut(W)| / |G|.                         (3)
```

In particular this number is odd. Formula `(1)` is a bijection of ordinary
isomorphism classes, **not** an equivalence preserving full automorphism
groups or inverse-automorphism weighted groupoid cardinalities.

**Small hostile.** At `n=1`, both input pairs have trivial automorphism
group. With `d=0` the lift is the transitive tournament on three vertices,
whose automorphism group is trivial. With `d=1` the lift is the directed
3-cycle, whose automorphism group has order three and whose three
anti-involutions are conjugate. Retaining the chosen `tau` recovers the
trivial centralizer in both cases. Forgetting it preserves the two class
identities but not their automorphism orders.

## 7. Enumeration and isomorphism consequences, carefully scoped

The inherited count identity now has a constructive witness:

```
SC(2n+1)=P_n=F_n(2,1).
```

A nonredundant aligned-pair generator can emit a nonredundant complete
self-converse bank by lifting each representative. No doubled-order
canonicalization is needed to establish pairwise distinctness; that is the
content of `(1)`, not a heuristic from matching counts. Computing an ordinary
canonical string for each output is a separate optional task.

The transformation takes `O(n^2)` elementary bit operations in explicit
adjacency storage. Two aligned input pairs are isomorphic iff their lifted
ordinary tournaments are isomorphic, even if the chosen anti-involutions
are omitted from the isomorphism test. This is a polynomial-time reduction
of aligned-pair isomorphism to ordinary tournament isomorphism with output
order `2n+1`.

Conversely, ordinary tournament isomorphism embeds in the aligned-pair
problem through `T -> (T,T,d=0)`. Its lift gives the explicit reduction

```
T -> L(T,empty_graph,all_zero_mask),
```

and two such odd self-converse tournaments are isomorphic iff their input
tournaments are. Thus the odd self-converse restriction is as hard as
ordinary tournament isomorphism under polynomial-time reductions; it is
not automatically an easier isomorphism subclass merely because it has
fewer classes. The reverse reduction is trivial inclusion into the ordinary
tournament-isomorphism problem.

This is **not** a claim that general graph isomorphism reduces to tournament
isomorphism, nor a polynomial-time solution of tournament isomorphism.
Historical self-converse-isomorphism reductions already exist; the session's
history/literature lane supplies their exact source and limits. No novelty
claim is made for the complexity equivalence. Indeed a simpler independent
reduction is the ordinal sum `T -> T + {v} + T^op`: all arcs point from
earlier blocks to later blocks. It is self-converse of order `2n+1`, and
an isomorphism preserves the ordered strong-component decomposition, so
its unique component-prefix of total size `n` recovers `T`. The paired lift
therefore adds a full class-bijection/decoder and symmetry sidecar, not a
previously unavailable proof of isomorphism hardness.

Also do not claim efficient decoding of an arbitrary unmarked self-converse
input without an anti-map algorithm. The singleton is trivial. Otherwise, if an anti-automorphism `f` is supplied,
its order is `2m` with `m` odd, and `f^m` is an anti-involution; its permutation
cycles allow that power to be computed efficiently. Locating such an `f`
from an unmarked tournament is a separate search problem.

## 8. Checks and exact scope

`check_self_converse_lift.py` implements the lift and a decoder accepting
an arbitrary supplied anti-involution. It exhausts:

- all **531** matrix inputs for `n=0,1,2,3` (531 includes `1+2+16+512`),
  verifying tournament validity, anti-invariance, involution, and round trip;
- all **3107** input-matrix/permutation pairs through `n=3`, verifying lift
  equivariance against direct adjacency relabeling;
- all anti-involutions of every lifted object through `n=2`, checking the
  **37** resulting decoded pair classes agree with the originals;
- every vertex relabeling of every lifted object through `n=2`, checking
  **1933** transported-mark decodings.

These are all-objects checks in their stated small universes, not random
samples. `check_self_converse_lift.out.json` records **23,210 active checks**.

The independent generation lane, `generation/check_lift.py`, compares
all lifted aligned-pair classes through `n=3` with complete ordinary
Tournament banks through order seven filtered by an independent converse
isomorphism test. Counts are `1,2,8,88`; the order-seven source bank contains
all 456 ordinary classes. It additionally checks `(2)` and `(3)` on every
lift. These independent checks are not implemented by the decoder above.

Reproduce from this directory:

```bash
python3 check_encoding.py > check_encoding.out.json
python3 -O check_encoding.py > check_encoding.optimized.out.json
cmp check_encoding.out.json check_encoding.optimized.out.json
python3 check_self_converse_lift.py > check_self_converse_lift.out.json
python3 -O check_self_converse_lift.py > check_self_converse_lift.optimized.out.json
cmp check_self_converse_lift.out.json check_self_converse_lift.optimized.out.json
```

Both checkers use explicit `need(...)` calls, so optimization does not erase
their controls. The normal/optimized receipts agree byte-for-byte: the
encoding checker executes 121,942 active checks and the lift checker 23,210.
The receipts are saved alongside their corresponding sources.
