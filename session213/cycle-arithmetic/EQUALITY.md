# All equality cases of the tournament endpoint half-bound

**Status: PROVED; independent adversarial proof audit PASS.** Finite
implementation tests are separate. The unrestricted Hamiltonian-path
spectrum conjecture remains open. No external literature-priority claim.

## 1. Inheritance, precise statement, and the zero boundary

Write `H(T)` for directed Hamiltonian path words and `e_a(T)` for those
ending at a specified vertex `a`. Let `TT_j` be a transitive tournament
of order `j`, allowing `j=0`. Let `R_m`, `m>=3`, have vertices
`v_0,...,v_(m-1)` with all arcs forward in this order except
`v_(m-1) -> v_0`. The ordinal sum `A => B` puts every arc from `A` to `B`.

The endpoint half-bound and its injection are inherited from THM-4876,
`01-canon/theorems/THM-4876-deletion-and-flip-increments-of-the-hamiltonian-path-count-an-increment-of-two-forces-at-most-three-paths-the-increment-ladder-thm-013-for-all-orders-and-what-extension-spectra-forbid-for-chan-pak.md`,
statement 7. Its audit
`05-knowledge/results/chanpak_endcount_audit_D_20261104.md`, under
“Observations (not corrections; candidate sharpenings)”, states the
equality classification only as a census observation through order nine.
That audit suggests the first step `e_O=1`, hence one exit. The missing
step supplied here is an explicit path outside the injection image when
the next vertex has two exits into the transitive remainder.

**Theorem 1.** Let `a` be a non-sink of a finite tournament `T`, and assume
`e_a(T)>0`. Then

    H(T)=2e_a(T)+1                                          (1)

if and only if the marked tournament `(T,a)` is isomorphic to

    TT_j => R_m,       j>=0, m>=3,

with `a=v_(m-2)` in the `R_m` component. For `m=3`, this describes any
marked vertex of the directed triangle, by cyclic relabeling. In every
case

    e_a(T)=2^(m-3),        H(T)=2^(m-2)+1.                  (2)

Equivalently, `a` has one out-neighbor `r`, that vertex `r` also has
outdegree one, and `T-{a,r}` is transitive. To see that this local
certificate gives the displayed structure, order the remainder as a
chain and let `y_j` be the unique out-neighbor of `r`. The chain before
`y_j` is a transitive prefix; its remaining chain vertices followed by
`a,r` form the stated `R_m`. Conversely that structure plainly has the
three certificate properties. This gives a recognition test independent
of counting Hamiltonian paths.

The assumption `e_a>0` matters. Among non-sinks with `e_a=0`, equation
(1) holds precisely when `T` is transitive, since it reduces to `H=1`.
These degenerate cases are not power-of-two endpoint counts.

For a sink, `e_a=H` and (1) cannot hold. Thus the positive-count
classification also applies without explicitly requiring non-sink if
equation (1) is assumed.

## 2. The inherited injection, including its unused paths

Let `O=N^+(a)` and let `E_a` be the paths ending at `a`. For

    P=(x_1,...,x_(n-1),a),

let `i` be the first index with `x_i in O`. Move the final `a` to just
before `x_i`. The result `phi(P)` is a Hamiltonian path: the predecessor
`x_(i-1)`, when present, beats `a` by minimality of `i`, and `a -> x_i`.
The last vertex `x_(n-1)` beats `a`, because it preceded `a` in `P`.
Consequently `phi(P)` ends neither at `a` nor in `O`. It is injective,
since deleting `a` and appending it recovers `P`.

Let `F` be all paths whose endpoint is in `V(T)\({a} union O)`.
This gives the inherited sharpened inequality

    2e_a + sum_{r in O} e_r <= H.                          (3)

It also tells us what equality must mean: when the other endpoints
contribute exactly one path, `phi` must exhaust `F`. A valid path of
`F` whose deletion of `a` is not a path cannot belong to this image.

## 3. Proof for a strong tournament

Suppose `T` is strong and satisfies (1). Every vertex ends at least one
Hamiltonian path: cut a Hamiltonian cycle ending there (the classical
Camion theorem). Since `O` is nonempty, (3) and (1) force

    O={r},                 e_r(T)=1.                      (4)

In particular `a` has the unique exit `r`, and every other vertex beats
`a`. In every Hamiltonian path ending at `r`, `a` is immediately before
`r`: it cannot be last and its successor must be `r`. Therefore deleting
the final `a,r` gives a bijection

    e_r(T) = H(T - {a,r}).                                 (5)

Every path of that deletion extends by `a,r`, since all its vertices
beat `a`. The remainder is nonempty, because a strong tournament here
has at least three vertices.

A tournament has exactly one Hamiltonian path iff it is transitive.
The forward implication follows because a nontransitive tournament
contains a directed triangle, and induced monotonicity of `H` gives
at least three paths. The converse follows from the unique chain order.
Induced monotonicity follows by inserting a new vertex somewhere in
each old path, choosing a position deterministically; deletion recovers
the old path.

By (4)-(5), write `X=T-{a,r}` as a transitive chain
`y_0 -> ... -> y_(ell-1)`.

**The decisive unused path.** Suppose `r` beats two chain vertices
`y_i,y_j` with `i<j`. The word

    (y_0,...,y_j with y_i omitted), a, r, y_i,
    y_(j+1),...,y_(ell-1)                                  (6)

is a Hamiltonian path. Its first block is nonempty and ends at `y_j`;
the connecting arcs are `y_j -> a`, `a -> r`, `r -> y_i`, and then
forward chain arcs. If `j=ell-1`, the last displayed suffix is empty,
and the path simply ends at `y_i`.

Thus (6) belongs to `F`, whose endpoints are in `X`. But deleting `a`
creates the forbidden consecutive step `y_j -> r`, since `r -> y_j`.
It is not in the image of `phi`. This contradicts equality: from
`H=2e_a+1` and `e_r=1`, one has `|F|=e_a`, and the injection from
the equally large set `E_a` must be onto.

Hence `r` beats at most one vertex of `X`. Strongness gives at least
one. That vertex must be `y_0`: otherwise the nonempty initial segment
before it beats all of its complement, so no path reaches that segment
from outside, contradicting strongness. We have proved precisely

    y_0,...,y_(ell-1),a,r

with all arcs forward except `r -> y_0`. This is `R_m`, `m=ell+2`,
with `a` penultimate. This proves necessity in the strong case.

## 4. Reduction to the final strong component

The strong components of a tournament are linearly ordered, all arcs
between two components pointing from the earlier to the later. Every
Hamiltonian path traverses them contiguously in that order. Since
`e_a>0`, `a` lies in the last component `C`; since it is not a sink,
`C` has at least three vertices.

Let `M` be the product of the Hamiltonian-path counts of the earlier
components, with empty product one. Then

    H(T)=M H(C),       e_a(T)=M e_a(C).

Equation (1) becomes `M(H(C)-2e_a(C))=1`. Both factors are integers
and `M` is positive, so `M=1` and `H(C)-2e_a(C)=1`. Every earlier
component has exactly one path, hence is transitive, hence is a
singleton component. Their union is `TT_j`. The strong-case proof
identifies `(C,a)` as the stated marked `R_m`.

## 5. Converse and exact counts

In `R_m`, the transitive word is one Hamiltonian path. Every other path
uses the single backward arc `v_(m-1) -> v_0`, so it has the form

    (S in increasing order), v_(m-1), v_0,
    (the complementary interior vertices in increasing order),

where `S` is any subset of `{v_1,...,v_(m-2)}`. All `2^(m-2)` choices
work. Such a path ends at `a=v_(m-2)` exactly when `a` is omitted from
`S`, giving `2^(m-3)` choices. This proves (2) and (1). A transitive
prefix multiplies both counts by one. The proof is complete.

## 6. Consequences for cycle-to-path conversion

Reverse all arcs to get the starting-endpoint version. If `s_r(T)=c>0`
and `H(T)=2c+1`, then `c` is a power of two; the marked tournament is
`R_m^op => TT_j`, with the corresponding marked vertex. Equivalently,
after reversing the chain labels, it is `R_m => TT_j` marked at the
second vertex of the near-transitive order. A source has `s_r=H`, so
there is no extra positive source exception.

For every positive non-power-of-two endpoint count `c`, the inherited
half-bound therefore sharpens to

    H(T) >= 2c+3

when the marked vertex is not a source (or, dually, not a sink). This
uses oddness of `H` after excluding equality. It is a lower bound,
not a statement that the bound is attained for every `c`.

In particular, take ANY tournament `G` with a sole-exit vertex `x -> r`.
Deleting `x` gives `T` with `HC(G)=s_r(T)=c` (the inherited sole-exit
cycle bijection). If `c>0`, this deletion can have exactly `2c+1`
Hamiltonian paths only when `c` is a power of two. Conversely those
values occur using the near-transitive cores. This is an all-tournament
obstruction to a count-preserving endpoint route, not merely a
statement about the Session212 binary comb.

It does not forbid realizing `2c+1` by an unrelated tournament or by a
conversion that does not preserve the input cycle count as one marked
endpoint count. In particular it does not rule out, or prove, every
odd path count at least 23.

## 7. Reproduction

`python3 equality_verify.py` and `python3 -O equality_verify.py`
check the equivalence on every labeled tournament of orders one
through six, including the zero boundary and the converse orientation.
The candidate structure is recognized directly through strong
components, transitivity, and the two sole exits, not by its counts.
The script separately checks the explicit unused path (6), and the
near-transitive family at larger small orders. These are finite controls;
the classification follows from the preceding argument.

The independent auditor used a separate reverse subset dynamic program
and the equivalent local degree certificate, rather than the author's
strong-component recognizer. Its checks covered all 33,867 labeled
tournaments through order six, all 456 and 6,880 isomorphism classes
at orders seven and eight, and 1,981 constructed unused paths. The
separate audit packet contains its source code and receipts; these counts
are verification scope, not additional theorem assumptions.
