# Arc-edit grammar on canonical weighted chains

**Session283, 2026-11-12. Status: PROVED statements with FINITE-EXACT controls;
independent review pending.** This is an exact local grammar for the
arc-reversal metagraph, not a finite classification of every pair-free core
and not a Chan–Pak coverage proof.

## 0. Inheritance and scope

[THM-4463](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
§2 proves that the homogeneous-pair graph is a union of directed paths,
whose vertex sets are precisely the maximal transitive modules. The
[weighted-core transport packet](../../weighted_core_transport_20261020/README.md)
already proves the unique representation

`T = Q[TT_(l_1),...,TT_(l_b)]`,

where `Q` has no homogeneous pair, each `l_i>=1`, and weights are taken
modulo `Aut(Q)`. It already gives the alternating one-vertex insertion
recurrence, the exact labelled `(n,b)` insertion kernel and the cycle-index
inversion. The [future-response packet](../../tournament_future_boundary_20261020/README.md)
already implements the cycle-index inversion without building all cores.
None of those mechanisms is claimed new here.

Our new operation is **one arc reversal**. The results below supply a sharp
rank jump, complete weighted-core rank-transition templates, the exact
separator polynomial and its drift sidecar. They contrast with the inherited
insertion law: edit transitions do not close on `(n,b)`, nor even on the
entire separator histogram.

The metagraph used here has tournament isomorphism classes as vertices.
An arc reversal is counted once for every unordered pair in a representative;
parallel transitions and loops are retained. Thus each row has total weight
`binom(n,2)`. This is not the simple unweighted quotient graph. All the
statistics below are invariant under isomorphism.

## 1. Sharp four-step rank bound

For distinct `u,v` let

`h_T(u,v)=#{w outside {u,v}: (u->w) differs from (v->w)}`.

Let `J(T)` join pairs with `h_T=0`, and let `b(T)` be its number of path
components, equivalently the order of the canonical pair-free core.
Since `J(T)` is a forest,

`b(T)=n-|E(J(T))|`.                                             (1)

**Theorem 1.** If `T^e` reverses one arc, then

`|b(T^e)-b(T)| <= 4`.                                          (2)

**Proof.** Write `e={u,v}`. The values `h(a,c)` are unchanged if
`{a,c}` avoids both `u,v`; `h(u,v)` itself is also unchanged. Thus all
changed homogeneous-pair edges touch `u` or `v`, excluding their mutual
pair. In either the old or new path forest there are at most four such
edges. The difference of their counts lies in `[-4,4]`; use (1). QED.

This is sharp already at order five. Reverse `1->3` in the chain
`0->1->2->3->4` with every forward arc. The old core rank is one and the
new one is five. Outgoing masks are respectively
`(30,28,24,16,0)` and `(30,20,24,18,0)`.

## 2. Exact local templates, independent of clone magnitudes

Any reversal in `Q[TT_l]` can be performed without expanding all its
vertices. Split each affected chain immediately before and after the chosen
endpoint(s), keeping the endpoints as singleton blocks, omit empty pieces,
and reverse the arc between those singleton blocks. The resulting template
has at most `b+4` vertices. Each template vertex carries the length of its
transitive piece. Canonical chain contraction of this weighted template gives
the exact child tournament, not only its count. Contraction is valid because
replacing two adjacent homogeneous quotient vertices by transitive chains
makes their union transitive and homogeneous; conversely the inherited
unique chain decomposition detects all possible merges.

This is an exact full-object rewrite **with the split-position sidecar**.
Dropping the positive piece lengths preserves the child's core rank, but
not its ordinary isomorphism class or Hamiltonian path count.

### Internal-chain reversals

**Theorem 2.** In a chain of length `l`, reverse positions `i<j` (zero-based).
The rank change is

- `0` when `j=i+1`;
- `2 + 1_[i>0] + 1_[j<l-1]` otherwise.

Therefore the internal rank-jump polynomial is

`I_l(y)=(l-1)_+ + 1_[l>=3] y^2 + 2(l-3)_+ y^3`
`       + 1_[l>=5] binom(l-3,2) y^4`.                           (3)

**Proof.** Adjacent reversal merely exchanges two consecutive vertices in
the transitive order. For nonadjacent endpoints the affected chain becomes
an ordinal sum of an optional prefix, the three-block cyclic tournament
`C3[1,TT_(j-i-1),1]`, and an optional suffix. The three cyclic blocks have
no homogeneous pair across blocks. No new homogeneous pair can cross into
another old core block: a distinguishing third core block for those two
old core vertices exists because `Q` is pair-free, and it is unaffected by
this internal reversal. The new number of chain blocks is exactly three
plus the two optional pieces. Subtract one. Counting endpoint locations
proves (3). QED.

### Cross-chain reversals

For a chain of length `l`, use endpoint types and multiplicities

- `l=1`: singleton, multiplicity one;
- `l=2`: first and last, each multiplicity one;
- `l>=3`: first and last, each multiplicity one, and interior,
  multiplicity `l-2`.

For each pair of distinct core blocks `i,j` and each pair of endpoint types,
form the split template above. Its core rank, denoted `b_ij(type,type')`,
depends on `Q` and these types but not on the sizes of nonempty split pieces.
To calculate it, replace each selected chain by a chain of length
`min(l,3)`, all other chains by singletons, choose the designated endpoint,
reverse, and contract. This uses at most `b+4` vertices.

**Theorem 3 (complete rank-transition compiler).** The exact Laurent polynomial of
all one-arc rank changes, counted with arc multiplicity, is

`F_(Q,l)(y) = sum_i I_(l_i)(y)`
` + sum_(i<j) sum_(types s,t) m_i(s)m_j(t) y^(b_ij(s,t)-b)`.      (4)

Thus at most `9 binom(b,2)` finite template calculations suffice. In
particular the number of templates is independent of the total size `n`;
integer arithmetic still incurs the bit complexity of the weights. An
ordinary cubic implementation of chain detection gives polynomial-in-`b`
work per template. No enumeration of `n`-vertex tournaments is involved.

**Proof.** Every arc has either two endpoints in one chain, counted by (3),
or endpoints in distinct chains, uniquely classified by their two types.
The split template is an exact transitive blow-up presentation. Its nonempty
pieces, not their magnitudes, determine which chains merge. The number of
arcs of a type is the product of the type multiplicities. QED.

### Universal bulk behavior

If both endpoints of a cross-chain reversal are interior to their chains,
the rank increases by exactly four. Each endpoint loses its two old
homogeneous neighbors and becomes isolated: an exterior prospective neighbor
is distinguished by at least one of the predecessor and successor in the
old chain. A same-chain prospective neighbor is now distinguished by the
other flipped endpoint. No other homogeneous edge changes.

Set `I=sum_i(l_i-2)_+`. Every pair of interior vertices gives a `+4` jump,
except consecutive vertices in the same chain. Thus

`[y^4]F_(Q,l) >= binom(I,2) - sum_i(l_i-3)_+`.                  (5)

For fixed `b` as `n` grows, `I>=n-2b`, and the right side is
`binom(n,2)-O(bn+n)`. Hence a uniformly random arc reversal raises the core
rank by four with probability tending to one, uniformly over all weighted
pair-free cores of fixed order `b`. This does not say repeated flips keep
that fixed core order; they generally do not.

## 3. Exact separator polynomial and mean drift

Let `P_T(z)=sum_(unordered {u,v}) z^h_T(u,v)`.
For a core pair `i,j` define its weighted separator size

`D_ij=sum_(k != i,j, Q(i,k)!=Q(j,k)) l_k`.

The pair-free property gives `D_ij>=1` for distinct core vertices.

**Theorem 4.**

`P_T(z) = sum_i sum_(d=1)^(l_i-1) (l_i-d) z^(d-1)`
` + sum_(i<j) z^D_ij (1+...+z^(l_i-1))(1+...+z^(l_j-1))`.       (6)

**Proof.** Within one transitive chain two vertices at distance `d` are
distinguished precisely by the `d-1` intermediate vertices. Across chains,
assume `i->j`. For endpoint positions `a` in `i` and `c` in `j`, the
separators are the third core blocks of total size `D_ij`, the `l_i-1-a`
vertices after the first endpoint, and the `c` vertices before the second.
Summing gives the two geometric factors. QED.

Write `e_r=[z^r]P_T`. In particular `e_0=n-b`, and

`e_1=sum_i(l_i-2)_+ + #{core pairs {i,j}:D_ij=1}`.              (7)

Equivalently the second term is `sum_(k:l_k=1) q_k(Q)`, where `q_k`
counts pairs of core vertices distinguished only by vertex `k`.

**Theorem 5 (drift).** For `n>=2`,

`sum_e (b(T^e)-b(T)) = 2[(n-2)(n-b)-e_1]`,                   (8)

so division by `binom(n,2)` gives the exact expected one-step rank change.

**Proof.** A homogeneous pair is destroyed by exactly `2(n-2)` possible
arc reversals: either of its two arcs to any exterior vertex. A pair with
exactly one separator becomes homogeneous under exactly two reversals, one
of its arcs to that separator. No pair with two or more separators can
become homogeneous after one reversal. Sum the losses minus gains in the
homogeneous-edge count and apply (1). QED.

When all chain lengths are at least two, (7) simplifies to `e_1=n-2b`.
Then the numerator in (8) is exactly

`2[n^2-(b+3)n+4b]`.                                           (9)

In the opposite extreme, a pair-free tournament has `b=n`, giving mean
change `-2e_1/binom(n,2)<=0`. A pair-free core can only remain at the top
rank or leave it downward in one step. Weighted low-rank families tend
upward; this is a precise asymmetry between inflation and editing.

## 4. Quotient information that is still lost

Unlike the inherited insertion kernel, edit rank does **not** form a Markov
quotient on `(n,b)` under uniform arc reversal. At `n=5,b=4`, outgoing masks

`A=(4,1,2,7,15)`, `B=(8,1,3,6,15)`

give jump histograms respectively

`A: {-3:3,0:4,1:3}`, `B: {-3:1,-1:2,0:4,1:3}`.

Their separator counts `e_1` are six and four, so even the mean differs.

Even the **entire polynomial** `P_T` does not determine the one-step rank
jump distribution. The pair-free tournaments

`C=(16,9,3,5,14)`, `D=(18,8,3,5,14)`

both have `P=3z+6z^2+z^3`, but

`C: {-2:3,0:7}`, `D: {-2:2,-1:2,0:6}`.

The mean agrees as required by (8). The loss is joint dependence: several
homogeneous pairs can be created by the same reversed arc. The endpoint
templates in (4) retain that dependence. Expected-observable closure is
therefore not a full stochastic quotient.

## 5. Connection to exact Hamiltonian counts and its boundary

For any split template `R` carrying positive chain lengths `m_i`, the
inherited [THM-3121 path-cover walk kernel](../../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)
gives an exact Hamiltonian-count evaluator. A transitive `m`-vertex factor
has exactly `S(m,c)` unordered `c`-path covers (Stirling numbers of the second
kind): a set partition specifies its uniquely increasing paths. Therefore

`H(R[TT_mi]) = sum_c [x^c]W_R(x) product_i c_i! S(m_i,c_i)`,
`W_R(x)=1^T X(I-A_R X)^(-1)1`.                               (10)

This supplies an exact joint *structural rewrite and count evaluator*.
It neither promises a fast evaluation at enormous weights nor supplies a
selector for an arbitrary requested odd count. Full core classification,
ordinary-child collision handling and Chan–Pak completeness remain open.

## 6. Reproduction

Run `python3 verify.py` from this directory (or its full path from any cwd).
It independently constructs vertex-level tournaments, recomputes homogeneous
pairs after every arc reversal, and checks the template compiler, separator
polynomial, drift and sharp bound for every labelled tournament through
order five, plus 36 weighted-core cases. It also computes a symbolic
three-million-vertex weighted-triangle example without expanding it.
Normal and `python3 -O` output agree; checks use explicit exceptions rather
than removable assertions. The finite controls support, but do not replace,
the all-order proofs above.
