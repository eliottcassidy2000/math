# Independent audit of both all-order spectral directions

**Verdict: VERIFIED, no proof gap found.** Relative to the recovered legal abcde
grammar, the bridge–sink proof establishes forbidden factor ⇒ singular M2;
the vertex-DAG proof establishes avoidance ⇒ permutation-unit-triangular M2.
Their combination proves singularity iff occurrence of be or ee. This audit
makes no additional completeness claim identifying the generated family with
other tournament classes.

Audited documents: `bridge/BRIDGE_SINK_STAR_LEMMA.md` and
`triangularity/AVOIDING_WORD_TRIANGULARITY.md`. No publication was performed.

## 1. Bridge–sink lemma: signs, rows, and boundaries

The kernel vector is precisely
`e_(a,r)-e_(r,a)-e_(s,r)+e_(r,s)` in ordered END-pair coordinates. Its four
coordinates are distinct, since a,s∈R and r∈Q. It is nonzero and integral.

If a has the unique in-neighbor s in R, any Hamilton path ending at s must
start at a; otherwise a needs the already-terminal s as predecessor. Removing
a gives the stated bijection. The condition |R|≥3 ensures the remaining path
has at least two vertices, so its initial vertex is not s and a dominates it.
The antisymmetric two-cover identity has sign `delta_(t=a)L-h(t)`: a singleton
in the second run gives the positive term; a singleton in the first gives the
negative term. Repeated starts are correctly zero.

For the one-start-in-P, one-start-in-R case, Q vertices cannot reappear after
entering R, and all P lies in one initial Hamilton path ending at its sink c.
The sum in equation (2) and its negative H term use the correct run order.
For q=a, comparison reduces both columns to `-F H(t,a)`; for q=s, both reduce
to `-F H(t,s)`. No sign reversal is missing.

When both starts lie in P, fix the two initial P-segments. Their endpoints each
dominate all of R, so the attachment sum really is unconstrained in the two
R starts. Its antisymmetry cancels the two-nonempty-R-run term. When only one
R-run exists, its Hamilton-path total is the same for endpoints a,s. Choices
of whether the other P-segment can end at r contribute coefficients independent
of a versus s. When a start is r, that run must be singleton and the other run
covers P then R. If that other start is in R, its contribution is zero because
P is unreachable; this zero case is implicit and harmless in the document.
Reversing the start order reverses signs. These cases exhaust all ordered rows.

No strongness or Hamilton-cycle uniqueness was silently used in this lemma.
In particular the |Q|=2 case is legitimate: P is a singleton c and the argument
still works. Independent computation below explicitly tests this boundary.

## 2. Birth and last-forbidden-pair persistence

For be, the type-A maximum r is the newest vertex, introduced by a or c just
before an older maximum a. Its parent was s=pred(a), giving bridge s→r and
leaving a with sole old in-neighbor s in R. For ee, the prior type-B pair is
a→r; removing its second/newest maximum r leaves a with its sole in-neighbor
s. The four-vertex seed is a valid ee boundary: a=0, r=1, s=3 and removing r
leaves the directed triangle 0→2→3→0. Thus |R|≥3 and every sink-star hypothesis
holds at either birth.

Freeze the triangle r→b→c→r at the **last** bad pair. Every R vertex has at
least the two incoming neighbors b,c; r has incoming q,c. Hence neither R nor
r contains a maximum, and no R vertex can be the cycle predecessor of one:
its only possible outgoing edge to Q is the bridge into r. All later legal
parents remain in Q.

At birth c has out-neighbors exactly R∪{r}. The first suffix operation cannot
be e (it would form a later ee), so it is c or d and does not choose vertex c.
The first added vertex therefore increases c's indegree to at least two. Before
any hypothetical first later choice of c, its out-neighborhood remains
R∪{r}, whose members all have indegree at least two. Thus c can be neither a
maximum nor a maximum's predecessor. This contradiction proves c is never
selected. The frozen sink remains a sink of Q−r, establishing the same star
kernel throughout the suffix. This does not claim monotonicity of the whole
kernel or of nullity.

## 3. DAG rule and legal avoiding induction

The transformed-arc definition gives no loop, since a loop would come from the
excluded cycle arc. The insertion rule correctly moves old row p from source
q to source x, and sources the new row x at q. All other old outgoing rows are
unchanged. For clarity “old sources other than q” in the rule means all old
arc-source vertices other than q, not merely indegree-zero vertices.

The preservation lemma is sound: after deleting source x, any new directed
cycle must use q. In case (i) q has no incoming edge; in case (ii) its sole
possible incoming neighbor p remains a source because the replacement row q
omits p. In neither case can q lie on a cycle.

The a,c,d cases use an old indegree-one maximum as q. A b step is preceded by
a or c; an avoiding e is either the seed e or follows d. In the nonseed cases,
the previous insertion gives exactly two T-in-neighbors r,p at q, of which p
is the cycle predecessor and is discarded in D_C, while r contributes the
sole D_C incoming edge from s(r)=p. The new p is a D_C source. The seed e has
the same configuration explicitly. Thus the induction covers every legal
avoiding word and every prefix, including the empty word.

## 4. Stronger all-k potential theorem

**PROVED.** Let a finite loopless digraph T on n vertices contain a directed
Hamilton cycle C, and suppose D_C is acyclic. For every 1≤k≤n, the native
ordered k-nonempty-path spanning-cover matrix M_k is permutation-equivalent
to a unit upper-triangular integer matrix. Its bipartite support has exactly
one perfect matching. In particular every M_k is unimodular.

Here rows and columns run over ordered k-tuples of distinct vertices. Choose
integer f strictly increasing along D_C arcs. For any cover with ordered
starts a_i and ends b_i, summing its arc increments gives

    sum_used_(u→v) [f(v)-f(s(u))]
      = sum_i f(s(b_i)) - sum_i f(a_i).

Every summand is nonnegative, and vanishes exactly for cycle arcs. For each
ordered start tuple a, delete the incoming cycle arc at every selected start.
The resulting paths are nonempty and uniquely determined; assign path i to
start a_i. Its endpoint is pred(next selected start after a_i in cyclic order).
Call the resulting ordered endpoint tuple Φ(a).

Φ is a bijection: from any ordered endpoint tuple b, delete all outgoing cycle
arcs at the b_i, and assign the resulting component ending at b_i to path i;
its ordered starts recover the unique inverse tuple. This covers k=1 and k=n,
and allows every singleton path without an exceptional convention.

For a cycle-only cover, `sum_i f(s(Φ(a)_i))=sum_i f(a_i)`. Reindex column b by
Φ^{-1}(b), then order rows and columns by the sum of their start weights, with
the same arbitrary tie order. An entry off the matching diagonal requires a
noncycle arc and therefore strictly increases weight. The diagonal entry is
exactly one, since the cycle-only cover is unique. The result is unit upper
triangular. Its sole supported perfect matching follows by repeatedly forcing
the last row/remaining last column. This proves the theorem without a
cancellation argument or determinant-only inference.

## 5. Independent finite controls

`check_all_k_dag_and_sink.py` imports only the audit's independent recursive
path-growth enumerator, not either production observer implementation.

- It exhausts all digraphs of orders 3 and 4 containing the fixed cycle, keeps
  the 110 with acyclic D_C, and checks all 433 M_k matrices for k=1,…,n.
  Canonical endpoint matching is checked bijective, every diagonal count is
  one, and every nonmatching supported entry has strictly increasing weight.
- It tests 20 arbitrary R configurations of orders 3 and 4 with a's sole
  in-neighbor s, both bridge choices q=a,s, and the boundary Q={r,c}, c→r.
  Every ordered row annihilates the claimed star.

```sh
python3 -B -O session103_words_necklaces/audit/check_all_k_dag_and_sink.py
```

Output: `all_k_dag_sink_receipt.json`. Checks remain active under optimization.
