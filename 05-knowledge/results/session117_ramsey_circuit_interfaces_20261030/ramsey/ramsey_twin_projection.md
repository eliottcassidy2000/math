# Skew doubling: exact twin projection, a 3-state counting algorithm, and complete alternating maxima

**Status:** proved below; independently checked by exact finite enumeration. No general Boolean circuit lower bound or new Ramsey-number bound is claimed. This closes the specific uniform-family gap in `01-canon/theorems/THM-483-zigzag-law-refutes-the-sandwich.md` (general `trans(D(A_l))` previously bracketed between `2l+1` and `2l+2`). No literature-priority assertion.

## Definitions and inherited source

For a nonempty tournament `T` on `V`, let `t(T)` be its largest transitive subtournament order. Its skew double has vertices `(v,0),(v,1)`, with these exact rules:

- copy 0 follows T and copy 1 reverses T;
- for distinct u,v, `(u,0)->(v,1)` and `(u,1)->(v,0)` iff `u->v` in T;
- `(v,0)->(v,1)` for every v.

A twin pair consists of both copies of one vertex. A transitive vertex subset has a unique linear order consistent with its arcs. The inherited THM-483 proves `t(D(T))=z(T)` for its zigzag invariant, at most one twin pair per transitive subset, and `t(T)+1 <= t(D(T)) <= 2t(T)`. We strengthen the twin statement and completely resolve its alternating-family gap.

## 1. Twin projection theorem

**Theorem 1.** Every transitive subset of `D(T)` containing a twin pair projects to a transitive subset of `T`. For each transitive k-vertex subset of T, exactly `k(k+1)/2` transitive subsets of D(T) project onto it and contain a twin pair. All those lifts have size k+1.

**Proof.** Let the twin be `(v,0)<(v,1)` in the doubled chain. Another plain `(x,0)` cannot lie between them: its two required arcs would demand both `v->x` and `x->v`. Another prime `(y,1)` can lie neither before both nor after both, by the same contradiction using the reversed copy-1 rule. Thus the chain has the shape

`P_before, (v,0), Q_primes, (v,1), P_after`.

This also excludes another twin pair, since its prime would be inside and plain outside, producing opposite cross/copy requirements. Project the chain and reverse the order of `Q_primes`. Then

`P_before, v, reverse(Q_primes), P_after`

is a transitive order in T: all cross-block arcs point forward by the cross rule, v dominates Q, and reversing Q corrects the copy-1 rule. This proves projection.

Conversely fix a transitive subset with order `x_1->...->x_k`. For each `1<=p<=q<=k`, take the doubled chain

`(x_1,0),...,(x_p,0),(x_q,1),...,(x_p,1),(x_{q+1},0),...,(x_k,0)`.

All four arc rules verify it. Its unique twin is x_p, and its primed vertices project onto the consecutive interval `[p,q]`. The preceding structural argument shows every twin-containing lift has exactly this form. The number of nonempty intervals is `k(k+1)/2`. QED.

**Corollary 1.1.** For every nonempty T of order n,

`t(D(T)) <= min(2t(T), max(n,t(T)+1))`.

In particular a nontransitive T satisfies `t(D(T))<=n`, while a transitive T satisfies `t(D(T))=n+1`.

Indeed a twin-free subset uses at most one vertex from each of the n pairs; a twin-containing one has size at most t(T)+1. The original `2t(T)` bound remains valid. The transitive case attains n+1 by Theorem 1. The nontransitive hypothesis cannot be dropped.

## 2. Exact generating polynomial and proof-carrying enumeration

Let `F_T(x)=sum_S x^|S|`, over all transitive vertex subsets S of T, including the empty set. Let `Z_T(x)` count only twin-free transitive subsets of D(T), including the empty set. Then Theorem 1 proves the uniform identity

`F_{D(T)}(x) = Z_T(x) + x^2 F'_T(x) + (x^3/2) F''_T(x)`.

Its k+1 coefficient from twin-containing lifts is exactly `binom(k+1,2)[x^k]F_T`. Thus a lossless algorithm is:

1. Enumerate 2^n subsets of T to calculate F_T.
2. Enumerate 3^n paired states: each original vertex is absent, plain, or prime. Test the corresponding doubled subset for transitivity and calculate Z_T.
3. Add the explicit derivative term.

Transitivity can be checked by induced outdegrees: on m vertices, a tournament is transitive iff its outdegree multiset is `{0,...,m-1}`. Consequently the algorithm uses `O(poly(n)3^n)` time and polynomial auxiliary space (apart from output coefficients), rather than naive enumeration of all `4^n` subsets of D(T). This is an explicit reduction relative to the naive doubled-subset baseline, NOT a claim of a best-known exact algorithm. For maximum order alone, return `max(t(T)+1, deg Z_T)`.

Each accepted state carries a directly checkable chain certificate: order vertices by decreasing induced outdegree and check all forward arcs. Rejected states can carry a directed 3-cycle. The arithmetic twin contribution needs only Theorem 1 and the base transitive-subset certificates. Enumeration is exponential; this does not yield a polynomial-size certificate of arbitrary maximum optimality.

## 3. Complete alternating-family classification

For integer l>=1, define A_l with U=`{u_1,...,u_l}` and W=`{w_1,...,w_{l+1}}`:

- u_i -> u_j for i<j;
- w_j -> w_i for i<j;
- w_i -> u_j iff i<=j.

**Theorem 2.** For every l>=1,

`t(A_l)=l+1`, and `t(D(A_l))=2l+1`.

For every l>=2, the unique maximum transitive vertex subset of D(A_l) is

`{(u_j,0):1<=j<=l} union {(w_i,1):1<=i<=l+1}`,

with chain order `w_1',u_1,w_2',u_2,...,u_l,w_{l+1}'`. At l=1, there are exactly 12 maximum transitive subsets, each of order 3.

**Proof of values.** The inherited argument for `t(A_l)=l+1` is short: if a transitive subset selects W-indices with minimum a and maximum b, it cannot select any u_j with `a<=j<b`, since `(w_a,u_j,w_b)` is a directed triangle. Thus it has at most `(l-(b-a))+(b-a+1)=l+1` vertices. If it selects no W, it has at most l. The set `{w_1,u_1,...,u_l}` attains l+1. A_l is nontransitive for l>=1; Theorem 1 gives `t(D(A_l))<=|A_l|=2l+1`. The displayed zigzag chain attains this bound.

**Proof of uniqueness.** When l>=2, a twin-containing chain has length at most l+2<2l+1. Hence a maximum chain selects exactly one copy of every original vertex. Write b_x in {0,1} for the selected copy. For every `a<=j<b`, `(w_a,u_j,w_b)` is a directed triangle in A_l. In the selected doubled triple, arcs are reversed precisely between pairs of prime vertices. A cyclic triple becomes transitive iff exactly two of its vertices are prime. Hence

`b_(w_a) + b_(u_j) + b_(w_b) = 2` for all `a<=j<b`.

Using j=1 and a=1 shows `b_(w_2)=...=b_(w_(l+1))`. Using j=l and b=l+1 shows `b_(w_1)=...=b_(w_l)`. Since l>=2 these equalities overlap, so every W label equals b. Then `2b+b_(u_j)=2` forces b=1 and every U label=0. The displayed chain proves sufficiency.

At l=1, A_1 is a directed triangle. Its three two-vertex transitive subsets each have three twin-containing lifts, giving nine maximum chains. The three full-support choices with exactly two primes give three more. Total 12. QED.

**Uniform-family closure.** THM-483's unbounded defect is now exactly `t(D(A_l))-t(A_l)=l` for every l>=1, rather than only a lower bound with an unresolved one-unit gap. The unique extremizer classification explains why the alternating witness is rigid, not just numerically optimal.

### Closed base-subset enumerator for the entire alternating family

There is also a uniform formula for **every** transitive-subset size in A_l:

`F_(A_l)(x) = (1+x)^(l-1) [1+(l+2)x+binom(l+2,2)x^2]`.

To prove it, select no W (contribution `(1+x)^l`), exactly one W (contribution `(l+1)x(1+x)^l`), or W with distinct extreme indices a<b and span d=b-a. In the last case those extremes are compulsory; the d-1 interior W positions are optional, and exactly l-d U positions remain available, all outside `[a,b-1]`. Every such selection is transitive, ordered U-left, W-descending, U-right. The contribution for fixed extrema is `x^2(1+x)^(d-1)(1+x)^(l-d)=x^2(1+x)^(l-1)`. There are `binom(l+1,2)` extrema pairs. Summing and factoring gives the formula.

In particular A_l has exactly `binom(l+2,2)` maximum transitive subsets. Combined with Theorem 1, the **entire twin-containing size distribution** of D(A_l) is explicitly

`x^2 F'_(A_l)(x) + (x^3/2) F''_(A_l)(x)`.

This does not purport to give the entire twin-free distribution; its top coefficient and degree are classified by Theorem 2.

## 4. Boolean encoding, scope, and non-transfer

For a fixed tournament T and a full-support twin-free assignment b in {0,1}^n, define `H_T(b)=1` iff the selected copy tournament is transitive. This predicate has a polynomial-size Boolean circuit: inspect every triple and reject exactly cyclic triples. For an originally cyclic triple, accept exactly the three assignments of Hamming weight two. For an originally transitive triple with order a->b->c, the only forbidden assignment is `(b_a,b_b,b_c)=(1,0,1)` (the reversed prime endpoints create a cycle). Thus H_T has an explicit conjunction of at most five 3-literal clauses per original triple. This gives O(n^3) gates over any standard fixed fan-in AND/OR/NOT basis up to constant translation.

For A_l with l>=2, Theorem 2 shows H_(A_l) is just the indicator of a single bit pattern, hence computable with O(n) gates. This family is therefore an easy explicit positive control, not a route to superlinear circuit lower bounds. The existential problem for arbitrary induced supports is distinct from evaluating H_T. No claim of NP-hardness, completeness, NP Boolean-family lower bounds, or Ramsey-number record improvement is made.

## 5. Reproduction and hostile controls

Run `python check_ramsey_twins.py` and `python -O check_ramsey_twins.py` (standard library). Exception gates are used instead of Python assertions.

The checker independently scans all doubled subsets, compares their exact transitive-subset polynomial against the 3-state transfer algorithm, checks twin projection and lift multiplicities, for all 1,099 labeled tournaments on 1 through 5 vertices. Every full-support assignment for those tournaments also checks the explicit 3-CNF predicate independently. It checks A_l doubled subsets through l=4, and every full-support copy assignment through l=8. Controls include transitive bases attaining n+1, and l=1 having twelve—not one—maximizers. Those controls forbid the false universal n bound and false l=1 uniqueness extrapolation. Finite checks corroborate, rather than replace, the uniform proofs.
