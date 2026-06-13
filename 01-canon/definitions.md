# Core Definitions

**All Claude instances must use these definitions exactly.** Do not redefine these terms locally in your work. If you believe a definition is wrong or incomplete, open a court case.

Last reviewed: SYSTEM-2026-03-05-S1

---

## Tournaments

**Tournament** T on vertex set V = {1, …, n}: a complete directed graph. For every pair {u,v} ⊂ V, exactly one of u→v or v→u holds. Encoded as T(u,v) = 1 if u→v, with T(u,v) + T(v,u) = 1 for all u ≠ v.

**Opposite tournament** T^op: reverses all arcs. T^op(u,v) = T(v,u).

**H(T)**: the number of directed Hamiltonian paths in T.

**Anti-automorphism** of T: a bijection α: V→V with T(α(u),α(v)) = T(v,u).

**Double-round-robin vertex doubling**: a 2-fiber lift of a tournament `T`
where every vertex `v` is replaced by `(v,0)->(v,1)`, and every base pair
becomes a `2 x 2` block in which the base winner wins one perfect matching
and the base loser wins the complementary perfect matching.  THM-378 proves
that such doublings are exactly voltage lifts `D_sigma(T)` up to sheet labels:
the all-zero voltage is the SC blowup, and sheet-flip gauge classes are
classified by triangle parities.

**SC blowup** `T_SC`: the all-zero double-round-robin vertex doubling.  For a
base arc `u->v`, same-sheet arcs follow `T` (`u_0->v_0`, `u_1->v_1`) and
cross-sheet arcs follow `T^op` (`v_0->u_1`, `v_1->u_0`), with internal arcs
`v_0->v_1`.

**Tournament Analysis**: the repo's general method of converting pairwise data
on `n` objects into a tournament-valued observable.  The data are:

1. a pairwise observable or metric on the `C(n,2)` unordered pairs, possibly
   depending on a parameter;
2. a switch functional deciding one orientation for each pair;
3. a fixed Hamiltonian path used to break ties or to orient threshold bits.

The resulting tournament, or trajectory of tournaments as the parameter
changes, is then studied using `H(T)`, score data, directed cycles, cut
structure, endpoint/packet incidence, and related invariants.

**Tournament Analysis switchboard**: a concrete Tournament Analysis object with
one switch value in `{-1,0,+1}` for every unordered pair, plus a fixed
Hamiltonian tie path.  THM-372 proves that every switchboard determines a
unique tournament, and a finite-wall switchboard determines a piecewise-constant
tournament movie.

**Half-turn circular tournament**: given distinct non-antipodal points on a
circle, orient `x -> y` iff `y` lies clockwise from `x` by distance in
`(0,1/2)`.  THM-374 proves that such a tournament is transitive exactly when
the point set lies in an open semicircle.

**Runner phase clock**: for integer runner speeds `s_i`, the Tournament
Analysis movie obtained from the half-turn comparator on positions
`s_i t mod 1`.  THM-373 proves that its wall times are exactly
`m/(2|s_i-s_j|)` and that it is a closed finite tournament walk.

**LRC observer-source tournament**: for a speed system
`{0,v_1,...,v_{n-1}}` with stationary observer `0` and threshold `1/n`, orient
`0 -> i` iff `||v_i t|| >= 1/n`; orient moving-runner pairs by a complete
tie-resolved comparator such as the half-turn phase clock.  THM-381 proves
that the observer is lonely at `t` iff the marked observer vertex is a source,
and that the source-marked target classes are counted by `A000568(n-1)`.

**LRC trienerment**: a three-state LRC pair comparator on the same circular
positions.  A strict tie means circular distance `< 1/n`; a weak tie means
circular distance `<= 1/n`; non-tied pairs are oriented by a complete circular
comparator such as half-turn order.  THM-389 fixes the boundary convention:
closed-threshold LRC loneliness is equivalent to observer strict-tie-degree
`0`, while the global pigeonhole statement that every `n`-point configuration
has a tie uses weak ties.  With strict ties, the only globally tie-free states
are regular `n`-gon wall states.

**LRC threshold-decorated class fiber**: a tournament-isomorphism class
together with LRC threshold colors, such as which circular gaps are at least
`1/n` or which pair-cells have zero danger deficit.  THM-382 records that in
the bounded `n=3,4` S512 audits, raw and rank-only fibers are mixed while
threshold-decorated gap and pair-cell fibers separate good from bad states.

**LRC pressure lift**: a Tournament Analysis switchboard on runners where an
arc records which runner is the more irreplaceable blocker under deletion.  The
current canon uses three finite variants: `k1` nearest-distance relief, `k2`
two-nearest-distance-sum relief, and two-neighbor threshold-deficit relief.
THM-377 is a selected-row exact certificate for these lifts at `n=14` and
`n=18`.

**LRC endpoint owner**: for an endpoint of the forbidden interval
`||v t|| < 1/n`, the owner is the speed `v` whose interval has that endpoint.
For an interval, the owner is its defining speed.  An interval never strictly
protects its own endpoint, because its own endpoints are boundary points.

**Endpoint-pressure owner graph**: given an endpoint protection core, draw an
edge `u -> v` when an interval owned by speed `u` strictly protects an endpoint
owned by speed `v`.  THM-379 proves that any nonempty owner-compatible core
has a directed owner cycle.

**Pressure-realized endpoint core**: an endpoint protection core whose
endpoint-pressure owner edges are contained in a chosen LRC pressure lift.
THM-380 packages the resulting proof rule: a full-open-cover counterexample
must be sieve-complete, have a nonempty endpoint core, and have a nontrivial
pressure SCC whenever the core is pressure-realized.

**LRC zero-branch star**: for `n >= 2`, `2 <= q <= n`, nonzero q-grid centers
`C subset {1/q,...,(q-1)/q}`, and speeds `S` all divisible by `q`, the local
interval family
`(u/q - 1/(n s), u/q + 1/(n s))` for `u/q in C` and `s in S`.  In the p-adic
application `q=p^d` and `C` is usually the unit set `(u,q)=1`.  THM-391 proves
that every such star has empty strict endpoint-protection core and explicit
peel layers.

---

## The Tiling Model

**Base path** P_0: the fixed path n → n−1 → ⋯ → 1.

**Tile**: a pair (a,b) with a,b ∈ V and a ≥ b+2. There are m = C(n−1, 2) tiles.

**Tiling** t ∈ {0,1}^m: assigns a bit to each tile. Bit 0 means a→b (forward arc); bit 1 means b→a (backward arc). Each tiling determines a unique tournament T_t containing P_0.

**Pin grid**: Grid(n) := {(r,c) ∈ Z² : r ≥ 1, c ≥ 1, r+c ≤ n−1}, where r = a−b−1, c = b. This is isomorphic to the staircase Young diagram δ_{n-2}.

**Strip**: Str(k) := {(r,c) : r+c = k, r,c ≥ 1}, containing k−1 tiles.

---

## The Odd-Cycle Collection Formula

**Conflict graph** Ω(T): vertices are directed odd cycles of T; two cycles are adjacent iff they share at least one vertex.

**Independence polynomial** I(G, x) of a graph G: I(G, x) = Σ_{k≥0} α_k x^k, where α_k = number of independent sets of size k in G. Note α_0 = 1.

**μ(C)** for an odd cycle C ∋ v: μ(C) = I(Ω(T−v)|_{avoid C\{v}}, 2), i.e., the independence polynomial at 2 of the conflict graph of T−v restricted to cycles vertex-disjoint from C\{v}.

---

## The Vertex Deletion Setup

**T−v**: the tournament T with vertex v and all incident arcs removed.

**Claim A** (PROVED — Grinberg-Stanley, arXiv:2412.10572; see CONJ-001): H(T) − H(T−v) = 2 Σ_{C∋v} μ(C), summing over all directed odd cycles C through v.

**Claim B** (proved): I(Ω(T), 2) − I(Ω(T−v), 2) = 2 Σ_{C∋v} μ(C).

---

## The inshat/insact Framework

**P'**: a Hamiltonian path of T−v (i.e., a path through all vertices of V\{v}).

**inshat(v, P')**: the number of ways to insert v into P' to get a Hamiltonian path of T, counting with sign. Specifically, inshat = boundary_term + #Type-I + #Type-II positions.

**Signature** s of P' with respect to v: a binary string of length n−1 where s[j] = 1 if v→P'[j] and s[j] = 0 if P'[j]→v.

**Type-I position** at index j: s[j] = 0, s[j+1] = 1 (arc from successor back to predecessor).

**Type-II position** at index j: s[j] = 1, s[j+1] = 0 (arc forward then backward).

**Boundary term** b = s[0] + (1 − s[n−2]).

**insact(v, P')**: the actual count of valid insertion positions for v into P'. Equal to B_v(P') + S_v(P'), where B_v counts "big" positions and S_v counts "small" positions (precise definitions in the paper).

---

## Key Identities (Reference)

- inshat(v, P') is always odd (proved for all n — inshat = 1 + 2*#TypeII, Lemma 5.3 in paper).
- (inshat(v,P') − 1)/2 = #{Type-II positions in P'} [algebraic identity, THM-004]
- #{Type-II positions in P'} = #{directed 3-cycles (v,a,b) : (a,b) consecutive in P'} [bijection, THM-005]
- insact(v, P') = B_v(P') + S_v(P') [proved, THM verified n≤6]
