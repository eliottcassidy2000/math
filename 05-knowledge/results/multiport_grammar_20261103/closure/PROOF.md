# Ordered-end extensions cannot give a reusable tripling grammar

2026-11-03. **PROVED — independent analytic/source/census audit passed.**
The finite controls below are **FINITE-EXACT**. The mathematical claim is
unbounded in suffix size, block size, number of ports, and the additive digit.
It concerns literal named tournaments and ports, not arbitrary signature vectors.

## 1. Inheritance, question, and exact scope

The closest proved mechanism is the seven-coordinate ordered-sum law and
terminal-SCC strongness criterion in
[`session182_physical_returns_20261104/induction/PROOF.md`](../../session182_physical_returns_20261104/induction/PROOF.md).
The count and root-defect coordinates used here are exactly those of that source.
Its infinite `+4` example is the canonical hostile to any claim that nonstrong
suffixes, changing defect, or positive tripling increments are impossible.
The corrected near miss is its explicit warning: an evaluable ordered-component
grammar does not supply the rewrites or return condition needed for both-digit
closure. The least-used relevant sidecar here is the actual parity of the
unrestricted two-slot count, together with Rédei's odd Hamilton-path count,
[THM-001, `redei`](../../../../01-canon/theorems/THM-001-redei.md).

The current concept board is:

1. the full named port state and its joint defect;
2. clock-uniform affine updates, including their intercept;
3. old-port-preserving ordered-end growth;
4. parity as a realizability constraint lost by formal signature algebra;
5. exact integer coverage versus modular reachability.

Other close inherited results were checked before deriving this obstruction:

- [`session183.../mining/PROOF.md`](../../session183_two_branch_frontier_20261104/mining/PROOF.md)
  computes arbitrary contexts for two fixed recipes, rather than arbitrary
  old-state-preserving ordered-end extensions.
- [`multiport_defect_grammar.../COMPONENT_LANGUAGE.md`](../../multiport_defect_grammar_20261103/COMPONENT_LANGUAGE.md)
  decides all word lengths in its fixed actual atom language. Its arbitrary-word
  acceptance test is different from closure by adding a block to a given word.
- [`session154.../closure/NOTE.md`](../../session154_gadgets_20261103/closure/NOTE.md)
  proves the two-child, growth-at-least-three density obstruction. The theorem
  below does not count levels or compare densities: it prevents even two
  consecutive transitions of a specified structural kind.
- The two actual positive controls in
  [`multiport_defect_grammar.../README.md`](../../multiport_defect_grammar_20261103/README.md)
  show exactly how internal insertion or port changes can escape this theorem.

The targeted survey covered those current proof files and their stated
continuation claims. No broad external novelty claim is made.

**Precise operation.** Start with a nonempty actual suffix tournament `B`, an
actual subset `W` of its named vertices, a fixed physical root, and a common
transitive prefix `TT_t`, `t >= 1`. The prefix beats `B`; the root beats the
prefix and `B \ W`; `W` beats the root. An allowed ordered-end extension chooses
an arbitrary nonempty new tournament `D` on disjoint names and arbitrary new
ports `V`, then takes either

\[
 (B',W')=(B\to D,W\cup V)\quad\text{or}\quad
 (D\to B,V\cup W).                                      \tag{1}
\]

All old suffix arcs and all old port memberships are preserved. The old root
and prefix retain their meaning. Both full attachments must be strong.
The new block need not be strong, transitive, bounded, or chosen from a finite
library. Its interior may contain arbitrarily many ports and arcs. A step chooses
its block once, independently of the common clock `t`.

The preserved predicate is the actual adjacency/port relation (1) and strongness
at both ends. Passing to seven coordinates loses individual local guards; these
are not silently recovered. A native Hamilton-complement program satisfying
(1) at its endpoints is covered by the obstruction, regardless of its
intermediate schedule. Conversely, no adjacency construction below is claimed
to be a legal native program without the separate native guards.

## 2. Main theorem and one-step classification

For the attachment described above write `H_t` for its actual Hamilton-path
count. A *clock-uniform tripling step* means

\[
 H'_t=3H_t+d\qquad\text{for every integer }t\ge1,           \tag{2}
\]

with one integer `d`. Equivalently it suffices that (2) holds at any two distinct
prefix lengths with the same `d`.

**Theorem.** No two consecutive extensions of type (1) are clock-uniform
tripling steps. The sides may be chosen adaptively, the two new blocks may be
unrelated, and the two additive digits may be arbitrary integers.

More precisely, for an actual starting signature `(h,p,f,g,j,e,s)`:

| Side | Necessary and sufficient condition | Only possible new block | Digit |
|---|---|---|---|
| Append | `p = 2f` | one new vertex, marked | `2j - 2h - e - g` |
| Prepend | `f = h` | one new vertex, marked | `-3h - s < 0` |

The theorem does not exclude a first append with a nonnegative digit; the finite
census below found none, but no universal append-sign assertion is needed or
made. A one-step tripling example exists, so the length-one boundary is genuine.

The [independent audit](AUDIT.md) checks the unbounded argument and reconstructs
the first-step census by a separate route. Root promoted only the status after
that audit; the producer's mathematical argument and source were unchanged.

The actual root defect is `j` and changes in these steps. It is neither assumed
zero nor replaced by an observer phase. The suffix `B` need not be strong.

## 3. Actual signature and elementary restrictions

For a nonempty tournament and port set, define:

- `h`: Hamilton paths;
- `p`: ordered spanning covers by two directed-path slots, either allowed empty;
- `f`: such covers with first slot nonempty and ending at a port;
- `g`: such covers with second slot nonempty and starting at a port;
- `j`: such covers satisfying both marked conditions;
- `e,s`: Hamilton paths ending, respectively starting, at a port.

The following restrictions are actual, not imposed on an abstract coordinate
cone:

\[
 h\ge1\text{ is odd},\qquad p\text{ is even},\qquad p\ge2h,
 \qquad f\le p-h<p.                                      \tag{3}
\]

Rédei gives the first statement. Swapping the two ordered slots is a
fixed-point-free involution because the tournament is nonempty; this gives even
`p`. The two empty-slot families give `p >= 2h`. The `h` covers with the first
slot empty cannot be counted by `f`, giving its upper bound.

For each marked vertex `w`, split any Hamilton path just after `w`. This is an
injective map into covers whose first slot ends at `w`; concatenation recovers
the original path on this image. Contributions for distinct marked endpoints
are disjoint. Thus

\[
 f\ge |W|h.                                               \tag{4}
\]

Strongness of the full attachment is equivalent to `W` meeting the terminal
SCC of `B`, by the cited S182 proof. In particular `W` is nonempty, so `f >= h`.
For completeness, every vertex of `B` reaches a terminal port, then the root,
then the prefix; the prefix reaches all of `B`. Conversely a terminal SCC
without a port cannot leave itself in the full attachment.

We also need two elementary rigidity statements.

1. `p = 2` forces a singleton tournament. If there are at least two vertices,
   splitting one Hamilton path into two nonempty segments, and swapping the
   slots, supplies two more covers besides the two empty-slot covers.
2. `h = 1` forces a transitive tournament. A nontransitive tournament contains a
   directed triangle, with three Hamilton orders. Each order extends to a full
   Hamilton path by successive vertex insertion without changing the old
   relative order. To insert `v`, prepend it if it beats the first vertex,
   append it if the last vertex beats it, or insert at a transition
   `x_i -> v -> x_{i+1}`. Such a transition exists in the remaining case.
   The three distinct triangle orders remain distinct. This also recovers the
   elementary argument used in the S183 mining proof.

For a transitive tournament on `n` vertices whose terminal vertex is marked,
that terminal vertex alone contributes `2^(n-1)` to `f`: assign every earlier
vertex to either path slot and preserve transitive order. Hence a terminal-port
state with `h = f = 1` must be a single marked vertex.

## 4. Exact one-step calculation

The realized ordered-sum law from S182 is, with upper-case coordinates for the
second block,

\[
\begin{aligned}
h'&=hH,&p'&=pP,&f'&=fH+pF,&g'&=gP+hG,\\
j'&=jH+hJ+eS+gF,&e'&=hE,&s'&=sH.
\end{aligned}                                             \tag{5}
\]

It counts actual ordered covers. In particular the `gF` contribution in the
joint coordinate is retained. The full attachment count and root defect are

\[
 H_t=h-j+2^t f,\qquad \text{root defect}=j.                 \tag{6}
\]

Consequently (2) is exactly

\[
 f'=3f,\qquad h'-j'=3(h-j)+d.                             \tag{7}
\]

**Append.** The terminal SCC of `B -> D` is in `D`, so strongness forces a
nonempty terminal-port set `V` in `D`. Thus `F >= H >= 1`. The slope equation is

\[
 (H-3)f+pF=0.                                             \tag{8}
\]

Since `H` is positive and odd, (8) forces `H = 1`. It then gives `pF = 2f`.
Because `0 < f < p`, the positive integer `F` is less than two, so `F = 1`
and `p = 2f`. The last rigidity statement in §3 makes `D` a single marked
vertex, whose signature is `(1,2,1,1,0,1,1)`. Conversely that block is strong at
the attachment level and satisfies the slope condition whenever `p = 2f`.
Its full state is

\[
 (h,2p,3f,2g+h,j+e+g,h,s),                               \tag{9}
\]

and the digit in the table follows from (7).

**Prepend.** The old terminal SCC is unchanged; new ports are not initially
required. The slope equation now reads

\[
 (P-3)f+Fh=0.                                             \tag{10}
\]

The positive even integer `P` cannot be at least four, so `P = 2`. The new
block is therefore a singleton. If it is unmarked, (10) is impossible; if
marked, it gives `f = h`. Conversely this marked singleton works when `f = h`.
The full state is

\[
 (h,2p,3f,p+g,j+s+f,e,h).                                \tag{11}
\]

Equation (4) and nonempty `W` show that `f = h` forces exactly one old port.
Then `j = 0`, since two disjoint nonempty slots cannot both use that vertex.
Its digit is therefore `2j - 2h - s - f = -3h-s < 0`.

## 5. Why no second step exists

Both possible first steps give

\[
 (h',p',f')=(h,2p,3f).                                   \tag{12}
\]

A second prepend would require `f' = h'`, that is `3f = h`, contradicting
`f >= h > 0`.

A second append would require `p' = 2f'`, that is `p = 3f`. If the first step
was an append, this contradicts its condition `p = 2f`. If the first step was
a prepend, its condition `f = h` makes `p = 3h`, impossible because `p` is even
and `h` odd. This exhausts both adaptive choices, with no size cutoff or digit
restriction, and proves the theorem.

The parity condition is essential to this proof's realized state space. The
formal vector `(1,3,1,1,0,1,1)` permits prepend then append with marginal
`1 -> 3 -> 9` under (5), but its odd `p` cannot arise from a nonempty tournament.
This is a hostile to treating all nonnegative algebraic signatures as actual.

## 6. Sharpness, escape routes, and coverage obligation

The single marked vertex has `(h,p,f,g,j,e,s)=(1,2,1,1,0,1,1)`.
Appending a marked vertex gives `(1,4,3,3,2,1,1)`, so for every `t >= 1`,

\[
 3\cdot2^t-1=3(2^t+1)-4.
\]

Both full attachments are strong. Thus a first uniform tripling step is real;
the theorem blocks its reuse. The defect changes from zero to two.

For the inherited positive controls, let the old suffix be `TT3` with all
three vertices marked. Its signature is `(1,8,7,7,6,1,1)`.

- `TT5` with ports at positions `{0,2,4}` has signature
  `(1,32,21,21,12,1,1)` and gives `3H_t+4`. It adds two terminal vertices
  but changes old vertex 1 from port to nonport, violating old-port preservation.
- `TT2 -> C3`, marking both vertices of `TT2` and one vertex of the cycle,
  has signature `(3,48,21,39,16,1,3)` and gives `3H_t+2`.
  Identify its cycle port with the old terminal vertex of `TT3`. Old ports and
  old suffix arcs are preserved, but the new vertices do not form an ordered-end
  block: a new cycle vertex beats the old terminal vertex.

These are actual adjacency controls, not asserted native-operation programs.
They show that the theorem cannot be enlarged to arbitrary port-changing or
internal constructions. Changing the root, prefix clock, or partition is also
outside (1)–(2).

Any reusable two-branch class with both digits `2,4` at every state necessarily
has a path of two consecutive steps. Consequently it cannot be implemented
solely by (1), even if arbitrary suffix blocks and arbitrary changing defects
are admitted. A successful such grammar must pay an escape coordinate: for
example internal insertion, an old-port change, or a different attachment frame.
This is a necessary design requirement, not a constructive sufficient closure
criterion for those larger languages.

Exact integer spectrum coverage is a separate obligation. An independently
proved closed class still needs seeds and a terminating inverse or another
surjectivity argument for its target integer set. Modular reachability does not
supply this. Even unrestricted two-branch growth from finitely many seeds, with
every child count at least three times its parent's, has the inherited density
obstacle; the present no-two theorem is a
structural obstruction before that counting question arises. Neither result
settles the global Hamilton-path spectrum or excludes other recursive grammars.

## 7. Reproducible controls and scope

The self-contained source is
[`closure.py`](../../../../04-computation/multiport_grammar_20261103/closure.py).
Run from the repository root:

```text
python3 04-computation/multiport_grammar_20261103/closure.py
python3 -O 04-computation/multiport_grammar_20261103/closure.py
```

Its explicit universes are all labeled tournaments of orders 1–5, all actual
port subsets including empty masks for interface controls, and all labeled new
blocks of orders 1–3 with all port subsets. There is no isomorphism quotient or
identification of port sets by a matching marginal. The theorem is analytic;
the finite universes only control definitions, formulas, and sharpness.

Direct permutation enumeration computes all seven coordinates. An independent
subset/end-point dynamic program checks Hamilton counts and the full attachment
formula at two distinct clocks. Reachability independently checks strongness.
The source verifies (5) against literal ordered unions, classifies every allowed
first step in the declared finite universe, retains its actual old-name embedding
and port mask, and rejects every second step with every declared new block.
It also checks the actual positive controls and the nonrealizable parity hostile.

Normal/optimized replay counts and hashes are recorded in `REPLAY.md` after
source freeze. In particular, the observed negative offsets for finite append
cases are reported only as a finite census, not promoted to an unbounded sign
theorem.
