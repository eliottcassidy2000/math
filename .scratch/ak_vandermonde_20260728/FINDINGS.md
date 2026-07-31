# Arithmetic-Kakeya / projective-matroid / Vandermonde audit

Date: 2026-07-28.  Status labels are deliberately separated below.

## Sources pinned (primary sources only)

- Epoch AI's public arithmetic-Kakeya problem page and PDF:
  `https://epoch.ai/frontiermath/open-problems/arithmetic-kakeya` and
  `https://epoch.ai/files/open-problems/arithmetic-kakeya.pdf`.
- Google DeepMind's public AlphaEvolve problem repository, current `main`
  tree queried through the GitHub API.  Its two similarly named arithmetic
  Kakeya notebooks optimize entropy ratios; they are **not** the six-line
  forcing-certificate verifier from the Epoch prompt.
- Beke--Goh--Hatami--Jaffe--Naylor,
  *A characterization of idempotent Schur multipliers*,
  arXiv:2607.14316v2 (2026-07-27).
- Le--Weber, *Powers of the Vandermonde determinant are eventually
  non-SNP*, and its Lean repository
  `https://github.com/steven-le-thien/vandermonde-snp`, commit
  `b9de864cbb13b5109790b5a8f66865f9632e9659` fetched 2026-07-28.

## Exact public semantics and the unavailable verifier boundary

The public Epoch prompt says that rule (1) may add the labelled difference
between any `e1,e2` with the prescribed first `i` coordinates.  It does not
require their suffixes to agree.  Thus the literal rule supplies the complete
bipartite bank between the two adjacent prefix fibres.  The intuitive
definition instead connects corresponding vertices of adjacent copies and
therefore supplies only the equal-suffix matching.

No source for the deployed six-line verifier is linked from the Epoch page or
PDF.  Consequently:

- **FINITE-EXACT:** the literal mathematical specification is decidable and
  accepts the certificate below;
- **FINITE-EXACT:** the equal-suffix interpretation rejects it;
- **UNKNOWN:** whether a private/deployed verifier follows the literal or
  intended interpretation.

The public AlphaEvolve notebooks do not resolve that unknown: they implement
a different entropy-distribution benchmark.

Literal-rule certificate:

```text
score=14/9 m=10 r=4 n=9 t=0
[(0, 0), (0, 1), (1, 0), (1, 1)]
[3, 3]
[{(1,): (1, 0), (2,): (0, 1)}, {(1, 1): (1, 1), (1, 2): (1, 1), (3, 1): (1, 1), (3, 2): (1, 1)}]
[]
[{(1, 1): (1, 0)}, {(1, 1): (0, 1)}, {(2, 1): (0, 1)}, {(3, 1): (1, 0)}]
```

It belongs to the exact family with dimensions `[D,q]` and score

```text
1 + 1/D + 1/q - 1/(Dq),
```

which tends to one.  This demonstrates why the suffix omission cannot be the
sound interpretation of a certificate-to-`AK(alpha)` theorem.

## New exact structural reduction: projective labels and coloop peeling

Let `L` be the rational row span of all initial seed vectors and rule-(1)
edge vectors.  Passing from the integer span to the rational span loses
nothing for rule (2): any rational witness can be multiplied by a common
denominator, and rule (3) generates all integer linear combinations.

At every vertex make the invertible coordinate change

```text
(a,b) -> (s,d) = (a+b,a-b).
```

For every allowed nonzero label, `s != 0`.  Scaling its generator row by
`1/s` shows that the label retains only

```text
rho = d/s = (a-b)/(a+b) in Q.
```

Conversely, if `rho=p/q`, then the integral label `(q+p,q-p)` realizes it.
Therefore searching arbitrary integer labels is exactly equivalent to
searching rational projective slopes.

Let `B` be the signed incidence matrix of the physical edges and seeds
(a seed is an edge to ground), and let `D_rho` put the slope of each
generator on the diagonal.  In the `(s,d)` columns, the generator matrix is

```text
A = [B | D_rho B].
```

For a current wildcard set `T`, delete both columns of every vertex in `T`.
Rule (2) fires `v` iff the remaining row space contains a nonzero vector
supported only on `d_v`: the target `(c,-c)` becomes `(0,2c)`.  By elementary
linear algebra this is equivalent to

```text
d_v is not in the span of all other live columns,
```

i.e. `d_v` is a **coloop** of the represented column matroid.  The forcing
closure is therefore exactly distinguished-coloop peeling.

Consequences:

1. Label magnitudes and common scaling are irrelevant; only rational
   projective slopes and their exact algebraic relations survive.
2. A forcing proof may be checked by one RREF per round, looking for
   singleton `d_v` pivot rows, or by rank-deletion tests.
3. The known floor `m+r >= n-t` follows again: order simultaneous firings
   arbitrarily.  At each of the `n-t` deletions, removing the coloop `d_v`
   lowers rank by one; hence the initial rank, at most the number `m+r` of
   generators, is at least `n-t`.
4. Search can separate topology from labels: enumerate admissible incidence
   patterns `B`, optimize a finite pool of rational `rho`, use fast
   finite-field matroid peeling as a filter, and replay every record over
   `Q`.  This is implemented experimentally in `ak_matroid_anneal.py`.

For the honest equal-suffix graph, this also proves a sharper boundary.
Every live connected component of a successful certificate must be pinned
by an initial seed or an edge into `T`; otherwise the constant `d`-motion on
that component lies in the kernel and no vertex there can ever fire.  Hence
the incidence block `B` has full column rank `u=n-t`.  If `f` vertices fire
in the first round, their `f` distinguished `d` columns are coloops.
Deleting them lowers rank by `f` and leaves a matrix containing the rank-`u`
block `B`.  Therefore

```text
m+r >= rank(A) >= u+f,     so     f <= (m+r)-(n-t).
```

In particular, **a certificate of score below two cannot fire every live
vertex in its first round**.  At target score `1.675`, no more than `0.675u`
vertices can fire initially.  The same argument applies at later stages
with the number of still-active generator rows.  Thus every genuine sub-two
construction must exploit multi-round grounding.  The 13/7 positive control
does exactly this: its round sizes are `2` and `5`.

The full round-by-round statement is sharper.  At the start of round `j`,
write `u_j` for the live count, `r_j=rank(A_j)`,
`H_j=r_j-u_j`, and `nu_j=2u_j-r_j`.  If the simultaneous firing set has
size `f_j`, first delete its `d`-columns (which loses exactly `f_j` rank)
and let `p_j` be the additional rank lost when its `s`-columns are deleted.
Then every successful paid physical-edge certificate satisfies

```text
rank(B_j)=u_j,                  0 <= p_j <= f_j,
r_(j+1)=r_j-f_j-p_j,           H_(j+1)=H_j-p_j,
nu_(j+1)=nu_j-f_j+p_j,         f_j <= H_j <= H_1 <= (m+r)-u.
```

At termination this telescopes to

```text
r_1=u+sum_j p_j,   H_1=sum_j p_j,
m+r >= u+max_j f_j,   and   score >= 1+1/q
```

for `q` nonempty rounds.  In exactly two rounds,
`p_2=f_2`, `r_1=u+f_2+p_1`, and the score is at least `3/2`.  Moreover any
two-round score at most `5/3` must have the balanced profile

```text
2u-(m+r) <= f_1,f_2 <= (m+r)-u.
```

Thus the candidate `5/3` search shapes `(g,u)=(10,6),(15,9),(20,12)` may
restrict both round sizes respectively to `[2,4]`, `[3,6]`, and `[4,8]`.
If `sigma=(m+r)-rank(A_1)` is row-rank slack, the upper endpoint sharpens
to `(m+r)-u-sigma`.

There is an exact topological form of this headroom.  If `e_j` generator
rows remain nonzero, `C_j=ker(B_j^T)` is their circulation space,
`c_j=e_j-u_j`, and `z_j` rows are bridges, then

```text
sigma_j=e_j-r_j=dim {y in C_j : D_rho y in C_j},
H_j=c_j-sigma_j,
sigma_j >= max(0,2c_j-(e_j-z_j)),
f_j <= H_j <= min(c_j,u_j-z_j).
```

Thus usable forcing headroom is **cycle rank minus slope-compatible
bicirculations**.  Bridges support no circulation and consume vertex-side
capacity; repeated/projectively aligned slopes can turn ordinary cycle
surplus into a bicirculation kernel.  This reframes the target search as a
graphic-matroid/cycle-space problem before exact slope optimization.

There is also an exact topology/slope separation at full row rank.
The polynomial matrix `[B | diag(rho)B]` can have independent rows iff the
paid row set partitions into two `B`-independent sets--two forests after
adjoining the ground vertex for seeds.  Necessarily every row subset `S`
satisfies `|S|<=2 rank(B[S])`.  Conversely, a two-forest partition supplies
a square minor whose determinant has a nonzero squarefree monomial indexed
by one forest, so generic rational slopes make it nonzero.  This gives a
cheap topology-only rejection test before slope annealing.

This is the precise useful Vandermonde connection.  The minor polynomial's
monomials encode two-forest partitions, while special projective slopes can
cancel them and create the bicirculation space.  Thus monomial/Newton support
is valuable for generic feasibility but exact rational replay remains
mandatory for the chosen labels--exactly the support-versus-coefficient
boundary exhibited by the Le--Weber Vandermonde non-SNP family.

Quantitatively, put
`delta_j=max_S(|S|-2rank(B_j[S]))`.  Always
`sigma_j>=delta_j`, while `epsilon_j=sigma_j-delta_j` measures only the
extra slope-specific cancellation.  The exact audits find `epsilon_j=0`
in every round of all four records.  Their rank drops are therefore forced
entirely by local two-arboricity defects created after peeling, not by
fortuitous chosen-label cancellations.  Minimal later dense cores
`(size/rank)` are `7/3` for strict `13/7`; `7/3,7/3` for `7/4`;
`5/2,7/3` for `12/7`; and `6/2` (defect two) for `9/5`.

The proof extends verbatim after vertex identifications if every seed/edge
row is pushed to the quotient and charged at least once: loops become zero
and parallel edges remain paid rows.  Every live quotient component must
still be pinned, or its constant-`d` right-kernel vector prevents any vertex
in that component from ever firing.  This extension audits the quotient
labelled multigraph only; it does **not** prove recursive same-`H` legality,
transport the prefix/suffix label dependence, or settle duplicate/loop cost.
The literal suffix-unconstrained rule lies outside the theorem precisely
because one charged `f_i` entry supplies an uncharged complete-bipartite
bank.

The reduction also identifies the literal exploit's mechanism precisely:
the omitted suffix equality turns an axis matching into a complete bipartite
generator bank.  Its within-fibre differences make all nine distinguished
`d` columns coloops at the initial round.  Under equal suffixes only the
three first-column vertices peel and the process stops.

## Independent computation

Script:

```text
.scratch/ak_vandermonde_20260728/ak_projective_matroid_audit.py
```

Reproduction:

```bash
python3 .scratch/ak_vandermonde_20260728/ak_projective_matroid_audit.py
python3 -O .scratch/ak_vandermonde_20260728/ak_projective_matroid_audit.py
```

The normal and optimized outputs are byte-identical.  Exact results:

- 3x3: literal forces `9/9` in one round; strict peels `(1,1),(2,1),(3,1)`
  and stops at `3/9`;
- the strict `13/7` positive control forces in two rounds of sizes `2,5`;
- family checks: `2x4 -> 13/8`, `3x3 -> 14/9`, `4x4 -> 23/16`,
  `6x6 -> 47/36`; every member forces literally and fails strictly;
- 160 deterministic hostile random instance/mode comparisons agree between
  the projective-coloop test and an independently coded untransformed rank
  criterion;
- the audit prints exact rational generator combinations for all nine
  literal firings.

Hashes at this checkpoint:

```text
audit script  f5159f97717a09dcd49fccd2a28e39a1c535c46309f23e0dbf483eed93fb7f05
normal output 1616fc1fd9d8101429ac1a00389289df915861477f0f5f40061bf6e2ac1f0605
```

The finite-field-filtered annealer was run from the lifted strict `13/7`
positive control:

```bash
python3 ak_matroid_anneal.py 2,2,2,2 25000 11
python3 ak_matroid_anneal.py 2,2,2,2 25000 12
python3 ak_matroid_anneal.py 2,2,2,2,2 12000 21
```

It visited `23,233`, `23,404`, and `11,258` distinct states respectively.
The best exact-replayed score remained `13/7`; no improvement was found.
This is heuristic negative evidence, not a lower bound.  Annealer hash:
`be38d1b846b2d06c89c250006364a1da1da73afa1d88fa13ccb4d74b38889786`.

The exact rank-profile replay is:

```bash
python3 .scratch/ak_vandermonde_20260728/ak_round_rank_profile_audit.py
python3 -O .scratch/ak_vandermonde_20260728/ak_round_rank_profile_audit.py
```

Normal and optimized outputs are byte-identical.  It verifies:

```text
strict 13/7:     r0=13, H0=6, f=(2,5),   p=(1,5)
per-suffix 7/4:  r0=14, H0=6, f=(1,4,3), p=(0,3,3)
quotient 12/7:   r0=12, H0=5, f=(1,3,3), p=(0,2,3)
quotient 9/5:    r0=9,  H0=4, f=(2,3),   p=(1,3).
```

All four frozen witnesses have zero paid-row/rank slack.  The recovered
`9/5` quotient has five classes, one of which contains five of the original
`3x3` grid vertices; its two-round profile `(2,3)` is an exact positive
control for the two-round identities.  All four are initially bridge-free
and bicirculation-free.  Their later active bicirculation dimensions are
respectively `(0,1)`, `(0,1,1)`, `(0,1,1)`, and `(0,2)`.

Rank-profile checkpoint hashes:

```text
profile audit   7009ac65e1eb4dcbafac664a36f60c452c91987a22f28e2e175d98af4984d0f0
profile output  7c7572b0a909a5672ece827839e62d2b6cad9ee25804a79d1029eac04ab7e28c
theorem note    0d590627b0d1f547bdc73d3d8ad5462be5912443cefb53fb34c49b4d28e92df1
```

## Vandermonde-SNP result and its legitimate AK use

The Le--Weber paper proves:

- odd powers are non-SNP from three variables onward by alternation;
- the quadratic power is non-SNP from five variables onward (cited there
  from Monical--Tokcan--Yong);
- for even `k=2m >= 4`, the seed already occurs in `n=k` variables, with
  explicit missing exponent

```text
H=(2m-1)^2,  L=(m-1)(2m-1),  alpha=(H,H,L,...,L).
```

The paper proves permutahedron membership by majorization and coefficient
vanishing by a Dyson constant-term specialization.  The cloned Lean project
has two entry points: `MainJack` is conditional on an explicit structure of
Macdonald--Jack facts, while `MainWard` proves
`even_seed_statement_unconditional` through an internal Ward recurrence.
No `sorry`, custom `axiom`, `admit`, or `unsafe` declaration was found by
source scan.

The exact AK lesson is a guardrail, not a reduction.  The minors controlling
the ranks of `[B|D_rho B]` are polynomials in the slope parameters.  Newton
polytope membership or generic monomial support is not a safe substitute for
exact evaluation: Vandermonde powers give a uniform family in which a lattice
point lies in the Newton polytope while its coefficient cancels.  Hence the
finite-field search filter above must retain exact rational replay and cannot
promote a support/polytope heuristic.

The paper and Lean source contain no `mu_3` fixed branch or `3+3+1=7`
recursion.  Occurrences of fibre/collision language are in the archived Codex
research transcript, not the proved argument.  This rules out the earlier
"Ward recurrence tree" reading of that fragment as a source-backed claim.

## Idempotent-Schur result: useful boundary, not a format compiler

arXiv:2607.14316v2 proves that a Boolean matrix of Schur multiplier
(`gamma_2`) norm at most `gamma` is a signed sum of at most
`2^(C gamma^6)` blocky Boolean matrices.  A blocky mask is a disjoint union
of all-one rectangles with row blocks and column blocks separately disjoint.

For a fixed AK axis/step, arrange physical slots as a prefix-by-suffix
Boolean mask.

- A verifiable/equal-suffix layer has prefix-controlled presence and is
  already rectangular/blocky at the mask level.
- An intuitive same-`H` layer has suffix-controlled presence and is likewise
  already rectangular/blocky.
- An unrestricted mode with both dependencies has a genuinely general mask.

Therefore the Schur theorem does **not** by itself compile the two published
AK formats: their individual support masks already have block complexity
one, while the missing information is which label occupies a slot, whether
copies are identical, and whether vertices are identified.  Moreover a
signed sum of masks is a linear identity, not a union of physical graph
edges; translating it into a forcing graph would need a cost-preserving
gadget for signs and overlapping labels.  The theorem becomes relevant only
after such a sidecar is supplied, or as an analytic complexity bound for a
general prefix-by-suffix mask.

This is a concrete no-go against treating "signed block decomposition" as
the missing equivalence lemma without a compilation construction.
