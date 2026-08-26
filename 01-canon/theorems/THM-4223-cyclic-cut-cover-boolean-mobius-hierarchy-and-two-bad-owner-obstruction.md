---
id: THM-4223
title: "Cyclic cut-cover Boolean Mobius hierarchy and two-bad-owner obstruction"
status: >
  PROVED / REFUTED / FINITE-EXACT / SUPERSEDED SIDECAR FIREWALL. Proves the all-order
  owner-refined Boolean zeta/Mobius dictionary between cyclic bad-owner sets
  and cyclically ordered path covers, including exact formulas for H,
  endpoints, two-path covers, one-defect counts, and q. Refutes the natural
  local product bound B_ij<=End_i End_j first at a strong tournament of order
  nine. The formerly OPEN repair B_ij<=End_i End_j+min(End_i,End_j) is
  FINITE-EXACT through order nine but REFUTED first at order ten by THM-4224,
  which also refutes every fixed multiple of the min repair. HYP-9081 remains
  OPEN.
source: codex-five-copy-switching-session-20260826
depends_on:
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
  - THM-4219-no-sink-endpoint-energy-floor-and-near-ordinal-sharpness
related:
  - THM-4224-order-ten-minimal-plus-min-two-bad-owner-obstruction
  - THM-3134-tournament-endpoint-jets-and-the-cubic-c3-gregory-newton-profile-transform
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - HYP-9081-strong-tournament-five-copy-endpoint-energy-inequality
script: 04-computation/tournament_cyclic_cut_cover_boolean_mobius_thm4223.cpp
output: 05-knowledge/results/tournament_cyclic_cut_cover_boolean_mobius_thm4223.out
script_sha256: 09d207d641ad0fe320878d2e8b51c4ba9ca209be27b8c14af355ef19e4146262
output_sha256: 90304fbd13dc84ce54a1aa63e6d507c564b4635266dc49309b5fcfa67513c191
hash_basis: raw LF bytes
audit: >
  PASS. A literal 9! permutation scan of the hostile order-nine witness
  independently agrees with subset Hamilton-path DP, checks every displayed
  low-layer formula, and checks Boolean Mobius inversion on all 512 owner
  subsets. Homebrew nauty gentourng -c supplies one representative of every
  strong tournament class through order nine (184,537 classes total); exact
  subset DP verifies H, End_i, N_ij, D_i, q and nonnegative B_ij, finds the
  two product failures displayed here, and finds no +min failure. Compiler
  cross-check details and reproduction commands are frozen in the output.
---

# THM-4223 -- cyclic cut-cover Boolean Mobius hierarchy and two-bad-owner obstruction

**PROVED / REFUTED / FINITE-EXACT / SUPERSEDED SIDECAR FIREWALL.**

The scalar cyclic defect profile in THM-4115 remembers how many bad cyclic
adjacencies occur. This theorem retains their owners. Cutting a cyclic word
at a prescribed owner set gives an exact Boolean transform, not merely a
first-order endpoint inequality. Its first two layers recover the Hamilton
endpoint data and identify the missing two-bad-owner sidecar for HYP-9081.

The hierarchy does not prove HYP-9081. In particular, the most immediate
pairwise injection suggested by the new sidecar is false in order nine. The
slightly enlarged bound that survived this theorem's census was recorded only
as an open question; THM-4224 subsequently refutes it first at order ten and
on an infinite family.

## 1. Owner-refined cyclic cuts

Let `X` be a tournament on a vertex set `V` of order `n>=2`. A **cyclic
ordering** is a cyclic word on `V`, modulo rotation. A cyclic adjacency
`x,y` is good when `x->y` and bad when `y->x`; its **owner** is its left
endpoint `x`. For `S subseteq V`, let

```text
C_S=# cyclic orderings whose bad-owner set is exactly S.       (1)
```

For nonempty `R subseteq V`, let `N_R` count spanning covers by `|R|`
nonempty directed paths, with terminal set exactly `R`, together with a
cyclic ordering of the path components. Put `N_empty=C_empty`. Thus `N_{i}`
is the number `End_i` of Hamilton paths ending at `i`. For two vertices,
cyclic order modulo rotation is unique, so `N_{ij}` is the ordinary number
of spanning unordered two-path covers whose two terminals are `i,j`.

> **Theorem 1 (Boolean cyclic cut-cover hierarchy).** For every tournament
> `X` and every `R subseteq V`,
>
> ```text
> N_R=sum_(S subseteq R) C_S,                              (2)
> C_R=sum_(S subseteq R)(-1)^(|R|-|S|)N_S>=0.              (3)
> ```
>
> In particular, for distinct `i,j,k`,
>
> ```text
> N_ij-N_i-N_j+N_empty>=0,                                (4)
>
> N_ijk-N_ij-N_ik-N_jk+N_i+N_j+N_k-N_empty>=0.            (5)
> ```

Equation `(5)` is the first genuinely higher recutting constraint. The
complete family `(3)`, rather than only the endpoint matrix, is the retained
sidecar.

### Proof

Fix a cyclic ordering `gamma` and cut it immediately after every vertex in
`R`. The result is a cyclically ordered collection of `|R|` nonempty blocks,
each ending at its associated vertex of `R`. Every block is a directed path
if and only if every uncut cyclic adjacency is good. Equivalently,

```text
BadOwners(gamma) subseteq R.                              (6)
```

Conversely, concatenate any cyclically ordered directed path cover counted
by `N_R`. The internal path adjacencies are good, and the only cyclic
adjacencies that can be bad are the component boundaries after vertices of
`R`. Cutting and concatenating are inverse operations. Summing `(6)` over
cyclic orders proves the Boolean zeta transform `(2)`. Boolean Mobius
inversion proves `(3)`; the inequality is the defining nonnegativity of the
literal count `C_R`. Equations `(4)--(5)` are its two- and three-owner
specializations. QED.

## 2. The first two owner layers

Write

```text
c=C_empty,
A_i=C_{i},
B_ij=C_{i,j}=B_ji,        B_i=sum_(j!=i)B_ij.              (7)
```

Thus `c` counts directed Hamilton cycles modulo rotation, `A_i` counts
cyclic orders with sole bad owner `i`, and `B_ij` counts cyclic orders with
exactly the two bad owners `i,j`. Retain THM-4219's notation: `H` is the
Hamilton-path count, `D_i` counts linear permutations with exactly one bad
adjacency, owned by `i`, and `q=sum_i D_i`.

> **Theorem 2 (exact low-layer formulas).** For every tournament,
>
> ```text
> End_i=c+A_i,                                             (8)
> N_ij=c+A_i+A_j+B_ij,                                    (9)
> H=nc+sum_i A_i,                                         (10)
> D_i=(n-1)A_i+B_i,                                       (11)
> q=(n-1)sum_i A_i+2sum_(i<j)B_ij.                        (12)
> ```
>
> Consequently,
>
> ```text
> B_ij=N_ij-End_i-End_j+c>=0,                             (13)
> N_ij>=End_i+End_j-c>=max(End_i,End_j).                  (14)
> ```

### Proof

Equations `(8)--(9)` are `(2)` for owner sets of size one and two. For
`(10)`, cutting an all-good cyclic order at any of its `n` vertices gives
`n` Hamilton paths. A cyclic order with sole bad owner `i` gives exactly one
Hamilton path, obtained by cutting after `i`. Closing any Hamilton path by
its missing cyclic adjacency gives exactly one of these two cases, so this
accounts for every Hamilton path exactly once.

Now rotate a cyclic order into a linear word. An exact-one-bad cyclic order
with owner `i` produces a linear word with sole bad owner `i` when the cut is
at any of the other `n-1` vertices. An exact-two-bad cyclic order with owners
`i,j` produces one such word owned by `i`, by cutting after `j`, and one
owned by `j`, by cutting after `i`. A cyclic order with zero or at least
three bad owners cannot produce a one-bad linear word after one cut. This
proves `(11)`, and summing it proves `(12)`. Finally `(13)` follows from
`(8)--(9)` and `B_ij>=0`; since `c<=min(End_i,End_j)`, it implies `(14)`.
QED.

## 3. Endpoint multigraph and the five-copy target

There is a useful exact packaging of the low layers. Form a weighted
loopless graph `Gamma_X` on `{0} disjoint_union V` with edge weights

```text
w_(0i)=End_i,                    w_(ij)=N_ij.              (15)
```

> **Corollary 3 (clique-star-edge decomposition).** Each all-good cyclic
> order contributes one copy of `K_(n+1)` to `Gamma_X`; each cyclic order
> with sole bad owner `i` contributes one full star centered at `i`; and
> each cyclic order with exactly two bad owners `i,j` contributes only the
> edge `ij`. The weighted degrees are
>
> ```text
> deg(0)=H,                         deg(i)=H+D_i=v_i.       (16)
> ```

### Proof

For one object counted by `c`, equations `(8)--(9)` add one to every edge
in `(15)`. One object counted by `A_i` adds one precisely to `0i` and every
edge `ij`; this is the star at `i`. One object counted by `B_ij` adds one
only to `ij`. This proves the decomposition. Equations `(8)--(11)` give

```text
deg(0)=sum_i End_i=H,

deg(i)=End_i+sum_(j!=i)N_ij
      =nc+sum_j A_j+(n-1)A_i+B_i
      =H+D_i.                                             (17)
```

THM-4219 identifies `H+D_i=v_i`. QED.

Since the total degree is

```text
sum_(x in {0} union V)deg(x)=(n+1)H+q=W+2H=b,             (18)
```

HYP-9081 is exactly the degree-collision inequality

```text
H^2+5sum_(i in V)deg(i)^2
 <=(sum_(x in {0} union V)deg(x))^2.                       (19)
```

The decomposition is stronger than an arbitrary endpoint matrix: it
retains a base clique, singleton-owner stars, and the two-bad-owner edge
measure. It is still not by itself a proof of `(19)`; the next section gives
a sharp obstruction to the most immediate local switching attempt.

## 4. The two-bad-owner product obstruction

The natural attempt to encode each exact-two-bad-owner cyclic order by one
Hamilton path ending at each owner would imply

```text
B_ij<=End_i End_j.                                        (20)
```

This is false even for a strong tournament. Use upper-triangle row-major
encoding, with bit `1` meaning `i->j` for `i<j`. The strong order-nine word

```text
110111111111111111111111011111111110                     (21)
```

has

```text
c=1,
End=(1,1,2,1,5,10,20,5,40),
A  =(0,0,1,0,4, 9,19,4,39).                              (22)
```

For its pairs `{6,7}` and `{7,8}` the exact data are

```text
pair     N_ij     B_ij     End_i End_j
6,7       128      104          100,
7,8       249      205          200.                       (23)
```

Thus `(20)` fails by `4` and `5`, respectively. The second pair attains

```text
B_78=End_7 End_8+min(End_7,End_8)=205.                    (24)
```

> **Finite-exact census 4.** Among one representative of every strong
> tournament isomorphism class through order nine, `(20)` has no failure
> through order eight. At order nine it has exactly the two pair failures
> in `(23)`, both in the class `(21)`. The enlarged inequality
>
> ```text
> B_ij<=End_i End_j+min(End_i,End_j)                       (25)
> ```
>
> has no failure through order nine. This makes nine the first failure
> order for `(20)` in the exact census. Statement `(25)` is **FINITE-EXACT
> THROUGH ORDER NINE / REFUTED ALL-ORDER**: THM-4224 gives its first failure
> at order ten. The finite statement here remains valid.

The strong-class counts for orders `3,...,9` are

```text
1, 1, 6, 35, 353, 6008, 178133,                            (26)
```

or `184,537` classes in total. The script first checks `(21)` by a literal
permutation scan and then uses exact subset DP for the census. It also checks
the formulas `(8)--(13)` on every census row. The stored output gives the
compiler cross-checks and exact reproduction commands.

## 5. Firewall and next switching target

The proved map in this theorem is

```text
cyclic word with exact bad-owner set S
  --cut after R--> cyclically ordered path cover,
```

with preserved predicate `S subseteq R`. The scalar cyclic profile loses
the owner set; the endpoint matrix loses all Boolean layers above two; and
the degree sequence loses the decomposition into clique, stars, and edges.
The exact sidecar restoring those losses is the tensor `(C_S)`, equivalently
`(N_R)` by `(2)--(3)`.

The order-nine witness proves that a proof of HYP-9081 cannot simply inject
each `B_ij` into `End_i times End_j`. THM-4224 further proves that no fixed
multiple of the smaller endpoint family repairs that local injection. A
viable five-copy switching must use the higher `N_R` coherently across
several owners or exploit tournament-specific compatibility among the
clique/star/edge layers. Nothing here or in THM-4224 closes HYP-9081.
