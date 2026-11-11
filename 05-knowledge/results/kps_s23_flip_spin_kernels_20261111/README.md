# Flip and spin for path counts (kind-pasteur-2026-11-11-S23)

Owner directive (2026-11-11): another session on the Chan-Pak conjecture (7 and 21 are the only odd numbers that are
not Hamiltonian-path counts of tournaments), with pasted texts from "a less accurate model" on "flippy" algebras
(split-complex, zero divisors, reflections) versus "spinny" ones (complex, quaternions, octonions), a table
"dimension 2^n <-> tournaments on n + 2 vertices", Cayley-Dickson doubling and the loss of structure across scales,
the scale 5, Fermat primes, A_5, as sources of topics, to be combined with results of the repo.

Canon: THM-5133, HYP-9664, MISTAKE-762. Pre-registration: `PREREGISTRATION.md` (pushed before the runs, commit
b91157e81d). Independent audit: `audit/AUDIT.md`, **PASS WITH CORRECTIONS** (5 majors, 22 minors, all applied; no
mathematical statement of the pre-audit text was wrong except one hypothesis and one boundary case, both listed in
section 8). Chan-Pak remains OPEN. `H(T)` = number of Hamiltonian paths, `S` = the set of values of `H`.

## 0. What comes out

The pasted dictionary is wrong as stated and the scale-5 coincidences are numerology already filed in canon
(section 6). Two statements in it can be made exact, and both lead somewhere:

1. **"Transitive = flip, cyclic = spin" is a theorem about lines.** An algebra with basis `1, e_a` (`a` a nonzero
   vector of `Z_2^n`), `e_a e_b = -e_b e_a = +-e_(a+b)`, *is* a tournament on the points of `PG(n-1, 2)`: `a -> b`
   iff the sign is `+`. The dimension `2^n` corresponds to `2^n - 1` vertices (1, 3, 7, 15), not to `n + 2` (this
   reading is THM-450 (5)). A line `{a, b, a+b}` is a cyclic triangle iff it spans Hamilton's quaternions, and a
   transitive triangle iff it spans the split quaternions, with the one unit of square `-1` as its middle vertex
   (section 3.1).
2. **"Loss of associativity = loss of transitivity" becomes a parity law for 4-sets.** For distinct units,
   `(e_a + e_b)(e_c +- e_d) = 0` for one sign iff `a + b = c + d` and the number of arcs between the two pairs is
   even; a 4-set with `a + b + c + d = 0` has such a pairing iff it contains an even number of cyclic triangles (it
   is transitive or strong), and then exactly two of its three pairings do. The octonion tables are the tournaments
   on the Fano plane with all lines cyclic in which every line complement holds exactly one cyclic triangle. On
   `PG(3, 2)` no tournament has all 4-sets of that kind (linear algebra over `F_2`), at most 8 of the 15 planes can
   have it, the sets of 8 are known exactly (135), the least number of bad 4-sets is exactly 21, and among
   tournaments with quaternionic lines it is 42, the value of the sedenions. On `PG(4, 2)` the doubling is beaten:
   71 octonion planes against the 50 of the Cayley-Dickson algebra of dimension 32 (sections 3.2-3.6).
3. **For path counts, flip versus spin is the sign of a determinant.** A tournament `K` with two vertices replaced by
   arbitrary tournaments `A`, `B` has `H = v_A^T N v_B` for an integer visit matrix `N`. When each port is visited at
   most twice, `N` is 2 x 2, and the sign of `det N` splits the kernels into elliptic (sum of two products),
   parabolic (one product; joins are the degenerate case) and hyperbolic (difference of two products). For a
   closing of THM-5130 with the two ports joined and a nonempty selected set, `det N = (offset)(slope)`: the seeds
   are elliptic, the anti-seeds hyperbolic, and the sign flip of THM-5130 (B) is a change of signature. The first
   hyperbolic kernel has order 5: `H = 3 Pi_2(A) Pi_2(B) - H(A) H(B)` (section 2).
4. **For Chan-Pak.** Every odd `N` has two addresses, `N = 1 + 2^s + d 2^(r+s)` and `N = d' 2^(r'+s') - 2^(s') - 1`.
   A kernel of form `(1, 2^u, 0, d 2^w)` covers a whole plus grid, one of form `(-1, -2^u, 0, d 2^w)` a whole minus
   grid. Among all two-port kernels of order at most 9 the plus carriers are `d = 3, 21, 35` (known) and the only
   minus carrier is `d = 31`: `31 2^(a+b+1) - 2^(a+4) - 1` for all `a, b >= 1` from one kernel of order 8. As path
   counts most of these numbers were already certified one line at a time in THM-5130 (E); new are the uniform
   kernel and the lines `(31 2^r - 1) 2^s - 1` with `r >= 8` (section 2.4). The module-side route to the open
   9-grid (section 4) fails: both pre-registered predictions F1, F2 are refuted, and a congruence modulo 3 explains
   part of it.

Pre-registered outcomes: **F1 FAILED, F2 FAILED, K1 CONFIRMED, K2 CONFIRMED, C1 CONFIRMED, C2 CONFIRMED** (details and
the caveats in `PREREGISTRATION.md`: K1, K2 and C1 were derivable from canon, half of C2 was already in canon).

## 1. Notation and what is known

- `C_T(k)` = number of ordered `k`-tuples of disjoint possibly empty paths covering `T` (cover polynomial; `C_T(1) = H`,
  `C_T(2) = Pi_2`). `P_i(T)` = number of ordered `i`-tuples of disjoint **nonempty** paths covering `T`
  (`P_1 = H`, `P_2 = Pi_2 - 2H`). KNOWN: Session197 bin polynomial, THM-3166, THM-3181; multiplicative under the
  ordered join `A => B`.
- Pair `(B, O)`, state `(r, M, h, E)`, closing `C(B, O; X) = (B => X => v, v -> O)`, closing law
  `H(C(B, O; X)) = r H(X) + M Pi_2(X)`: THM-5130; the law for an arbitrary module is eq. (2) of Session193
  `profile_substitution/PROOF.md` (2026-11-05) and Session216.
- **KNOWN, credited (found by the sweep after the pre-registration).** Working claim 1 of the pre-registration (the
  closing acts on cover polynomials by `k [u_k E + r_k + s_k E^(-1)]`) is Theorem 1 of Session193
  `profile_substitution/PROOF.md` (multi-port transfer, eq. (1)); the 3-cycle operator and the sequence
  1, 3, 15, 123, 1479, 23955 are in Session159 `profile_operator/NOTE.md`, Session197 `path-polynomials/PROOF.md`
  (two positive Jacobi operators, eq. (5)-(13)) and THM-4900. Working claim 2 (`H(C_3[A, B, 1]) = x x' + 2 h h'`) is
  eq. (1) of Session260 `arithmetic/QUASIPERFECT_NORMS.md`, not only its diagonal. The sign tournament of the
  octonions, the action of `e_u -> -e_u` on line orientations, the 16 alternative orientations of the Fano plane
  and the `H` values 189 / 171 are in `arithmetic_braids2_20260917_fano_code.md` and S13-S16 of
  `collatz_mod6_20260922_paley_fano_octonion_design.md`. None of these is claimed here.
- **KNOWN, credited (found by the audit).** The substitution formula of Theorem 2.1 is THM-3121 eq. (4) with `m - 2`
  singleton factors; the two-port expansion and the factorisation with the determinant for one kernel are in
  Session261 `constructors/BILINEAR_KERNEL.md` (`37 H = x x' + 329 H H'`, `329 = det`); `de = N_22` and the
  two-visit criterion are THM-5126 (E). The doubly regular tournaments on `2^k - 1` vertices, `T_15` with
  `H = 198,335,025` and 21 automorphisms, and the reading "a Cayley-Dickson tower at `2^k - 1`" are THM-447, THM-448,
  THM-450 (5), THM-585 (kind-pasteur-2026-06-09-S1). The block recursion of the Cayley-Dickson sign table is eq. (1)
  of Session168 `cayley_dickson/NOTE.md`. The unique grid address is R9 of `seed_recursion_births_20261111/PROOFS.md`,
  which also closes the grids 3 and 21. The anti-seed lines `c 2^s - 1`, `c <= 4097`, are THM-5130 (E).
- Other bilinear laws in canon: s91 `s91_f1_local_densities_20261110/README.md` (one-arc joins, signature (2,1)),
  THM-4875 (Lorentzian gadget adjoint), `multiplicand_saturation_chanpak_20261111/chanpak/OFFSET_TRANSPORT.md`
  section 4 (6).
- **CONCURRENT.** Session263 (`session263_spin_flip_20261111`, same pasted prompt) and Session264
  (`session264_constructive_20261111`) worked the same day on physical source algebras, selectors, Fourier filters
  and multi-return seams; neither I nor the audit found a mathematical overlap with this packet.

## 2. The visit matrix of a two-port kernel

`K` a tournament of order `m`, `x != y` two vertices (ports). `K[A, B]` = the tournament obtained by substituting
`A` for `x` and `B` for `y`.

**Theorem 2.1 (visit matrix; KNOWN: THM-3121 eq. (4)).** Let `N_ij` be the number of walks in `K` that visit every
vertex other than `x, y` exactly once, `x` exactly `i` times and `y` exactly `j` times. For all tournaments `A, B`

    H(K[A, B]) = sum_{1 <= i, j <= m - 1} N_ij P_i(A) P_j(B),

and `N_ij = 0` if `i >= m` or `j >= m`.

*Proof.* Cut a Hamiltonian path of `K[A, B]` into maximal segments inside `A`, inside `B`, and single other vertices.
Contracting the segments gives a walk in `K`; the segments inside `A`, in the order of the path, are an ordered
tuple of disjoint nonempty paths covering `A`, and likewise for `B`. Conversely a walk and two such tuples of the
right lengths glue to a Hamiltonian path, because consecutive symbols of a walk are distinct vertices of `K` joined
by an arc of `K`, and that arc gives every arc between the corresponding parts. The correspondence is bijective.
Between two visits of `x` the walk passes at least two symbols (`x -> u -> x` is impossible in a tournament), at
least one of which is not `y`, and those are used once: at most `m - 2` returns. QED

Check: `ker_verify.py` V1 (220 random cases); audit: 900 cases, walks by plain enumeration; the bound `m - 1` is
attained at every order 3..7.

**Corollary 2.2 (at most two visits; PROVED).** If no walk visits a port three or more times, put
`G = [[N_11, N_12], [N_21, N_22]]`. Then with `h = H`, `pi = Pi_2`:

    H(K[A, B]) = al h h' + be pi h' + ga h pi' + de pi pi',
    de = N_22,  be = N_21 - 2 N_22,  ga = N_12 - 2 N_22,  al = N_11 - 2 N_12 - 2 N_21 + 4 N_22,  al de - be ga = det G,
    N_22 H(K[A, B]) = L(A) L'(B) + det(G) H(A) H(B),   L(A) = N_22 P_2(A) + N_12 H(A),  L'(B) = N_22 P_2(B) + N_21 H(B),

and for transitive modules `H(K[TT_a, TT_b]) = al + be 2^a + ga 2^b + de 2^(a+b)` (`a, b >= 1`). (The four-term
shape for transitive modules is the "module family" of the S21/S22 packets; the factorisation for one kernel is
Session261's.) Checks: V5 (160 cases); audit 282 cases.

**Definition.** A two-port kernel with at most two visits is *elliptic*, *parabolic*, *hyperbolic* as
`det G > 0, = 0, < 0`. When `N_22 > 0`, `N_22 H` is a product plus, or minus, a multiple of `H(A) H(B)`. (Audit: no
kernel of order at most 7 has `N_22 = 0` and `N_12 N_21 > 0`.) Specific to this note: the identity with `det G` for
every kernel, the trichotomy and the census.

**Proposition 2.3 (the sign of the offset is the signature; PROVED).** Let `(B_0, O)` be a pair with state `(r, M)`
and let `K = C(B_0, O; x => y)`, with ports `x, y`. Then both ports are visited at most twice,
`(al, be, ga, de) = (r, 0, 0, M)` and `det G = r M`:

    H(K[A, B]) = r H(A) H(B) + M Pi_2(A) Pi_2(B).

If `O` is nonempty (so `M >= 1`), `K` is elliptic iff the offset `r` is positive and hyperbolic iff it is negative
(`r` is odd, so never 0). If `O` is empty, `M = 0` and `K` is a join.

*Proof.* Closing law with `X = A => B`; `H` and `Pi_2` are multiplicative under the ordered join. No arc goes from
`X` to `B_0` and `v` is used once, so a Hamiltonian path leaves `X` at most once: it meets `X` in at most two
segments, each of which runs through `A` and then through `B`. QED

Examples (V2, V3 of `ker_verify.py`, 120 random module pairs; audit 300):

- **Gaussian kernel** (order 4, the strong tournament `C_3[x => y, 1, 1]`, pair = a point with itself selected,
  `(r, M) = (1, 1)`): `H = Pi_2(A) Pi_2(B) + H(A) H(B)`; on the diagonal `Pi_2^2 + H^2`, discriminant -4. (It follows
  from the two working claims of the pre-registration; so K2 was derivable by hand.)
- `C_3` with one spare vertex: `x x' + 2 h h'`, discriminant -8 (KNOWN, Session260).
- **First hyperbolic kernel** (order 5, `(B_0, O) = (TT_2, both vertices)`, `(r, M) = (-1, 3)`, the anti-seed of
  THM-5130): `H = 3 Pi_2(A) Pi_2(B) - H(A) H(B)`; `G = [[11, 6], [6, 3]]`, `det = -3`, discriminant +12.

So the reversal of the offset under a transitive pair at distance 1 or 2 (THM-5130 (B)) is the passage from a
definite to an indefinite form: in the vocabulary of the pasted text, the closers that reverse the offset are
exactly the ones that turn a "spin" kernel into a "flip" kernel.

### 2.4 Census and grid carriers (FINITE-EXACT)

`kernels.c` (exact visit matrix, orders 3-7) and `kernels_sat.c` (counters saturating at three visits; identical
output at orders 5-7, used for orders 8 and 9; `runs/ker_census_le8.txt`, `runs/ker_census_9.txt`). Counts are over
(class representative, ordered port pair). The audit reproduced every row with its own class lists and programs.

| order | all | at most two visits | join-like (`N_12 = N_21 = N_22 = 0`) | elliptic | parabolic, not join-like | hyperbolic |
|---|---|---|---|---|---|---|
| 3 | 12 | 12 | 6 | 6 | 0 | 0 |
| 4 | 48 | 46 | 12 | 22 | 12 | 0 |
| 5 | 240 | 160 | 26 | 76 | 54 | 4 |
| 6 | 1,680 | 630 | 60 | 290 | 218 | 62 |
| 7 | 19,152 | 2,892 | 162 | 1,242 | 946 | 542 |
| 8 | 385,280 | 16,962 | 602 | 6,732 | 5,008 | 4,620 |
| 9 | 13,790,592 | 148,596 | 3,804 | 53,804 | 39,306 | 51,682 |

No kernel of order at most 4 is hyperbolic; the four hyperbolic kernels of order 5 are the two port orders of two
port pairs of one class, all with `G = [[11, 6], [6, 3]]`. (The order-8 row agrees with the table of S21's
`modfam.c`.) The share of hyperbolic kernels among those with at most two visits is 0, 0, 2.5, 9.8, 18.7, 27.2,
34.8 per cent at orders 3..9. Discriminants `4 (N_12^2 - N_11 N_22)` of the symmetric kernels (`N_12 = N_21`,
`N_22 > 0`) of order at most 6: negative `-4, -8, -12, -16, -20, -24, -32, -36, -44, -72`; zero; positive
`8, 12, 20, 24, 88, 108, 140, 148` (`runs/e2_kernel_signatures.txt`).

**Grids.** Every odd `N > 1` that is not `2^s + 1` is `N+(d, r, s) = 1 + 2^s + d 2^(r+s)` for unique `d` odd,
`r, s >= 1` (R9 of `seed_recursion_births_20261111/PROOFS.md`; the grid of Sessions 260/261), and every odd `N >= 1`
is `N-(d, r, s) = d 2^(r+s) - 2^s - 1` for unique `d` odd, `r, s >= 1` (write `N + 1 = 2^s q`, `q + 1 = d 2^r`). If
`N = 1 (mod 4)` its minus address has `s = 1` and its plus address `s >= 2`; if `N = 3 (mod 4)` it is the other way
round. A kernel of form `(al, be, ga, de) = (1, 2^u, 0, d 2^w)` gives `N+(d, r, s)` for all `s >= u + 1`,
`r >= max(1, w - u + 1)`; a kernel of form `(-1, -2^u, 0, d 2^w)` gives `N-(d, r, s)` for the same ranges.

**Proposition 2.5 (FINITE-EXACT search; each identity PROVED by Corollary 2.2).** Among all two-port kernels of
order at most 9 with at most two visits, the values of `d` with a carrier are:

- plus: `d = 3` (order 8, code 172278145, ports (0, 7): `1 + 2^(a+4) + 3 2^(a+b+3)`), `d = 21` (order 8, code
  71647617, ports (7, 0): `1 + 2^(a+2) + 21 2^(a+b+1)`), `d = 35` (order 9, code 17332994050, ports (8, 1):
  `1 + 2^(a+5) + 35 2^(a+b)`). All three are KNOWN: they are the carriers (16, 24), (4, 42), (32, 35) of
  `OFFSET_TRANSPORT.md` with a one-vertex tail.
- minus: `d = 31` only: order 8, code 17082497, ports (0, 7), `G = [[215, 124], [108, 62]]`, `det = -62`,

      H(K[TT_a, TT_b]) = 31 2^(a+b+1) - 2^(a+4) - 1      for all a, b >= 1,

  so `N-(31, r, s)` is a path count for all `s >= 5`, `r >= 1`: the lines `c 2^s - 1` with `c = 31 2^r - 1`
  (`61, 123, 247, 495, ...`), each for all `s >= 5`; the rows `b = 1, 2, 3` give `27 2^s - 1` (`s >= 3`),
  `29 2^s - 1` (`s >= 4`), `15 2^s - 1` (`s >= 6`). Recount: V4; audit: 25 direct counts and 12 non-transitive
  module pairs.

At order 8 each of the three forms occurs in exactly two classes, a class and its converse (audit). There is no
carrier for `d = 9` on either side.

**What is new in the minus grid (audit, M3).** THM-5130 (E) certifies an anti-seed for every odd `c <= 4097`. The
lines `c = 61, 123, 247, 495, 991, 1983, 3967` (`r <= 7`) are there from the levels 5, 3, 2, 1, 1, 1, 1, and the
rows `c = 27, 29, 15` from the levels 1, 4, 1 (15 also in Session261). So as path counts these were known, from
levels not higher than here. New: one kernel of order 8 with two free transitive modules that carries the whole
grid, and the lines `c = 31 2^r - 1` for `r >= 8`.

## 3. Sign tournaments of twisted group algebras

`V` = the nonzero vectors of `Z_2^n` = the points of `PG(n-1, 2)`. A *sign pattern* is a function
`g(a, b) = -g(b, a) in {+1, -1}` on pairs of distinct points; it is the tournament `a -> b` iff `g(a, b) = +1`. It
is the multiplication table of the imaginary units of any algebra with basis `1, e_a` (`a in V`) and
`e_a e_b = g(a, b) e_(a+b)` for `a != b`; the squares `e_a^2 = +-1` are extra data. The Cayley-Dickson algebras and
their split forms are of this kind. (For the octonions this tournament and its gauge group are in the notes of
2026-09-17 and 2026-09-22 cited in section 1; the tower at `2^k - 1` vertices is THM-448 / THM-450.)

### 3.1 Lines: transitive = split, cyclic = definite

**Proposition 3.1 (line rule; PROVED).** If `1, e_a, e_b, e_(a+b)` span an associative subalgebra, then the line
`{a, b, a+b}` is a cyclic triangle iff the three units have square `-1`; otherwise exactly one of them has square
`-1`, the triangle is transitive, and that unit is its middle vertex.

*Proof.* Let `c = a + b`, `e_a e_b = e_c` (replace `e_c` by `-e_c` if necessary; this reverses all three arcs of
the line and keeps the statement), `e_a^2 = al`, `e_b^2 = be`. Then `e_c^2 = e_a e_b e_a e_b = -al be`,
`e_b e_c = e_b e_a e_b = -be e_a`, `e_c e_a = e_a e_b e_a = -al e_b`. So `a -> b`; `b -> c` iff `be = -1`; `c -> a`
iff `al = -1`. For `(al, be) = (-1, -1)` the triangle is cyclic and `e_c^2 = -1`. For `(-1, +1)`: `c -> a -> b`,
`c -> b`, the unit of square `-1` is `a`. For `(+1, -1)`: `a -> b -> c`, `a -> c`, it is `b`. For `(+1, +1)`:
`a -> c -> b`, `a -> b`, it is `c`. QED

For a generalized Cayley-Dickson algebra (doubling parameters `g_1, ..., g_n`) the map `a -> -e_a^2` is a character
of `Z_2^n` (`e_(i+h)^2 = -g e_i^2`), so the units of square `-1` together with 0 form a subgroup. FINITE-EXACT
(`cd_lines_aut.py`; audit): for every such algebra of dimension 8, 16, 32 the line rule holds on every line. So the
definite algebra has all lines cyclic and every split form has exactly the lines of one hyperplane cyclic (1 of 7,
7 of 35, 35 of 155). The split quaternions are the transitive triangle, the quaternions the cyclic one; the sign
tournament of the split octonions has `H = 45`.

### 3.2 Zero divisors: the parity of cyclic triangles in a 4-set

A *zero-sum 4-set* is `Q = {a, b, c, d}` of distinct points with `a + b + c + d = 0` (an affine plane missing 0; the
complement of a line in the plane it spans). Call `Q` **odd** if it contains an odd number of cyclic triangles
(exactly one: `Q` is a cyclic triangle with a source or a sink) and **even** otherwise (`Q` is transitive or
strong).

**Proposition 3.2 (PROVED; the algebraic half is classical).** For distinct `a, b, c, d in V`,
`(e_a + e_b)(e_c + s e_d) = 0` for some sign `s` iff `a + b = c + d` and the number of arcs from `{a, b}` to
`{c, d}` is even. For a zero-sum 4-set `Q` this holds for no pairing if `Q` is odd and for exactly two of the
three pairings if `Q` is even.

*Proof.* The product is `g(a,c) e_(a+c) + s g(a,d) e_(a+d) + g(b,c) e_(b+c) + s g(b,d) e_(b+d)`. The four indices
can cancel in pairs only as `a + c = b + d`, `a + d = b + c`, that is `a + b = c + d`; then the product vanishes iff
`s = -g(a,c) g(b,d) = -g(a,d) g(b,c)`, that is iff `g(a,c) g(a,d) g(b,c) g(b,d) = 1`: an even number of the four
arcs go from `{c, d}` to `{a, b}`. The number of arcs from `{a, b}` to `{c, d}` is `s_a + s_b - 1` with `s` the
scores inside `Q`. So the pairing is a zero divisor iff `s_a + s_b` is odd. Scores `(1,1,1,3)` and `(0,2,2,2)`
(one cyclic triangle): never. Scores `(0,1,2,3)` (none): the pairings `{0,1}|{2,3}` and `{0,3}|{1,2}`. Scores
`(1,1,2,2)` (two): the two pairings that mix the scores. QED

Checked against the algebra for the Cayley-Dickson algebras of dimension 16 and 32: 84 and 1,260 zero products,
exactly the even pairings and none outside zero-sum 4-sets (`runs/e3_cd_signs.txt`; audit, own doubling).

**Proposition 3.3 (the group of moves; PROVED + FINITE-EXACT).** "Every zero-sum 4-set of a given family is odd" is
a system of affine equations over `F_2` in the arc variables (two independent equations per 4-set). The parity of
every zero-sum 4-set is unchanged by switching at a vertex (reverse all arcs at `u`) and by a line coboundary
(reverse all arcs `{v, w}` with `v + w = f`); rescaling a basis vector, `e_u -> -e_u`, is the product of the two at
`u`. For `n = 3, 4` these moves span exactly the kernel of the parity map, of dimension `2 (2^n - 1) - (n + 1)`
(10 and 25; `cd_orbit.py`): two tournaments have the same parities iff they differ by moves.

*Proof of invariance.* The product of the four cross signs of a pairing changes by
`s_a^2 s_b^2 s_c^2 s_d^2 = 1`, respectively `t_(a+c) t_(a+d) t_(b+c) t_(b+d) = t_(a+c)^2 t_(a+d)^2 = 1`. QED

### 3.3 The Fano plane (FINITE-EXACT, exact linear algebra and full enumeration)

- The tournaments on `PG(2, 2)` whose seven zero-sum 4-sets (line complements) are all odd are exactly `2^10 = 1024`:
  the orbit of the Paley tournament under the moves. Exactly 16 of them have all seven lines cyclic; these are the
  octonion tables (KNOWN: the Hamming coset of 16 alternative orientations). The other 1,008 have 1, 2 or 4 cyclic
  lines (560, 336, 112 of them; audit).
- The other 112 orientations with all lines cyclic have exactly 6 even 4-sets (and `H = 171`, KNOWN). This is the
  second half of prediction C2; it follows in two lines from the gauge action of the cited notes (reversing one
  line of an octonion orientation makes the six other line complements even).
- Call the number of even 4-sets of a plane its **defect**. Over all `2^21` tournaments on the Fano plane the
  defect takes the values 0, 2, 3, 4, 5, 6, 7 on 1024, 21504, 143360, 322560, 688128, 637952, 282624 tournaments.
  **It is never 1.**

### 3.4 PG(3, 2): forced zero divisors and two extremal theorems (FINITE-EXACT)

15 points, 35 lines, 15 planes, 105 zero-sum 4-sets (seven in each plane, each in exactly one plane), 105 arcs.
A plane is **odd** if its defect is 0.

**Theorem 3.4.** (`cd_signs.py`, `cd_more.py`, `cd_extremal.py`, `cd_min.py`, `cd_f5_orbits.py`; audit: own programs
for every item)

1. No tournament on the points of `PG(3, 2)` has all 105 zero-sum 4-sets odd: the 210 equations have rank 80 and
   are inconsistent. So every sign pattern on `Z_2^4` has a zero divisor `(e_a + e_b)(e_c +- e_d)`, whatever the
   squares of the units.
2. At most 8 of the 15 planes are odd in one tournament. The sets of 8 planes that are odd in some tournament are
   exactly 135: the 15 sets "all planes off a point" (type A) and the 120 sets "a plane `P` and the seven planes
   through a point not on `P`" (type B).
3. In type B the other seven planes have defect 6 each, in every solution: 42 even 4-sets. In type A the total
   defect of the other seven planes takes the values 21, 27, 32, 33, 36, 37, 39, ..., 45 (2^19 violation patterns;
   64 of them have total 21, with defect 3 in each of the seven planes).
4. **The least number of even zero-sum 4-sets of a tournament on `PG(3, 2)` is exactly 21** (42 zero products), and
   it is attained only in type A (the `15 x 64` violation patterns of 3).
5. **If all 35 lines are cyclic, the number of even 4-sets is `6 (15 - number of odd planes) >= 42`**, with equality
   iff 8 planes are odd; it is attained by the sign tournament of the sedenions (type B: the plane `{1, ..., 7}` and
   the planes through the point 8) and by tournaments of type A.

*Proof of 4.* Let `o` be the number of odd planes and `E` the number of even 4-sets, the sum of the defects of the
`15 - o` other planes, each at least 2 (section 3.3). If `o <= 4` then `E >= 22`. If `5 <= o <= 7`, let `k` be the
number of planes of defect exactly 2; the exhaustive search of `cd_f5_orbits.py` (the 2,835 sets of five planes
that can be odd together fall into 3 orbits under linear maps; in each, every consistent choice of further planes
of defect exactly 2, with all 21 possible pairs of lines, has at most three planes) gives `k <= 3`, so
`E >= 2k + 3 (15 - o - k) >= 24 - 3 = 21`. If `o = 8`: type A gives `E >= 21` and type B gives `E = 42` by 3.
`o >= 9` is impossible by 2. The value 21 is attained in type A. QED

Second proof (audit, independent of the lemma on planes of defect 2 and of the orbit argument): all violation
patterns next to every one of the 2,835 feasible sets of five odd planes were enumerated (4.6 x 10^11 patterns);
the least number of even 4-sets with at least five odd planes is 21, with `o = 7` forcing at least 42, `o = 6` at
least 25 and `o = 5` at least 27. Hence the "only in type A" of statement 4.

*Proof of 5.* A plane with all lines cyclic has defect 0 or 6 (section 3.3), and at most 8 planes are odd. Both
types have solutions with all lines cyclic (solution spaces of dimension 17 and 11). QED

A minimiser of 4 (`runs/cd_minimum.txt`) has only 6 cyclic lines, irregular scores and `H = 46,700,719`. So on
`PG(3, 2)` the doubling is extremal among sign patterns with quaternionic lines (tied with type A), and not among
all sign patterns.

### 3.5 The Cayley-Dickson sign tournaments CD_n

Standard doubling `(a, b)(c, d) = (ac - conj(d) b, d a + b conj(c))`; `CD_n` = the sign tournament on `2^n - 1`
vertices.

**Proposition 3.5.**

1. (KNOWN: Session168 `cayley_dickson/NOTE.md` eq. (1), in Boolean form; checked `n <= 5`, audit `n <= 6`) Block
   recursion of the sign table `G_n = (g(i, j))` of order `2^n`: `G_(n+1) = [[G, G^T], [G D, -G^T D]]` with
   `D = diag(1, -1, ..., -1)`.
2. (PROVED by induction from 1) Every line of `CD_n` is cyclic; so `CD_n` is regular.
3. (KNOWN: B. Kivunge, "Sign matrices for frames of 2^n-ons under Smith-Conway and Cayley-Dickson multiplication",
   Mile High Conference on Nonassociative Mathematics 2013, for this doubling; short proof below) `G_n` is a
   Hadamard matrix; hence `CD_n` is a doubly regular tournament for every `n >= 2`. Prediction C1 was the case
   `n = 4`, which is KNOWN in canon: `CD_4` is isomorphic to `T_15` of THM-448 (audit: explicit isomorphism), whose
   double regularity, `H = 198,335,025` and 21 automorphisms are THM-447 / THM-448 / THM-585.
4. `CD_2 = C_3`, `CD_3` = the Paley tournament on 7 points (KNOWN). The planes of `CD_4`: 8 odd planes, each a Paley
   tournament (`H = 189`), and 7 planes of defect 6, each the regular tournament with `H = 171`; 42 even 4-sets,
   84 zero products on 42 pairs `{a, b}`. `CD_5`: 50 odd planes, 105 of defect 6, 630 even 4-sets. `CD_6`: 6,510
   even 4-sets. CITED for every `n`: Wilmot, arXiv:2505.11747, Table 5 and Theorems 8, 9, 11 (8, 50, 310, 2046, ...
   octonion and 7, 105, 1085, 9765, ... quasi-octonion subalgebras; `(1/2) N (N-1)(N-3)` zero-divisor pairs with
   `N = 2^(n-1) - 1`, which is `84 [n-1 choose 3]_2`); here odd planes = octonion subalgebras and planes of defect
   6 = quasi-octonion subalgebras. Recomputed here for `n <= 6`. Specific: the tournaments of the planes (Paley and
   the `H = 171` class).
5. (FINITE-EXACT) Split forms of dimension 16 have 42, 78 or 90 even 4-sets.
6. (audit) `CD_5` is not isomorphic to `T_31` of THM-448 nor to its converse, although both are doubly regular with
   21 automorphisms: the two towers part at 31 vertices.

*Proof of 2.* Lines inside the lower half are cyclic by induction. Let `h = 2^n`. For `i != 0`: `g(i, h) = g(0, i)
= 1`, `g(h, i + h) = -eps_i g(i, 0) = 1`, `g(i + h, i) = g(i, i) eps_i = 1`: `i -> h -> i + h -> i`. For a cyclic
line `i -> j -> k -> i` of the lower half: `g(i, j + h) = g(j, i) = -1`, `g(j + h, k + h) = -eps_k g(k, j) = -1`,
`g(k + h, i) = g(k, i) eps_i = -1`: `i -> k + h -> j + h -> i`; the other lines `{j, k + h, i + h}`,
`{k, i + h, j + h}` follow by rotating the names. QED

*Proof of 3.* If `G G^T = G^T G = N I` then the two block rows of `G_(n+1)` have norms `G G^T + G^T G = 2N I` and
`G D D G^T + G^T D D G = 2N I`, and their product is `G D G^T - G^T D G`, which vanishes because
`D = 2 E_00 - I` and the first row and the first column of `G` are all ones: both terms are `2J - N I`. For
`i != j` nonzero, orthogonality of rows `i, j` and regularity give `N/4 - 1` common out-neighbours. QED

Other literature: 84 zero products in the sedenions (Moreno; de Marrais's 42 assessors, arXiv math/0011260); 8
octonion and 7 quasi-octonion subalgebras of the sedenions (Cawagas; Wilmot). Of these only Wilmot's tables and
theorem statements and Kivunge's abstract were read. **Not found in the repo; literature not searched:** the
tournament formulations of Propositions 3.1-3.3 beyond the octonion case, and Theorem 3.4.

### 3.6 PG(4, 2): the doubling is not extremal (added after the pre-audit checkpoint; refutes the first HYP-9664 (C))

**Proposition 3.6 (pullback; PROVED).** Let `n >= 4`, `W` a subspace of `Z_2^n` of codimension 3 and
`p: Z_2^n -> Z_2^n / W = Z_2^3` the projection. Take an octonion tournament `O` on the seven nonzero vectors of the
quotient (all lines cyclic, all line complements odd), put `a -> b` whenever `p(a) -> p(b)` in `O` (`p(a), p(b)`
distinct and nonzero), and orient the remaining arcs so that every remaining line is cyclic. This is possible; the
tournament has all lines cyclic and its `8^(n-3)` planes complementary to `W` are odd.

*Proof.* A line whose three points have distinct nonzero projections projects onto a line of the quotient and is
cyclic because `O` has cyclic lines. Every other line has a point in `W`, and then its other two points have the
same projection (or all three are in `W`): none of its three arcs has been fixed. Two points lie on one line, so
these lines have pairwise disjoint sets of arcs and can be oriented cyclically one by one. A plane complementary to
`W` is mapped by `p` bijectively and linearly onto the quotient, so its tournament is a copy of `O`. The complements
of `W` are the graphs of the linear maps from one complement to `W`: `2^(3(n-3))`. QED

For `n = 4` these are the tournaments of type A of Theorem 3.4: the pullback ties with the doubling. For `n >= 5`
it wins: 64 odd planes against 50 at `n = 5`, 512 against 310 at `n = 6` (with all lines cyclic: 546 and 5,298
even 4-sets against 630 and 6,510; audit, explicit solutions recounted).

**Proposition 3.7 (FINITE-EXACT; explicit certificate `data/pg42_71_odd_planes.txt`).** On `PG(4, 2)` there is a
tournament with all 155 lines cyclic in which the 64 planes disjoint from a line `W` and the 7 planes containing `W`
are odd: 71 odd planes and `6 * 84 = 504` even zero-sum 4-sets. `CD_5` has 50 odd planes and 630 even 4-sets, and no
further plane can be made odd next to its 50. So among sign patterns with quaternionic lines the Cayley-Dickson
algebra of dimension 32 is not extremal: there is one with 71 octonion planes and 1,008 zero products
`(e_a + e_b)(e_c +- e_d)`, against 50 and 1,260.

The seven extra planes are `W` together with one fibre of `p`; all their arcs are free in the pullback. `cd_pg42.py
build` solves the `F_2` system (rank 420 of 465; no further plane can be added); `cd_pg42.py verify` recounts lines,
planes and 4-sets directly from the arc word, with no linear algebra and no import; the audit recounted the
certificate with its own program. Random greedy orders of the planes reach at most 63 (`runs/cd_pg42.txt`).
Bounds for the largest number of odd planes on `PG(4, 2)`: at least 71; at most 82 (audit: each of the 31 solids
holds at most 8 odd planes and each plane lies in 3 solids).

## 4. Fermat modules (E1): two refuted predictions

A *module-side source*: for a seed `(B, O)` of slope `M` and any tournament `X`,
`H(C(B, O; X => TT_t)) = H(X) + M 2^t Pi_2(X)`. With `H(X) = 2^s + 1` and `Pi_2(X) = p 2^e` (`p` odd) this would
serve the grids `d = p c` at level `s`.

- **Congruence (KNOWN, THM-4756 steps 2-3).** `Pi_2(X) = (-1)^n H(X) (mod 3)`. So `3 | H` iff `3 | Pi_2`: for
  `H = 2^s + 1` and odd part of `Pi_2` dividing 9, `s` even forces `Pi_2` to be a power of two and `s` odd forces
  the odd part 3 or 9. Checked on all 198,948 classes of order at most 9.
- **F1 FAILED.** Among the classes of order at most 9 such `X` exist only at levels `s = 1, 2, 3, 7`
  (`C_3`; the strong 4-tournament; `C_3 => C_3`; four classes of order 7 with `H = 129`, `Pi_2 = 3 2^9`): 2 of the 9
  levels 3..11, against the predicted "at least 7".
- **F2 FAILED.** No one-apex word has `H = 2^s + 1` (`3 <= s <= 21`) and odd part of `Pi_2` dividing 9. (An
  irreducible word of length `L` has `H >= 2^(L-2) + 1`, so `H = 2^s + 1` forces `L <= s + 1`: the run to length 22
  covers every one-apex word at these levels.) There are only 2 to 6 words per level and the least odd parts grow
  (21, 7, 153, 107, 1221, 221, ...). For a closing with a transitive module and `O` nonempty,
  `Pi_2 = 2 [u_2 3^t + r_2 2^t + s_2]` with `u_2 >= 1`, `s_2 <= -1`, so "odd part dividing 9" is an S-unit equation
  in 2 and 3; it has finitely many solutions `t` for each pair (CITED: the S-unit theorem, which bounds the
  solutions without vanishing subsum; the subsums `u_2 3^t + s_2 = 0` and `r_2 2^t + s_2 = 0` vanish for at most one
  `t` each; not made effective here).
- **Fans (PROVED, immediate; three special cases are `OFFSET_TRANSPORT.md` (6)).** For every tournament `X`, every
  pair with state `(r, M)` and every `t >= 0`: `r H(X) + M 2^t Pi_2(X)` is a path count. With the certificates of
  THM-5130 (E): for every odd `c <= 4097` and every `X`, both `H(X) + c p 2^i` and `c p 2^i - H(X)` are path counts
  for all large `i`, `p` the odd part of `Pi_2(X)`. The classes with `Pi_2` a power of two (all slopes) have
  `H = 1, 5, 19, 25, 217, 223, 259, 599, 695, 707, 713, 761, 773, 803, 815, 1945, 1957, 1969` (orders at most 9).
- The least `Pi_2 / H` at orders 3..9 is 4, 6.4, 7.2, 9.51, 10.67, 12.98, 14.51; among orders at most 9 the bound
  `n + 1` of Session260 / THM-5130 (H) is attained only by the point and `C_3`.

## 5. What this does and does not give for Chan-Pak

- A uniform carrier for the minus grid of `d = 31` (one kernel, two free modules). As path counts, its lines with
  `r <= 7` and the rows `c = 27, 29, 15` are in THM-5130 (E); new lines: `c = 31 2^r - 1`, `r >= 8`.
- No uniform family for the 9-grid: the module-side route is closed at orders at most 9 and for one-apex words,
  and no two-port kernel of order at most 9 is a carrier for `d = 9`.
- A reformulation: plus and minus addresses, carriers of both signs, and the signature as the invariant that tells
  them apart (HYP-9664).
- Nothing here bounds the exceptional set. Chan-Pak is OPEN.

## 6. The pasted texts, statement by statement

| statement | verdict |
|---|---|
| dimension `2^n` <-> tournaments on `n + 2` vertices; 4 classes at size 4 <-> quaternion-type algebras, 12 at size 5 <-> octonions | numerology; the exact dictionary is `2^n - 1` vertices (THM-450 (5); section 3) |
| transitive = flip (split), cyclic = spin | true for lines (Proposition 3.1) |
| loss of associativity = loss of transitivity; zero divisors at 16 | exact form: parity of cyclic triangles in zero-sum 4-sets; forced at 16 (Theorem 3.4) |
| the doubling as the canonical way to lose structure | extremal at 16 among quaternionic-line patterns, beaten at 32 (section 3.6) |
| adjacency matrices as generators: transitive -> nilpotent, circulant -> Fourier | standard linear algebra; not used |
| "self-converse ratio peaks at 5" | false as stated: the ratio is 1 for `n <= 3` and decreases after; `n = 5` is an interior local maximum (canon: `two_edge_transport_20261102/ARITHMETIC_AND_SOURCE_AUDIT.md`) |
| unit-ball volume is largest in dimension 5 | true; the link to tournaments is an analogy (canon: `arc_connectivity_route...20261104.md`) |
| "5 Ogg primes 2, 3, 5, 7, 11" | false: `X_0(p)` has genus 0 for `p = 2, 3, 5, 7, 13` (canon: `GENUS_BRIDGE.md`, MISTAKES) |
| MRB constant, "eighth residue singularity", `D = -8` | no theorem behind it; `D = -8` does occur as the discriminant of the `C_3` kernel (Session260), unrelated to MRB (canon: `PRIME-MRB.md`) |
| `A_5`, Petersen graph, `K_5`, `K_(3,3)` | canon: THM-436, THM-4529; nothing added |

## 7. Reproduce

    cd scripts
    gcc -O2 -o classgen classgen.c ; ./classgen 9 cls          # class lists cls_1.txt ... cls_9.txt (orders 1..9)
    gcc -O2 -o pairspec pairspec.c ; for m in 1..9: ./pairspec m < cls_m.txt >> pairspec_le9.txt
    python e1_summary.py pairspec_le9.txt                       # runs/e1_pair_spectrum.txt   (F1)
    gcc -O2 -o words_f2 words_f2.c ; ./words_f2 22              # runs/e1_words_f2.txt        (F2)
    gcc -O2 -o kernels kernels.c ; ./kernels m < ../data/cls_m.txt > ../runs/ker_m.txt   (m = 3..7)
    python ker_stats.py ; python ker_verify.py                  # runs/e2_*.txt               (K1, K2, identities)
    gcc -O2 -o kernels_sat kernels_sat.c ; ./kernels_sat m < codes_m | python ker_census.py  # runs/ker_census_le8.txt, ker_census_9.txt
                                                                # (order 9: 191,536 classes, about 15 minutes on 11 cores)
    python cd_signs.py ; python cd_more.py ; python cd_lines_aut.py ; python cd_orbit.py
    python cd_extremal.py 30 fano ; python cd_min.py 2 ; python cd_f5_orbits.py
    python cd_pg42.py build ../data/pg42_71_odd_planes.txt ; python cd_pg42.py compare 2000
    python cd_pg42.py verify ../data/pg42_71_odd_planes.txt

`data/cls_m.txt` (orders 1..7): canonical code and `H` of every class; code = bit `b` for the pair `(i < j)` in
lexicographic order, 1 iff `i -> j`. `runs/ker_m.txt` is kept for `m = 3..6` (the order-7 file, 1.5 MB, is not; its
summary is in `runs/e2_kernel_signatures.txt` and `runs/ker_census_le8.txt`). `opcalc.py`, `towers.py`: the
derivation check of working claim 1 (KNOWN). `grid_kernels.py` reads the order-8 table of `modfam.c` (S21 packet)
and is superseded by `ker_census.py`. `audit/`: the auditor's report, programs and outputs.

## 8. Errors of this session (MISTAKE-762)

- The pre-registration was pushed before the prior-art sweep returned: three of its four working claims and half of
  prediction C2 were already in canon.
- The pre-audit text claimed as specific `H(CD_4) = 198,335,025`, its 21 automorphisms and the double regularity of
  the sedenion tournament: all three are in my own THM-447 / THM-448 of June (`CD_4 = T_15`). It presented the minus
  grid of `d = 31` as new path counts although THM-5130 (E), audited the day before, certifies its lines with
  `r <= 7`. It missed Wilmot's theorems for all `n`, Session168's recursion and Kivunge's Hadamard statement.
- The first HYP-9664 (C), pushed with the pre-audit checkpoint, said that the doubling is extremal in every
  dimension; it is false from `n = 5` on (section 3.6). I found the 71 an hour after pushing it and the audit found
  the 64 and 512 independently.
- Proposition 2.3 said "elliptic iff `r > 0`" without the hypothesis `O` nonempty.
- A data slip caught before publication: the first run of E1 read `cls9_1.txt ... cls9_9.txt` (orders 1..9 of one
  `classgen` run) as nine pieces of order 9, which added 7,412 duplicate rows and the wrong class count 206,360. The
  tables did not change; `runs/e1_pair_spectrum.txt` is the corrected run (198,948 classes).
