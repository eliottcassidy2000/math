# Tournament Clock Prime Collatz (TCPC): clock digraphs and their combination laws; two-sheet clocks give tournaments with every arc on an odd number of Hamiltonian paths (the open cases N = 14, 34 and 38 of HYP-9167(b) are settled); the triplet principle made precise; the square-clock law of the Syracuse measure

**Session:** collatz-procgen-20260922, lane `tcpc`, 2026-10-01. Research note, not a canon theorem; the orchestrator audits and promotes. Final lane version (v3); independent audit owed.

**Owner's prompt (core):** a triplet with one member primarily differentiated and the other two differentiated by even/odd halving/doubling (hypotenuse versus legs); `F = U + S` and the 3-4-3 split of the ten proper divisors of `p^2qr`; squares mod 9 `{0,1,4,0,7,7,0,4,1}` versus cubes `{0,1,8}`; "Tournament Clock Prime Collatz": tournaments with missing, doubled or self-looped edges on partitions of a modulus, laws for how they combine, and grounding in open problems.

## Status

| # | Statement | Label |
|---|---|---|
| D | Clock digraph `Clk(G,C)`, two-sheet clock `T(S,X)`; Paley minus a vertex is the two-sheet clock of its discrete-log clock | DEFINITION + PROVED |
| L1 | Reflection law: in every clock digraph on an abelian group of odd order `n` (loops, doubled arcs, missing pairs, multiplicities allowed) every arc lies on an even number of Hamiltonian paths, and `H = n F (mod 2n)`, where `F` counts antipodal half-paths. With `M` the multiplier group, this refines to mod `2n|M|` when `M` has no involution other than `-1` (always true for tournaments), and for tournaments `F/|M|` is odd. For even order, `c(u -> v)` has the parity of the number of `rho_(u+v)`-symmetric HPs with middle arc `u -> v` | PROVED (generalizes THM-4524 C1) + FINITE-EXACT (all 5460 clock sets on `Z/m`, `m` odd `<= 13`; all 256 on `Z/3 x Z/3`; 300 random multisets; even order: every clock set on `Z/4, Z/6, Z/8, Z/10`) |
| L2 | Lexicographic law: connection sets that are unions of `K`-cosets give lexicographic products; `H(T[S_1..S_n]) = sum_walks prod_v k_v! pc_{k_v}(S_v)`; mod-2 and arc-parity corollaries (the cross-arc formula needs visit counts 1 and 2; the draft's one-term version is REFUTED) | PROVED (elementary; formula classical in spirit) + VERIFIED |
| L3 | Power-residue clocks mod `3^k` are blow-ups: squares `= C_3[3^(k-1) K_1]`, `2 3^j`-th powers `=` directed `C_(3^(j+1))` blown up, `3^j`-th powers `=` undirected `C_(3^(j+1))` blown up; squares `x` cubes `= U(3^k)`, a direct product only for `k = 2`; `x = x^3 x^4 (mod 9)` | PROVED |
| L4 | Multiplicative law `Clk(G, A B) = union_a a Clk(G,B)`; spectra combine by the twisted sum `lambda_j(AB) = sum_a lambda_(aj)(B)`; mod 9: unit clock = three Hamiltonian 9-cycles = square clock plus its reverse | PROVED |
| L5 | Complement duality `H(D) = sum_k (-1)^(n-k) k! pc_k(D-bar)`; mod 2, missing pairs and doubled pairs trade places | PROVED; mod 2 this is Berge's Stronger Theorem (CITED) |
| T1 | Vertical-arc theorem: in a two-sheet clock with reversed sheets and symmetric cross set, every vertical arc is on an odd number of HPs | PROVED |
| T2 | **All-odd tournaments on 14, 34 and 38 vertices exist** (HYP-9167(b): the first three open existence cases; `15, 35, 39` are not prime powers): `D_7` (`H = 24540117`, `Aut = Z/7`; the unique all-odd class among the 8192 two-sheet clocks on `Z/14`), `D_17`, `D_19` | FINITE-EXACT (`N = 14`: three independent engines; `N = 34`: two inclusion-exclusion engines, on `T` and on `T^op`; `N = 38`: karp4 on `T` and karp3 on `T^op` (orchestrator run), 14467258263 subset representatives each) |
| T3 | The family `D_r` (rotational sheets, distance-threshold cross arcs) is all-odd for every odd `r` from 3 to 19 (`N = 6..38`); `D_3 = QR_7 - v`, `D_5 = QR_11 - v`; `D_7, ..., D_13` are not Paley minus a vertex | FINITE-EXACT (two independent computations for every `r <= 17`, one for `r = 19`; for `r <= 13` a subset DP and an inclusion-exclusion engine agree) |
| C1 | Conjecture: `D_r` is all-odd for every odd `r`; this would prove the existence half of HYP-9167(b) for every `N = 2 (mod 4)` | OPEN (conjecture; verified `r <= 19`) |
| T4 | One-sheet control: among all circulant tournaments on `Z/m`, `m = 7, 11, 15, 19, 23`, the vertex-deleted tournament is all-odd only for Paley (none for `m = 15`); naive Paley-type products `chi_p(a) chi_q(b)` on `Z/3 x Z/5` are not doubly regular and not all-odd after deleting a vertex | FINITE-EXACT |
| G1 | Syracuse clock coordinates: `x = (-1)^c 4^s (mod 9)`; one Syracuse step is `(c,s) -> (v mod 2, 1 + c + v mod 3)` (the `+1` erases the square clock); mod `3^k` the same with the 3-adic log | PROVED (elementary) + VERIFIED on `10^6` integers |
| G2 | Doubling law: `mu(2E) = 2 mu(E)` for every `E` inside `1 + 3Z_3`, `mu` = Tao's Syracuse measure; the class `0 mod 3` is null | PROVED + FINITE-EXACT to `3^8` |
| G3 | Square-clock IFS: `3 mu` restricted to `1 + 3Z_3` is the invariant measure of `x/4, (3x+1)/4, (6x+1)/4` with weights `1/4, 1/4, 1/2`; this re-derives Tao's mod-9 law `(8,16,11,4,2,22)/63` from a 3-state clock | PROVED + FINITE-EXACT |
| TP | Triplet principle: exact core = an involution with exactly one fixed point (on `Z/3`, doubling = negation since `2 = -1 (mod 3)`) plus a halving/doubling cocycle on the swapped pair; theorems: Pythagorean legs (`c +- b` squares, `c +- a` twice squares), Syracuse classes mod 3 (G2), square-clock carries (G3), `F = S + U` as a palindromic split | PROVED for the listed instances; ANALOGY for "p^2qr is p and p^3 combined" |
| - | Any statement about the Collatz conjecture itself | none claimed |

## 0. Answers in brief

* **What TCPC is, concretely.** A *clock digraph* is a Cayley digraph `Clk(G, C)` of a finite abelian group (default `Z/m`) with a connection *multiset* `C`. It has a loop at every vertex iff `0 in C`, a doubled pair iff `c` and `-c` are both in `C`, a missing pair iff neither is, and it is a tournament iff `C` and `-C` partition `G - 0`. A *two-sheet clock* adds a sheet bit: vertex set `G x {a, b}`, with connection sets depending on the sheets. The sheet bit is exactly the owner's halving/doubling bit: Paley minus a vertex is the two-sheet clock of the discrete-log clock, the sheet being "even or odd exponent of the primitive root".
* **How clocks combine (§2).** Five exact laws: reflection (L1), lexicographic (L2), power-residue blow-ups (L3), multiplicative unions (L4), complement duality (L5). The intuitive ones are L2–L5. The surprising one is L1: the arc-parity theorem THM-4524 C1 for Cayley *tournaments* holds verbatim for *every* clock digraph of odd order, including loops, doubled arcs and missing pairs.
* **The main new result (§3).** Two-sheet clocks produce tournaments in which every arc lies on an odd number of Hamiltonian paths ("all-odd", Rédei-rigid). This includes the first three open existence cases of HYP-9167(b), `N = 14, 34, 38`, where Paley is unavailable because 15, 35 and 39 are not prime powers. The family `D_r` gives such a tournament on `2r` vertices for every tested `r` (`r <= 19`) and conjecturally for all odd `r`. Beyond `N = 26` the verification uses a new `O(1)`-memory inclusion-exclusion parity engine.
* **Collatz grounding (§5).** In clock coordinates the Syracuse step is a two-hand clock. The cube hand (`±1`, the sign mod 3) is set by the parity of `v`. The square hand is erased by the `+1` and then rotated `v` ticks. Tao's 3-adic Syracuse measure `mu` then obeys two exact laws:
  * **Doubling law (G2).** The class `2 (mod 3)` carries exactly the doubled image of the class `1 (mod 3)`, with twice the mass, and the class `0` is null. This is the owner's triplet, as a theorem.
  * **Square-clock IFS (G3).** On the class `1 (mod 3)`, `mu` is generated by one pure rotation of the square clock (`x -> x/4`) and two carries (`(3x+1)/4`, `(6x+1)/4`), with weights `1/4, 1/4, 1/2`. The two carries again differ by doubling.

  These reformulate known structure (Tao 2019; the repo's frequency recursion THM-4519/4520). They do not touch the conjecture.
* **Honest verdicts.** The triplet principle is a theorem wherever the triplet is `Z/3` with doubling (`2 = -1 mod 3`), or wherever a geometric weight `2^-v` splits by parity. It is ANALOGY elsewhere. In particular "`p^2qr` is `p` and `p^3` combined" is ANALOGY: the precise content is the exponent identity `(2,1,1) = (3,0,0) + (-1,1,1)`, repeated for every `k`-free level.

## 1. Definitions

**1.1 Clock digraphs.** `G` a finite abelian group, `C` a multiset of elements of `G` with multiplicity `mu_C`. `Clk(G, C)` has vertex set `G` and `mu_C(c)` arcs `x -> x + c` for every `x` and `c`. Each unordered pair class `{c, -c}` with `c != -c` is
* *missing* if `mu(c) = mu(-c) = 0`;
* *forward* if only `mu(c) > 0`, *backward* if only `mu(-c) > 0`;
* *doubled* if both are positive.

`0 in C` gives loops ("selfie" in the sense of THM-4524). An element with `2c = 0, c != 0` gives a digon automatically. `Clk(G, C)` is a tournament iff `C` is a set with `C ⊔ -C = G - {0}`, which forces `|G|` odd.

Natural families:
* `k`-th power residues `P_k(m) = {x^k : x in U(m)}`, as a set or as the multiset `(x^k)_(x in Z/m)` (§2.3);
* cyclotomic classes `g^i <g^e>` of a primitive root `g` (§2.3, `Z/19`);
* Paley `QR_p` (§2.1, §3);
* prime residue classes. By Dirichlet's theorem every unit class contains primes, so the "prime clock" `Clk(m, {p mod m})` is the unit clock `Clk(m, U(m))`, the unitary Cayley graph. Mod 9 it is `K_(3,3,3)` (§2.4). TCPC adds nothing beyond L2–L4 there;
* divisor-lattice classes, treated as a triplet in §4.3, not as clocks.

**TCPC dictionary (owner's words -> objects in this note).**

| owner's word | precise object |
|---|---|
| clock | a cyclic group `Z/m`, or a discrete-log clock (`F_q^* = Z/(q-1)`, `U(3^k) = Z/(2·3^(k-1))`) |
| partition of a numeric modulus | the pair classes `{c, -c}` of `C`, or the coset partition by a subgroup `K` (L2) |
| self-looped / doubled / missing edge | `0 ∈ C` / `c, -c ∈ C` / `c, -c ∉ C` |
| prime | the field case (Paley, cyclotomic clocks on `F_q`); where `q` is not a prime power, the two-sheet clock substitutes (§3) |
| how they combine | L1–L5 |
| triplet | an involution with one fixed point plus a factor-2 cocycle (§4) |
| Collatz | the Syracuse 3-adic clock (§5) |

**1.2 Two-sheet clocks.** For `r` odd, `S` a tournament set on `Z/r` (`S ⊔ -S = Z/r - 0`) and `X ⊆ Z/r` arbitrary, `T(S, X)` is the tournament on `{a_j, b_j : j in Z/r}` with
* `a_j -> a_k` iff `k - j in S` (sheet A is `Clk(r, S)`);
* `b_j -> b_k` iff `j - k in S` (sheet B is the reverse);
* `a_j -> b_k` iff `k - j in X`, and `b_k -> a_j` otherwise.

Translation `j -> j + 1` on both sheets is an automorphism. The general two-sheet clock (sheet-dependent connection sets on `Z/2r`, `x -> x + d` iff `d in C_(x mod 2)`) is the same object with independent sheets. Its tournaments number `2^((r-1)/2) 2^((r-1)/2) 2^r = 2^(2r-1)` (8192 for `r = 7`).

**Proposition 1.1 (Paley minus a vertex is a two-sheet clock; PROVED).** Let `q = 2r + 1 = 3 (mod 4)` be a prime power, `g` a primitive root of `F_q`, `chi` the quadratic character. Put `a_j = g^(2j)` and `b_j = g^(2j+1)`, so that the sheet is the parity of the discrete log. Then `QR_q - 0 = T(S_A, X)`, where
* `S_A = {d : chi(g^(2d) - 1) = 1}`;
* `X = {e : chi(g^(2e+1) - 1) = 1}`.

*Proof.* `chi(g^(2k) - g^(2j)) = chi(g^(2(k-j)) - 1)` and `chi(g^(2k+1) - g^(2j+1)) = -chi(g^(2(k-j)) - 1)`. Also `chi(g^(-2d) - 1) = chi(-1) chi(g^(2d) - 1) = -chi(g^(2d) - 1)`, so `S_A` is a tournament set and sheet B is its reverse. Finally `chi(g^(2k+1) - g^(2j)) = chi(g^(2(k-j)+1) - 1)`. ∎

Isomorphism checked with nauty for `q = 7, 11, 19, 23`. So the discrete-log clock of `F_q^*`, cut into its two square classes, *is* the Paley-minus-a-vertex tournament. The owner's "doubled" clock is the sheet bit.

**1.3 Invariants.** For a digraph with arc multiplicities:
* `H` = Hamiltonian paths (HPs);
* `HC` = directed Hamiltonian cycles;
* `c(e)` = HPs through the arc `e`;
* `start(v)`, `end(v)`;
* "all-odd" = every `c(e)` odd (Rédei-rigid in THM-4524's language);
* `F` = antipodal half-paths (L1);
* the DFT spectrum `lambda_j = sum_(c in C) zeta^(jc)`;
* `|Aut|` (nauty `dreadnaut`).

Three independent HP engines are used:
* the C program `procgen_tcpc_20261001_hp.c` (128-bit forward×backward convolution for `n <= 18`, and a bitset parity mode, used up to `n = 22`);
* a pure-Python subset DP;
* a Python memoized recursion computing `c(e) = H(T) - H(T - e)`.

For `n = 26` a fourth, memory-lean C engine `procgen_tcpc_20261001_par26.c` computes `H(T - e) mod 2` by a forward-only DP. It was validated against the full engine on every orbit representative for `r <= 11`.

Beyond `N = 26` the subset DPs need too much memory. Two inclusion-exclusion engines, `procgen_tcpc_20261001_karp3.c` and `_karp4.c`, use `O(1)` memory and an algorithm independent of the DPs. They work for two-sheet clocks with translation `tau`, where every `tau`-orbit `O_t` of arcs has odd size `r`.

*The identity.* Over `GF(2)[e_1, ..., e_K]/(e_i e_j)`,

```text
sum_(S ⊆ V)  sum_(walks with N vertices inside S)  prod_(arcs) (1 + sum_t e_t [arc in O_t])  =  H + sum_t e_t r c(e_t)   (mod 2).
```

Inclusion-exclusion over `S` kills every walk that repeats a vertex, and its signs vanish mod 2. The sum over `S` can be restricted to one subset per `tau`-orbit, because orbit sizes divide `r`, which is odd.

*How the walks are counted.* For each subset the engines build forward and backward walk-parity vectors.
* karp3 extracts the marked-arc counts orbit by orbit.
* karp4 extracts all of them at once, as cyclic cross-correlations computed by carry-less multiplication (ARM64 `PMULL`).

*Validation.* karp3 and karp4 agree with each other and with the subset DPs on `D_r`, on the interval controls and on random two-sheet clocks with mixed parities, for every `r <= 13`.

## 2. Combination laws

### 2.1 L1, the reflection law (PROVED)

**Theorem L1.** Let `D = Clk(G, C)` with `|G| = n` odd and `C` any multiset (loops, doubled arcs and missing pairs allowed). Let `M = {u in Aut(G) : u C = C}` be the multiplier group, acting by `x -> u x`.
1. Every arc of `D` lies on an even number of Hamiltonian paths.
2. `H(D) = n F(D) (mod 2n)`, where `F(D)` is the number of directed paths `0 -> w_1 -> ... -> w_((n-1)/2)` of `D` whose vertices contain exactly one element of every antipodal pair `{x, -x}`, `x != 0` (weighted by multiplicities). Always `|M|` divides `F(D)`. If `M` has no involution other than possibly `-1`, then `H(D) = n F(D) (mod 2n|M|)`. This holds, for example, for `G` cyclic of prime-power order, and whenever `|M|` is odd.
3. If `D` is a tournament, then `|M|` is odd and `F(D)/|M|` is odd.

*Proof.* For every `g in G` the reflection `rho_g(x) = g - x` maps the arc `x -> x + c` to the arc `g - x -> g - x - c`. That is an arc of the reverse digraph with the same label `c`, so `rho_g` is an isomorphism `D -> D^op` for **every** connection multiset (not only for tournaments). Hence `P -> reverse(rho_g(P))` maps HPs to HPs. Together with translations and multipliers this gives an action of `G ⋊ (M x {±1})` on the HPs.
* **Free part.** Translations and multipliers act freely: an automorphism fixing a path as a sequence fixes every vertex.
* **Reflection fixed points.** A path fixed by `rho_g` has middle vertex `g/2`, since `n` is odd.
* **At most one reflection.** No path is fixed by two reflections, because their product is a non-trivial translation.
* **Conjugacy.** All reflections are conjugate (`rho_g = tau_(g/2) rho_0 tau_(-g/2)`), so each fixes `F` paths.
* **Item 2.** An element `x -> g - ux` (with path reversal) can fix a path only if its square is the identity, i.e. `u^2 = 1` and `g(1-u) = 0`.
  * `u = -1` (possible when `C = -C`) acts as a translation with reversal, which fixes no path.
  * If `M` has no other involution, only the reflections (`u = 1`) have fixed paths. Then orbits have size `2n|M|` (free) or `n|M|` (reflection-stabilised), and the `nF` reflection-fixed paths form `F/|M|` orbits. Hence `H = 2n|M| (#free orbits) + n F`.
  * The weaker congruence mod `2n` uses translations and reflections only, so it always holds.
  * `|M|` divides `F` because `M` commutes with `rho_0` and permutes the `rho_0`-fixed paths freely.
* **Item 1.** For the arc `e = u -> v`, the reflection `rho_(u+v)` maps `e` to itself. A path through `e` fixed by `rho_(u+v)` would need `e` in its middle position, which does not exist for `n` odd. So the involution has no fixed points on the HPs through `e`.
* **Item 3.** `H` is odd (Rédei) and `n` is odd. Let `u != 1` be an involution of `G`. Since `|G|` is odd, `G = ker(u-1) ⊕ ker(u+1)`, so some `x != 0` has `ux = -x`, and a tournament set cannot be `u`-invariant. So `M` has no involution, and `|M|` is odd. ∎

*Even order (PROVED).* If `|G|` is even, the same involution `P -> reverse(rho_(u+v)(P))` acts on the HPs through `e = u -> v`. Its fixed points are the `rho_(u+v)`-symmetric HPs with `e` as the middle arc. So `c(e) = M(e) (mod 2)`, where `M(e)` counts the half-paths ending at `u` (built backwards) that, together with their mirror images, partition `G`. If `rho_(u+v)` fixes a vertex, `M(e) = 0`. VERIFIED on every arc out of 0 of every clock set on `Z/4, Z/6, Z/8, Z/10` (2844 arcs).

*Remarks.*
* Item 1 for Cayley tournaments is THM-4524 C1; the proof is the same, and the content of L1 is that the tournament hypothesis is never used.
* A `rho_0`-fixed path is determined by a "half-path" from `0` picking one of each `{x, -x}`. This is the antipodal triplet `{-x, 0, x}` of §4: `0` fixed, `±x` swapped.

**Data (FINITE-EXACT).** L1 holds for all 5460 sets `C ⊆ Z/m - 0` (`m = 3, 5, ..., 13`), and for 300 random multisets with loops (`m <= 11`).

Paley `F` values:

| `p` | `F` | `F/((p-1)/2)` |
|---|---|---|
| 3 | 1 | 1 |
| 7 | 9 | 3 |
| 11 | 185 | 37 |
| 19 | 573057 | 63673 |
| 23 | 63871533 | 5806503 |

All quotients are odd, as item 3 requires. Example: `H(QR_7) - 7F = 189 - 63 = 3 (2 x 7 x 3)`.

### 2.2 L2, lexicographic decomposition and the composition formula (PROVED)

**L2a (coset unions are lexicographic products).** Let `K <= G` and let `mu_C` be constant on every `K`-coset outside `K`. Then
`Clk(G, C) = Clk(G/K, C-bar)[Clk(K, C ∩ K)]`,
where `C-bar` is the image multiset of `C - K`, with the common multiplicity per coset. Here `D[E]` replaces each vertex of `D` by a copy of `E`, and arcs between copies follow `D`.

*Proof.* Arc multiplicities between different cosets depend only on the coset difference; inside a coset they form `Clk(K, C ∩ K)`, identified by any section. ∎

**L2b (composition formula).** For digraphs `T` on `[n]` and `S_1, ..., S_n`,

```text
H(T[S_1, ..., S_n]) = sum over block walks w of T  prod_v  k_v(w)! pc_(k_v(w))(S_v)
```

Here:
* a block walk is a sequence of vertices of `T` with consecutive entries joined by arcs, weighted by the product of their multiplicities, and visiting every `v` exactly `k_v >= 1` times;
* `pc_k(S)` is the number of covers of `V(S)` by `k` vertex-disjoint directed paths.

*Proof.* Cut a Hamiltonian path of `T[S]` into maximal runs inside one block. The runs of block `v` form an ordered cover of `S_v` by `k_v` paths, and the block sequence is a walk of `T`. This decomposition is reversible. ∎

**Corollaries (PROVED).**
1. **Mod 2.** `H(T[S_1..S_n]) = H(T) prod H(S_v) (mod 2)`, since `k! = 0 (mod 2)` for `k >= 2`.
2. **Inner arcs.** For an arc `e` inside block `v`, `c(e) = H(T) c_(S_v)(e) prod_(u != v) H(S_u) (mod 2)`: the `k_v!` symmetry survives the constraint "uses `e`".
3. **Cross arcs.** For an arc `(u,x) -> (v,y)` between blocks the visit symmetry is broken at `u` and `v`. Only `(k_u - 1)! (k_v - 1)!` survives there, so visit counts 1 *and 2* contribute mod 2:

   `c = prod_(w != u,v) H(S_w) · sum_(a,b in {1,2}) N_ab(u -> v) pc_a(S_u; ends at x) pc_b(S_v; starts at y) (mod 2)`,

   where `N_ab` counts the walks of `T` that visit `u` `a` times, `v` `b` times and every other vertex once, with a marked step `u -> v`.
4. **Mod 3.** Only `k_v in {1, 2}` survive (in H; for cross arcs `k_u, k_v <= 3`).

*Corrected near miss (REFUTED by the runner).* The draft kept only the `(a, b) = (1, 1)` term, `c = c_T(u->v) end_(S)(x) start_(S)(y) (mod 2)`. The runner's random lexicographic products refute it (section C of the output), and the four-term formula above holds on every tested cross arc. The draft's two consequences are withdrawn:
* "all-even is closed under lexicographic products" survives only where it follows from L1: `T[S]` is all-even when `T` and `S` are Cayley tournaments of odd abelian groups, because `T[S]` is then a Cayley tournament of `G_T x G_S`.
* "no lexicographic all-odd tournament below 13" is replaced by a FINITE-EXACT statement: no `T[TT_2]` and no `TT_2[T]` with `|T| <= 7`, and no `T[QR_7 - v]` with `|T| <= 4`, is all-odd. Two-sheet clocks are not lexicographic products: they reverse one sheet.

**Mod-9 instances (FINITE-EXACT, L2b reproduces every value).**

| `C ⊆ Z/9` | structure | `H` | `HC` |
|---|---|---|---|
| `Sq = {1,4,7}` (unit squares) | `C_3[3K_1]` (missing inner pairs) | 648 | 72 |
| `{2,5,8}` (non-squares) | `C_3^op[3K_1]` | 648 | 72 |
| `Sq ∪ {3}` | `C_3[C_3]`, tournament | 3159 | 207 |
| `Sq ∪ {3,6}` | `C_3[K_3*]` (doubled inner pairs) | 14256 | 1152 |
| `Cu = {1,8} = {±1}` | undirected `C_9` | 18 | 2 |
| `U(9) = Sq·Cu` | `K_3*[3K_1]` = doubled `K_(3,3,3)` | 37584 | 3168 |
| `Sq ∪ Cu = {1,4,7,8}` | mixed (doubled `±1`, one-way `4, 7`) | 2268 | 154 |
| `{1,2,3,4}` | rotational `R_9` | 3267 | 222 |

The formula is explicit for these rows.
* **Lexicographic products over `C_3`.** `H(C_3[S]) = 3 sum_(q>=1) sum_(rho=0..2) P_(q+1)^rho P_q^(3-rho)`, with `P_k = k! pc_k(S)`. It gives:
  * `648` for `S = 3K_1` (only `P_3 = 6`);
  * `3159` for `S = C_3` (`P = 3, 6, 6`);
  * `14256` for `S = K_3*` (`P = 6, 12, 6`).
* **Smirnov-word counts.** `H(K_3*[3K_1]) = 6^3 x 174`, where 174 is the number of words with three letters each three times and no two equal neighbours.

So the inner "multiples of 3" clock modulates `H` through the path-cover polynomial of the inner clock alone.

### 2.3 L3, power-residue clocks mod `3^k` (PROVED)

`U(3^k) = {±1} x (1 + 3Z/3^k) ≅ Z/2 x Z/3^(k-1)`, cyclic of order `2·3^(k-1)` and generated by 2. Let `d = 2^a 3^j` with `j <= k-1`.
* **The `d`-th powers.** If `a >= 1` they are `1 + 3^(j+1) Z/3^k`; if `a = 0` they are `±1 + 3^(j+1) Z/3^k`.
* **As clocks (by L2a).**
  * `Clk(3^k, P_d) = vec-C_(3^(j+1)) [3^(k-j-1) K_1]` for `a >= 1`: a directed cycle, blown up.
  * `Clk(3^k, P_d) = C_(3^(j+1)) [3^(k-j-1) K_1]` for `a = 0`: an undirected cycle, blown up.
  * So squares give `C_3[3^(k-1) K_1]` and cubes give `C_9[3^(k-2) K_1]`.
* **Fractal reading.** Squares see only the residue mod 3, cubes only mod 9, `3^j`-th powers only mod `3^(j+1)`. The cube map zooms the 1-unit clock out by one 3-adic level: `(1 + 3x)^3 ∈ 1 + 9Z`. As a group, the image of cubing in `U(3^(k+1))` is isomorphic to `U(3^k)`.
* **Squares times cubes.** `Sq·Cu = U(3^k)` for all `k >= 2`, but `Sq ∩ Cu = 1 + 9Z/3^k` has order `3^(k-2)`. So the product is direct (`Z/6 = Z/3 x Z/2`) only for `k = 2`.
* **Mod 9.** The CRT idempotents of `Z/6` are `3` and `4`, so `x = x^3 · x^4` with `x^3 ∈ {±1}` (the cube, or sign, part) and `x^4 ∈ {1,4,7}` (the square part). For example `2 = 8 · 7`. The owner's square sequence `{0,1,4,0,7,7,0,4,1}` is the multiset `(x^2)_(x ∈ Z/9)` = `{0^3, 1^2, 4^2, 7^2}`. As a clock it is `2·Clk(9, Sq)` plus 3 loops per vertex, so its HP count is `2^8 · 648`. The cube sequence `(x^3) = {0^3, 1^3, 8^3}` is `3·C_9` plus loops.

**Field versus ring: mod 7 and mod 9 (PROVED).**
* **The same clock group.** `U(7)` and `U(9)` are both `Z/6 = Z/3 x Z/2` = squares x cubes:
  * mod 7: squares `{1,2,4}`, cubes `{1,6}`;
  * mod 9: squares `{1,4,7}`, cubes `{1,8}`.
* **Field.** On `Z/7` every nonzero residue is a unit, so the square clock is a tournament: the Paley heptagon `QR_7 = {1,2,4} = <2>`, the trivial-cycle code of the Paley-bridge note.
* **Ring.** On `Z/9` the non-units `3, 6` are the missing pairs, so the square clock is the blow-up `C_3[3K_1]`. Here `<2>` is the whole unit group and the squares are `<4>`.

So the Collatz-relevant triplets `{1,2,4} mod 7` (trivial cycle) and `{1,4,7} mod 9` (the coset of `3A + 1`) are the same `Z/3` square clock, once on a field and once on a ring.

**Cyclotomic unions: which combination wins (FINITE-EXACT).** On `Z/19`, the 8 tournaments that are unions of sextic cyclotomic classes (one class from each pair `{C_i, C_(i+3)}`, since `-1 ∈ C_3`) fall into two isomorphism classes:

| class | unions | `H` | `|Aut|` | minus a vertex |
|---|---|---|---|---|
| Paley | `C_0 ∪ C_2 ∪ C_4` and its reverse | 1172695746915 | 171 | all-odd |
| mixed | the other 6 unions | 1167595581285 | 57 | 117/153 arcs odd |

The alternating union (the squares) maximises `H` and is the only one whose vertex-deleted tournament is Rédei-rigid.

### 2.4 L4, multiplicative combination (PROVED)

For `A` a set of units and `B` any multiset, `Clk(G, A·B) = union_(a ∈ A) a·Clk(G, B)` as arc multisets. Here `a·D` is the image of `D` under `x -> ax`, and `A·B` is the product multiset. The DFT eigenvalues combine by the twisted sum `lambda_j(A·B) = sum_(a ∈ A) lambda_(aj)(B)`. These are Gauss-period identities.

Spectra mod 9 (`zeta = e(1/9)`, `omega = zeta^3`):
* `lambda_j(Sq) = 3 zeta^j` for `3 | j` and 0 otherwise (`3·spec(C_3)` plus six zeros: the blow-up);
* `lambda_j(Cu) = 2 cos(2 pi j/9)`;
* `lambda_j(U(9))` is the Ramanujan sum: `6, -3, -3` at `j = 0, 3, 6` and 0 otherwise.

Mod 9, with `Sq·Cu = U(9)`, the multiplicative law reads in both orders:
* `K_(3,3,3) = C_9 ∪ 4C_9 ∪ 7C_9`: three Hamiltonian 9-cycles with steps `±1, ±4, ±2`;
* `K_(3,3,3) = C_3[3K_1] ∪ (-1)·C_3[3K_1]`: the square clock and its reverse.

So the cube factor `{±1}` acts on clocks by doubling every arc (arc reversal), and the square factor `{1,4,7}` rotates the three blocks. That is the precise form of "`Z/6 = Z/3` (squares) `x Z/2` (cubes, `±1`)".

### 2.5 L5, complement duality (PROVED; mod 2 CITED as Berge's Stronger Theorem)

For any digraph `D` on `n` vertices with complement `D-bar` (in the complete loopless digraph),

```text
H(D) = sum_(k=1..n) (-1)^(n-k) k! pc_k(D-bar),   hence   H(D) = H(D-bar) (mod 2).
```

*Proof.* `H(D) = sum_(pi) prod_i (1 - y_i)`, where `y_i` indicates that the `i`-th step of `pi` is an arc of `D-bar`. Expand, and group the factors by the runs they create. ∎

The complement turns missing pairs into doubled pairs and back, and reverses one-way arcs. So mod 2, "missing" and "doubled" are dual. For tournaments this gives nothing new (`D-bar = T^op`). It is the parity backbone used in §3.

In the language of mixed graphs, where a doubled pair is a non-oriented edge and a missing pair is a non-edge, the mod-2 statement is Berge's Stronger Theorem. Rédei's own Stronger Theorem (1934) gives `c_T(e) = c_(T')(e') (mod 2)` when the arc `e` is reversed to `e'`. Both are CITED from Schweser, Stiebitz and Toft, "The Tournament Theorem of Rédei revisited", arXiv:2510.10659 (2025), Theorems 1.1–1.3. The integer identity is the inclusion-exclusion behind them. No literature search found the all-odd property of §3 studied; priority is not claimed.

### 2.6 CRT and tensor products (classical; recorded for completeness)

* `Clk(m_1, C_1) ⊗ Clk(m_2, C_2) = Clk(m_1 m_2, CRT(C_1 x C_2))` when `gcd(m_1, m_2) = 1`, with spectra multiplying.
* The directed cycles satisfy `C_(m_1) ⊗ C_(m_2) ≅ gcd(m_1, m_2) · C_(lcm)`: two clocks tick as one exactly when their periods are coprime.

This is the prime-gear model of HYP-2081 in digraph form.

## 3. Two-sheet clocks and Rédei-rigid tournaments

THM-4524 calls a tournament *all-odd* if every arc lies on an odd number of HPs. All-odd needs `N = 1, 2 (mod 4)`. It exists at `N = 6` (unique: `QR_7 - v`) and `N = 10` (two classes), and does not exist at `N = 3, 4, 5, 9`. HYP-9167 conjectures:
* (a) `QR_q - v` is all-odd for every prime power `q = 3 (mod 4)`;
* (b) all-odd tournaments exist iff `N = 2 (mod 4)`.

Its first open cases were `N = 13` (non-existence) and `N = 14` (existence). Proposition 1.1 shows `QR_q - v` is a two-sheet clock, so the natural TCPC move is to search two-sheet clocks where no field exists.

### 3.1 The vertical-arc theorem (PROVED)

**Theorem T1.** Let `T = T(S, X)` with `X = -X`. Then every vertical arc (between `a_j` and `b_j`) lies on an odd number of HPs.

*Proof.*
* **An anti-automorphism.** The sheet swap `psi: a_j <-> b_j` maps `a_j -> a_k` to `b_j <- b_k`, and `b_j -> b_k` to `a_j <- a_k`. It maps `a_j -> b_k` (`k - j in X`) to the pair `{b_j, a_k}`, where `a_k -> b_j` because `j - k in X`. So `psi` is an anti-automorphism.
* **The group action.** `P -> reverse(psi(P))` is a bijection of HPs. With translations, the group `<tau, psi> ≅ Z/2r` acts on arcs, and `c` is invariant: `c(u -> v) = c(psi v -> psi u)`.
* **Orbit sizes.** `psi` fixes every vertical arc. So the vertical arcs form one orbit of odd size `r`, and every other orbit has size `2r`.
* **Parity.** `sum_e c(e) = (2r - 1) H` is odd by Rédei. Hence `c(vertical) = 1 (mod 2)`. ∎

This is the two-sheet analogue of the antipodal half of THM-4524 C2. It was verified on all 5688 vertical arcs of all symmetric two-sheet clocks with `r <= 9`.

### 3.2 Exhaustive two-sheet census (FINITE-EXACT)

| `N` | two-sheet clock tournaments | all-odd (labelled) | all-odd classes | the class |
|---|---|---|---|---|
| 6 | 32 | 12 | 1 | `QR_7 - v`, `H = 45`, `|Aut| = 3` |
| 10 | 512 | 40 | 1 | `QR_11 - v`, `H = 15745`, `|Aut| = 5` |
| 14 | 8192 | 84 | **1** | **`D_7`, `H = 24540117`, `|Aut| = 7`** |

The labelled counts equal `(r - 1)·N`. For `r` prime, that is the number of two-sheet labellings of a tournament with `Aut = Z/r`: choose a generator of `Aut` (`r - 1` ways), a base vertex (`N` ways) and a base vertex on the other orbit (`r` ways), then divide by `|Aut| = r`.

**One-sheet control (FINITE-EXACT).** Circulant tournaments are the one-sheet clocks. For every circulant tournament on `Z/m` with `m = 7, 11, 15, 19, 23` (2, 4, 16, 30, 94 classes up to multipliers), the vertex-deleted tournament is all-odd exactly when the circulant is Paley; for `m = 15` it never is. The naive Paley-type product `chi_3(a) chi_5(b)` on `Z/3 x Z/5`, completed on `0 x Z/5` by either tournament set, is not doubly regular. Its `H` is `197728485` or `197094945`, and after deleting a vertex only 73 or 41 of the 91 arcs are odd. So within one-sheet clocks only the field works. The second sheet is what replaces the missing field at `N = 14`.

**Theorem T2 (FINITE-EXACT).** There are tournaments on 14, on 34 and on 38 vertices in which every arc lies on an odd number of Hamiltonian paths. This settles the first three open existence cases of HYP-9167(b). An explicit one is `D_7`. In two-sheet clock form on `Z/14` (`x -> x + d` iff `d ∈ C_(x mod 2)`):
* `C_0 = {2, 4, 6, 7, 9, 11}`;
* `C_1 = {1, 8, 9, 10, 11, 12, 13}`.

Its numbers:
* `H = 24540117` and `HC = 1001369`;
* the 91 arc counts take 7 values, all odd, between `3085307` and `4087295`;
* `start = 1641842 / 1863889` and `end = 1863889 / 1641842` on the two sheets.

All three engines agree on all 91 arcs.

The same holds on 34 and on 38 vertices: `D_17` and `D_19` (§3.3) are all-odd, which settles the next two open existence cases (`35 = 5·7`, `39 = 3·13`).
* `D_17` was verified by the inclusion-exclusion engines karp3 on `T` and karp4 on `T^op`: all 33 arc orbits odd, `H` odd.
* `D_19` was verified by karp4: all 37 arc orbits odd, `H` odd, 14467258263 subset representatives, 4352 s. This is a single computation. A confirming run on `T^op` was started and then stopped at wrap-up; the auditor should repeat it.

### 3.3 The family `D_r` (FINITE-EXACT for `r <= 19`)

Let `m = (r-1)/2` and `S = {1, ..., m}`, so that both sheets are the rotational tournament `R_r` and its reverse. Let
* `Y = {0, ±1, ..., ±(m-1)/2}` if `m` is odd (`r = 3 (mod 4)`);
* `Y = {±1, ..., ±m/2}` if `m` is even (`r = 1 (mod 4)`).

Put `D_r = T(S, Y)`. In words: `a_j -> b_k` iff the circular distance `|k - j|` is at most `floor(m/2)`, with distance 0 included exactly when `r = 3 (mod 4)`. So the orientation of the vertical arc is set by `r mod 4`.

**Quadrant form (PROVED; a restatement).** Let `omega = e^(2 pi i/r)` and `z = omega^(k-j)`. Then:
* `a_j -> a_k` iff `Im z > 0`;
* `b_j -> b_k` iff `Im z < 0`;
* for `k != j`, `a_j -> b_k` iff `Re z > 0`;
* the vertical arc goes `a_j -> b_j` iff `r = 3 (mod 4)`.

This holds because `r/4` is never an integer, so `|k - j| <= floor((r-1)/4)` iff `cos(2 pi (k-j)/r) > 0`.

Compare Paley minus a vertex in the discrete-log coordinates of Proposition 1.1. There the three relations are `chi(z - 1) = +1`, `chi(z - 1) = -1` and `chi(z - 1) = +1`, with `z` a square, a square and a non-square respectively. So the field's quadratic character of `z - 1` is replaced by the quadrant of `z` on the circle. That is the precise sense in which `D_r` is "Paley minus a vertex without a field".

**Circle form (PROVED; a restatement).** Put `x_j = omega^j` and `y_j = i omega^j`: the `r`-th roots of unity and their quarter turn, `2r` distinct points with no antipodal pair. Let `Circ` be the circular tournament on them: `u -> v` iff `v` lies in the open counterclockwise half-circle from `u`, i.e. `Im(v/u) > 0`. Then `D_r` is obtained from `Circ` (with `a_j = x_j`, `b_j = y_j`) by
* reversing every arc among the `y_j`;
* and, when `r = 1 (mod 4)`, also reversing the `r` vertical arcs `x_j -> y_j`.

The reason is that `Im(i z) = Re z`. In the owner's words, the second clock is the first one turned a quarter and run backwards. (Checked for every odd `r < 60` in scratch.)

**Symmetries (PROVED; checked `r < 30`).**
* The translations `j -> j + 1` form `Z/r <= Aut(D_r)`, and nauty gives equality for `r <= 11`.
* The sheet swap `psi` is an anti-automorphism, so `D_r` is self-converse.
* `chi: a_j -> b_(-j), b_j -> a_(-j)` is an isomorphism from `D_r` onto its switching with respect to sheet A, which reverses all arcs between the sheets.

So `D_r` is a fixed point, up to isomorphism, of the loop gauge of THM-4524 applied to a whole sheet.

| `r` | `N` | all-odd | identity | evidence |
|---|---|---|---|---|
| 3 | 6 | yes | `QR_7 - v` | three subset engines; karp3; karp4 |
| 5 | 10 | yes | `QR_11 - v` | three subset engines; karp3; karp4 |
| 7 | 14 | yes | new (T2) | three subset engines; karp3; karp4 |
| 9 | 18 | yes | new: not `QR_19 - v` (nauty); `H = 116670839805` | exact 128-bit engine, full parity engine, lean engine; karp3; karp4 |
| 11 | 22 | yes | new: not `QR_23 - v` | full parity engine, lean engine; karp3; karp4 |
| 13 | 26 | yes | new: not `QR_27 - v` (`|Aut| = 13` vs `39`) | lean engine on `T` and on `T^op`; karp3; karp4 |
| 15 | 30 | yes | (`31` is prime; `D_15` is the interval clock) | karp3 and karp4 |
| 17 | 34 | **yes** | **new existence case: `35 = 5·7` is not a prime power** | karp3 on `T` (1010580543 subset representatives, 802 s) and karp4 on `T^op` |
| 19 | 38 | **yes** | **new existence case: `39 = 3·13` is not a prime power** | karp4 on `T` (14467258263 subset representatives, 4352 s) |

Controls:
* the plain interval `X = {0..m-1}` is all-odd exactly for `r = 3 (mod 4)`; it fails at `r = 5, 9, 13`;
* with the interval, `r = 9` has 135 of 153 arcs odd, and `r = 13` has 23 of 25 arc orbits odd.

**Census of `X` for rotational sheets (translation classes).**

| `r` | all-odd translation classes of `X` |
|---|---|
| 5 | 2 |
| 7 | 2 |
| 9 | 2 |
| 11 | 6 |
| 13 | 2 (of 632; scratch run with karp4, not in the runner) |

The classes come in pairs under the converse-with-sheet-swap involution `X -> -(Z/r - X)`. At `r = 5, 7, 9, 13` the only class is `D_r`. At `r = 11` the six `X`-classes form three isomorphism classes (nauty):
* `D_11`, from `{0..4}` and `{0..5}`;
* `{0,1,4,5,9}` together with `{0,1,2,4,7,8}`;
* `{0,1,3,7,8}` together with `{0,1,2,5,6,9}`.

None of the last four sets is centrally symmetric (checked), and neither class is `QR_23 - v`.

**Other sheets at `r = 9` (FINITE-EXACT).** Over all 16 tournament sheets `S` and all `X` (176 orbit representatives), exactly two all-odd classes occur: `D_9` and `QR_19 - v` (sheet `{3,5,7,8}`). So at `N = 18` there are at least two all-odd classes, and at `N = 22` at least four (`QR_23 - v`, `D_11`, and the two further rotational classes, each with `|Aut| = 11`).

**The apex completion.** Adding a vertex `∞` with `a_j -> ∞ -> b_j` turns `D_r` into a regular tournament `Q_r` on `2r + 1` vertices.
* `Q_3 = QR_7` and `Q_5 = QR_11`.
* For `r = 7, 9, 11`, `|Aut(Q_r)| = r`, and `Q_r` is not doubly regular.

So `D_r` is "Paley minus a vertex" with the field replaced by the rotational clock. For `r = 3, 5` the two constructions coincide.

*Observations on `Q_r` (FINITE-EXACT, not explained):*
* `Q_r` is all-even for `r = 3, 5, 7, 9`, although `Q_7, Q_9` are not vertex-transitive. Only the vertical arcs are forced even, by the anti-automorphism `psi` extended with `∞ -> ∞`.
* `Q_7 - v` is all-odd for every vertex `v` (`∞`, `a_0`, `b_0`), but `Q_9 - a_0` is not (97 of 153 arcs odd).
* "`Q` all-even" and "`Q - v` all-odd" are not equivalent in general. Over all `(Q, v)` with `|Q| = 7` (456 classes) the four combinations occur 7 / 189 / 23 / 2973 times (both / even only / odd only / neither).

**Conjecture C1 (OPEN).** `D_r` is all-odd for every odd `r >= 3`. Consequently all-odd tournaments exist for every `N = 2 (mod 4)`, which is the existence half of HYP-9167(b).
* **Evidence.** All nine cases `r <= 19`, including the non-prime-power sizes `N = 14, 34, 38`; the proved vertical-arc part T1.
* **What is missing.** A mechanism for the non-vertical arcs. The `Z/2r` action has orbits of size `2r` there, so the counting argument of T1 gives nothing.
* **Not affine.** The arc parities are not affine in `X` over `GF(2)` (tested at `r = 3, 5`), so no linear-algebra shortcut is available.

**What this means for HYP-9167.**
* (b) existence is settled at `N = 14, 34, 38`. With Paley, the remaining open existence cases below 60 are `N = 50, 54` (`51, 55` are not prime powers). C1 covers them, but they are beyond computation: `N = 50` needs about `4.5·10^13` subset representatives.
* (a) is untouched. `D_9, D_11, D_13` are different from `QR_19 - v, QR_23 - v`, and from `QR_27 - v`, which the selfie lane verified.
* **Non-existence at `N = 13` (EMPIRICAL, negative).** Structured families were searched: tournaments with a `Z/11` action (orbits 11+1+1; 256 candidates), a `Z/7` action (7 + six fixed vertices; 28672), and a `Z/5` action (5+5+1+1+1; 262144). None is all-odd. The best has 72 of 78 arcs odd (in the `Z/5` family). This is consistent with (b) but proves nothing: circulants on `Z/13` are all-even by L1, and there is no odd-order analogue of the two-sheet clock (`Z/6` carries no tournament).

## 4. The triplet principle

**4.1 The exact core.** A "triplet with one member primarily differentiated and the other two subtly differentiated by halving/doubling" is, precisely:
* a 3-element structure `{z, x, y}`;
* an involution fixing `z` and swapping `x <-> y` (a transposition in `S_3`: one fixed point, one 2-cycle);
* a `Z/2`-valued cocycle telling `x` from `y`, given by a factor 2.

The universal instance is `(Z/3, x2)`. Since `2 = -1 (mod 3)`, doubling is negation: it fixes `0` and swaps `1 <-> 2`. `p = 3` is the only prime with `2 = -1`, so it is the only modulus where halving/doubling *is* the sign involution. Every exact instance below is this `Z/3` or a geometric weight `2^-v` split by parity.

**4.2 Instances that are theorems.**

| triplet | differentiated member | the pair and its halving/doubling | status |
|---|---|---|---|
| primitive Pythagorean `(a, b, c)`, `b` even | `c` (odd; `c^2` is the sum) | `c ± b = (m ± n)^2` are squares; `c ± a = 2m^2, 2n^2` are twice squares. The legs differ by the square class of 2 | PROVED (classical parametrisation) |
| `Z/3` under `x2` | `0` | `1 <-> 2`, doubling = negation | PROVED |
| unit squares mod 9 `{1,4,7}` under squaring | `1` (fixed) | `4 <-> 7` (`2^2 <-> 2^4`, exponent doubling mod 6) | PROVED |
| Syracuse measure on `Z/3` classes (G2) | class `0` (null: no Syracuse value is a multiple of 3) | class 2 = doubled image of class 1, with doubled mass `2/3 : 1/3` | PROVED (§5.3) |
| square-clock IFS (G3) | the rotation `x/4` (isometry, no carry) | the carries `(3x+1)/4` and `(6x+1)/4 = (3(2x)+1)/4`: argument doubled and weight doubled (`1/4` vs `1/2`) | PROVED (§5.4) |
| `F = S + U` | `F` (the total) | `S = 2^U - 1` for every non-squarefree composite (squarefree divisors are subsets of the primes, "doubling per prime"), `U = r` linear; `F = S + U` iff `tau(N) = 2^U + U + 1` (DB2) | PROVED (inherited DB1/DB2) |
| antipodal pairs in a clock of odd order (L1) | `0` | `±x` swapped by the reflection; the half-path picks one of each | PROVED |
| two-sheet clocks (§3) | the vertical arcs (fixed by `psi`) | the other arcs come in `psi`-pairs | PROVED (T1) |

**4.3 The divisor triplet.** Write `N = p^2 q r = s^2 a b`. Its 10 proper nontrivial divisors fall into three classes:
* primes `{p, q, r}`, which are in `S` and in `U`;
* squarefree composites `{pq, pr, qr, pqr}`, in `S` only;
* square-containing divisors `{p^2, p^2q, p^2r}`.

That is the owner's 3-4-3 split.

*Precise statements.*
* **Palindrome (tautology).** For every `N >= 2`, `F = S + U` iff the split `(U, S - U, F - S)` is a palindrome, i.e. #square-containing = #primes. The three solutions give `0-0-0` (`p`), `1-0-1` (`p^3`) and `3-4-3` (`p^2qr`).
* **The matching (inherited from the divisor-balance note §3.1).** For `p^2qr` the square-containing divisors are paid by the primes through `d -> p·lcm(p, d)`: `p -> p^2`, `q -> p^2q`, `r -> p^2r`.
* **A second 3-4-3 (numerical coincidence).** The rank (`Omega`) distribution of the proper divisors of `p^2qr` is also 3-4-3: `(1+t+t^2)(1+t)^2 = 1 + 3t + 4t^2 + 3t^3 + t^4`. But it is a different partition: `p^2` and `pqr` trade places. Only the rank split is invariant under `d -> N/d`. For `p^3` the rank split is `1-1`, not `1-0-1`.
* **The five pairs `(d, N/d)`.** By class:
  * `(p, pqr)`: prime + composite;
  * `(q, p^2r)`, `(r, p^2q)`: prime + square;
  * `(p^2, qr)`: square + composite;
  * `(pq, pr)`: composite + composite.

  The swap `q <-> r` fixes three pairs and exchanges two. The divisor lattice has no halving/doubling cocycle telling `q` from `r`. They are differentiated only by external data, e.g. their sides mod 6 in the sandwich census: on `6k - 1`, `q` and `r` have opposite residues `±1`. That is a "sign" differentiation, and it is the mod-3 one.

*"`p^2qr` is `p` and `p^3` combined": ANALOGY.* The exact content found:
1. `Omega(p^2qr) = Omega(p) + Omega(p^3)`.
2. The exponent identity `(2,1,1) = (3,0,0) + (-1,1,1)`. This is the owner's `s^3(ab/s)`: `p^3` times the "virtual prime" `qr/p`, with `Omega = 1`.
3. The same identity relates the solutions at every `k`-free level. The solution set `{p, p^(k+1), p^k qr}` (plus `p^k q^2` for `k >= 3`) of `F = S_k + U` satisfies `p^k qr = p^(k+1)·(qr/p)` for all `k`. This is the precise "fractal" (`k`-indexed self-similar) structure. It is a pattern in the solution list, not a mechanism: `D = F - S - U` is not additive in exponent vectors.
4. Mod 9, for primes `p >= 5`: `p^3 = (p/3) = ±1` (the cube is the sign) and `p^2 ∈ {1,4,7}` (the square hand). So `p^2qr = (p^2)(qr) (mod 9)` splits as square hand times `qr`.

No law was found in which `F, S, U` of `p^2qr` are computed from those of `p` and `p^3`.

**4.4 Verdict.** The triplet principle is exact wherever it reduces to `Z/3` with `x2 = -1`, or to a geometric weight `2^-v` split by the parity of `v`. Both hold in the Syracuse 3-adic structure (§5) and in Pythagorean triples. Outside that it is a useful heuristic, labelled ANALOGY. The mod-3 coincidence `2 = -1` is why mod 9 matters for Collatz: halving and sign reversal are the same operation mod 3. So the parity of the number of halvings *is* the residue mod 3 of the next odd iterate, `S(A) = (-1)^v (mod 3)`.

## 5. Grounding: the Syracuse clock

Notation: `S(A) = (3A + 1)/2^v`, `v = v_2(3A+1)`, for `A` odd. Tao (arXiv:1909.03562, Forum Math. Pi 2022) defines `Syrac(Z/3^n)` by (1.22) with iid `Geom(2)` valuations, shows the consistency (1.23), and gives the recursion Lemma 1.12. He computes the mod-9 law

`(0, 8/63, 16/63, 0, 11/63, 4/63, 0, 2/63, 22/63)`

and in Remark 1.13 identifies the law of `Syrac(Z_3)` as the stationary measure `mu` of `x -> (3x+1)/2^a` with probability `2^-a` (CITED). By (1.21) and Proposition 1.9, `Syr^n(N) mod 3^k` is close in total variation to `Syrac(Z/3^k)` for typical `N` (CITED).

### 5.1 Clock coordinates (PROVED)

Every unit mod `3^k` is uniquely `x = (-1)^c (-2)^l`, with `c ∈ Z/2` (the cube/sign hand, `x = (-1)^c mod 3`) and `l ∈ Z/3^(k-1)` (the square hand; `-2` generates the 1-units `1 + 3Z/3^k`). Mod 9 one may use `4 = (-2)^2` instead:

| `x mod 9` | 1 | 4 | 7 | 8 | 5 | 2 |
|---|---|---|---|---|---|---|
| `(c, s)` | (0,0) | (0,1) | (0,2) | (1,0) | (1,1) | (1,2) |

Here `x = (-1)^c 4^s`.

### 5.2 The step law (PROVED; VERIFIED)

For `A` a unit mod 3: `3A + 1 = 4^(1+c) (mod 9)` and `2^-1 = 5 = -4 (mod 9)`. Hence

```text
S(A) = (-1)^v 4^(1 + c + v)   (mod 9),   i.e.  (c, s) -> (v mod 2, 1 + c + v mod 3).
```

* **Erasure and lockstep.** The new state does not depend on `s`: the `+1` erases the square hand, and the `v` halvings then advance both hands in lockstep. (For `3 | A`, `(c', s') = (v mod 2, v mod 3)`.)
* **Mod `3^k`.** The general form is `sign(S(A)) = (-1)^v` and `l(S(A)) = lambda(A) - v`, with `lambda(A) = log_(-2)(1 + 3A) = -A (mod 3)`. Here `lambda(A) mod 3^(k-1)` depends only on `A mod 3^(k-1)`: the hand loses one 3-adic digit per step.
* **Consequence (classical).** `A_(n+1) mod 9` is a function of `(v_(n-1) mod 2, v_n mod 6)`. More generally `A_(n+1) mod 3^k` is a function of the last `k` valuations. This is the formula behind (1.21).
* **The transition digraph.** The mod-9 transition digraph of the residues is an out-twin blow-up over the sign bit. The three states with the same `c` have identical out-rows, as in `C_3[3K_1]`. So the transition matrix has rank 2, and `P^2 = 1 pi` exactly: every row of `P^2` is Tao's law (checked in exact arithmetic). The chain is exactly mixed after two steps.

VERIFIED: the step law on all 1,000,000 odd `A < 2·10^6`, and the two-hand formula on 2,200,000 orbit transitions. Empirically the law of `S^6(A) mod 9` for 200000 random 80-bit `A` matches Tao's law within sampling error (TV distance `0.0048` in the runner's seeded sample).

**Coupling of the two hands (PROVED).** Tao's mod-9 law is not a product of its marginals: `P(1) = 8/63` while `P(c=0) P(s=0) = 10/63`. The coupling is one tick: `P(s | c = 1) = P(s + 1 | c = 0)`, i.e. `(11, 2, 8)/21` versus `(8, 11, 2)/21`. The reason: the odd valuations `2i - 1` carry exactly twice the weight of the even valuations `2i`, so `Law(v mod 3 | v odd) = Law(v + 1 mod 3 | v even)`.

### 5.3 The doubling law (PROVED + FINITE-EXACT)

**Theorem G2.** Let `mu` be the law of `Syrac(Z_3)`. Then
1. `mu(3Z_3) = 0`;
2. `mu(1 + 3Z_3) = 1/3` and `mu(2 + 3Z_3) = 2/3`;
3. for every Borel `E ⊆ 1 + 3Z_3`, `mu(2E) = 2 mu(E)`.

Equivalently, `(x2)_* mu = (mu + nu)/2` with `nu = Law(3·Syrac + 1)`.

*Proof.* Stationarity gives `mu = sum_(a>=1) 2^-a (x 2^-a)_* nu`. Shifting the index, `(x2)_* mu = nu/2 + mu/2`. On `2 + 3Z_3`, `nu` vanishes, and `x2` maps `1 + 3Z_3` onto `2 + 3Z_3`. ∎

FINITE-EXACT: `P(2x) = 2 P(x)` for every `x = 1 (mod 3)` modulo `3^k`, `k <= 8`, using the exact laws from Lemma 1.12 in integer arithmetic. For example, mod 27: `P(1) = 1376/37449` and `P(2) = 2752/37449`.

This is the owner's triplet as a theorem about the Syracuse measure. The class `0` is primarily differentiated (it is never hit). The two unit classes are the doubled images of each other, with masses in ratio `1 : 2`.

### 5.4 The square-clock IFS (PROVED + FINITE-EXACT)

**Theorem G3.** `mu_1 := 3 mu|_(1+3Z_3)` is the unique invariant probability measure of the iterated function system on `1 + 3Z_3`:

```text
x -> x/4   (prob 1/4),     x -> (3x+1)/4   (prob 1/4),     x -> (6x+1)/4   (prob 1/2).
```

Equivalently, `4 (x4)_* mu_1 = mu_1 + (3x+1)_* mu_1 + 2 (6x+1)_* mu_1`.

*Proof.* Restrict `2(x2)_* mu = mu + nu` (G2) to `1 + 3Z_3`, and use G2 again to replace `mu|_(2+3Z_3)` by `2 (x2)_* mu|_(1+3Z_3)`. Then `nu = (3x+1)_* mu|_1 + 2 (6x+1)_* mu|_1`. Dividing by 4 and pushing forward by `x -> x/4` gives the IFS equation.

Uniqueness, by induction on the level `k`. Let `mu', mu''` be two invariant probability measures that agree on the classes mod `3^(k-1)`, and let `E` be a class mod `3^k`. The preimages of `E` under `(3x+1)/4` and `(6x+1)/4` are classes mod `3^(k-1)`. So `delta = mu' - mu''` satisfies `delta(E) = delta(4E)/4`, hence `|delta(E)| <= 4^-j |delta(4^j E)| <= 2 4^-j -> 0`. ∎ (Identity and uniqueness PROVED; the identity is FINITE-EXACT to `3^8`.)

**Reading.** The class `1 + 3Z_3 ≅ Z_3` (log base 4) is the square clock.
* `x/4` turns it back one tick without carrying: a 3-adic isometry.
* The two carries `(3x+1)/4` and `(6x+1)/4` land in `1 + 9Z_3` and `4 + 9Z_3` and contract by 3.

At level 9 the IFS is the 3-state clock

`mu_1(1) = mu_1(4)/4 + 1/4`, `mu_1(4) = mu_1(7)/4 + 1/2`, `mu_1(7) = mu_1(1)/4`,

whose solution `(8, 11, 2)/21`, doubled onto class 2 by G2, is exactly Tao's `(8, 16, 11, 4, 2, 22)/63`.
* `7` is fed only by rotation (smallest mass).
* `1` and `4` are fed by the carries with weights `1/4` and `1/2`: one tick of doubling apart.

At level `3^k` the IFS reads `mu_1(x) = mu_1(4x)/4 + g(x)`, where `g(x) = mu_1(f_B^-1 x)/4 + mu_1(f_C^-1 x)/2` depends only on level `k - 1`. Going once around the square clock (4 has order `L = 3^(k-1)` on `1 + 3Z/3^k`) gives `mu_1(x) = sum_(j=0..L-1) 4^-j g(4^j x) / (1 - 4^-L)`: a geometric average around the square clock, driven by the carries of the previous level. This is the measure-side twin of THM-4520's circulant `(S/2)(I - S/2)^-1`. On the Fourier side G2 is the stationarity equation in the form `mu-hat(2 xi) = (mu-hat(xi) + e(xi/3^n) mu-hat(3 xi))/2` (`mu-hat(xi) = E e(xi Syrac/3^n)`). That is equivalent to Tao's Lemma 1.12 and to the repo's frequency recursion (THM-4519/4520). G2 and G3 add an exact reading; they do not add a bound.

### 5.5 Where TCPC touched open problems

1. **HYP-9167(b) (tournaments).** The two-sheet clock lens produced the first all-odd tournaments on 14, 34 and 38 vertices (T2), and a conjectural uniform construction for every `N = 2 (mod 4)` (C1, verified to `N = 38`). This is the lane's genuinely new result on an open problem.
2. **Collatz 3-adic structure.** G1–G3 are exact, verified reformulations of the Syracuse measure. They are not progress on Collatz: the fine-scale mixing (Tao's Proposition 1.14) and H1 (HYP-9166) are what matter, and the clock picture does not bound them. What the clock picture *does* explain exactly:
   * why the mod-9 law has the shape it has: rotation plus two carries, doubled;
   * why class 0 is empty;
   * why classes 1 and 2 are related by doubling.

## 6. Plain-language summary for the owner

* **Clocks.** Take the numbers 0 to `m-1` around a clock face and a list of allowed jumps; draw an arrow for each allowed jump. Squares mod 9 (`1, 4, 7`) give three groups of three, each pointing at the next group. Cubes mod 9 (`±1`) give a ring where every arrow goes both ways. Allowing both gives everything except jumps by 3.
* **Five combination rules.** Squares rotate, cubes double, multiples of 3 nest one clock inside another. Every clock with an odd number of positions has a mirror symmetry that makes every arrow lie on an even number of complete routes, whatever the jump list.
* **Two-sheet clocks and the open problem.** A clock with two sheets is a clock where every number carries an extra even/odd tag, which is your halving/doubling bit. A two-sheet clock on 14 points turns out to have every arrow on an odd number of complete routes. Nobody had such an example: the classical Paley construction needs 15 to be a prime power, and it is not. The same recipe works for every size tested up to 38 points, including 34 and 38, where 35 and 39 are not prime powers either.
* **How to draw it.**
  1. Take the 7 hour marks of a 7-hour clock, then the same 7 marks turned a quarter turn: 14 points on a circle.
  2. Draw an arrow from each point to every point less than half a turn ahead (counterclockwise).
  3. Reverse the arrows among the turned points.

  For 34 points use a 17-hour clock, and because 17 leaves remainder 1 when divided by 4, also reverse the 17 arrows joining each mark to its own turned copy.
* **Mod 9 and Collatz.** After each Collatz step:
  * the residue mod 3 is set by whether you halved an even or odd number of times;
  * the finer mod-9 digit is wiped clean by the "+1" and then turned by the number of halvings.

  The long-run frequencies of `1, 2, 4, 5, 7, 8 mod 9` (Tao computed them: 8, 16, 11, 4, 2, 22 out of 63) come from a three-position clock fed by two "carries", one twice as heavy as the other. That is your triplet: one position fed only by turning, two fed by carries that differ by a doubling.
* **What is solid and what is not.** The tournament result and the Collatz clock laws are exact. The claim that `p^2qr` is "`p` and `p^3` combined" stays an analogy: the arithmetic is `p^2qr = p^3 × (qr/p)`, and no counting law was found behind it. Nothing here proves Collatz.

## 7. Open problems and next steps

1. **C1:** prove `D_r` all-odd for all odd `r`.
   * The first unverified case is `r = 21` (`N = 42`, about `2·10^11` subset representatives, roughly a day for karp4).
   * The next new existence case, `N = 50`, is out of computational reach.

   A proof of C1 would settle the existence half of HYP-9167(b).
2. Explain the extra all-odd classes at `r = 11` (`X = {0,1,4,5,9}`, `{0,1,3,7,8}`); at `r = 5, 7, 9, 13`, `D_r` is the only rotational-sheet class. Is the number of all-odd two-sheet classes on `Z/2r` unbounded?
3. **HYP-9167(b) non-existence at `N = 13`.** The `Z/11`, `Z/7` and `Z/5` families are searched here. The order-3 family (orbits 3,3,3,3,1; about `6.7·10^7` members) remains, as does a proof.
4. **L1 for even order.** The middle-arc count `M(e)` is proved and verified. Find the even-order clocks where it is computable in closed form, e.g. two-sheet clocks viewed on `Z/2r`.
5. **G3.** Does the IFS form (one isometry, two contractions) give a cheaper proof of any piece of Tao's fine-scale mixing? The isometry `x/4` is exactly what prevents a naive contraction argument at fine scales.
6. **`Q_r`.** Why is the apex completion `Q_r` all-even for `r <= 9` without vertex-transitivity?

## 8. Reproduction

```bash
python3 04-computation/experiments/procgen_tcpc_20261001_run.py --long --n26 --n34 > 05-knowledge/results/procgen_tcpc_20261001.out
python3 04-computation/experiments/procgen_tcpc_20261001_run.py --only-n38 > 05-knowledge/results/procgen_tcpc_20261001_n38.out
```

The runner compiles its C engines into `scratch/procgen_tcpc/` and re-verifies every FINITE-EXACT claim above. It prints to stdout only and ends with `ALL CHECKS PASSED`.

| mode | adds | time | checks |
|---|---|---|---|
| default | everything except the items below | about 4 min | about 380000 |
| `--long` | the `r = 11` cross-set census, the `N = 13` structured searches, `F(P_23)`, the `m = 23` circulant scan, `D_15` by karp3 and karp4 | | |
| `--n26` | `D_13` and its interval control by the lean subset DP, on `T` and on `T^op` | | |
| `--n34` | `D_17` by karp3 on `T` and karp4 on `T^op` | | |
| `--long --n26 --n34` | all of the above (the committed `.out`) | 3698 s on a loaded machine | 671807 |
| `--only-n38` | only `D_19` by karp4 (`procgen_tcpc_20261001_n38.out`) | 4352 s | 2 |

Peak memory is below 300 MB (the lean subset DP uses 257 MB; the inclusion-exclusion engines use about 1 MB). Engines:
* `procgen_tcpc_20261001_hp.c`: modes `exact`, `par`, `rooted`, `batchpar`;
* `procgen_tcpc_20261001_par26.c`;
* `procgen_tcpc_20261001_karp3.c`;
* `procgen_tcpc_20261001_karp4.c` (ARM64 with the crypto extension, for `PMULL`).

Library: `04-computation/experiments/procgen_tcpc_20261001_lib.py`.

SHA-256 (raw bytes):

| file | SHA-256 |
|---|---|
| `procgen_tcpc_20261001.out` | `5a76471fbb6d432e6103e5194812f33007f083d70435d1f45d42ad47869d4d02` |
| `procgen_tcpc_20261001_n38.out` | `983db6ecd24ca20fc3b996ae03bfd479b0642e9fd89ae6abf0c9d41cc30788f4` |
| `procgen_tcpc_20261001_run.py` | `7ac530938f29a0cfc5707569b8e17b560b343bfc816de1553ec2f09aff304add` |
| `procgen_tcpc_20261001_lib.py` | `1e639a74c3ae21ec4dc2f0b347ed78feee82a682eb8b058abe53f31216802247` |
| `procgen_tcpc_20261001_hp.c` | `3e5225e75e6dd86ab28a3e1df140a72d1081b2bf5c25391dde3860cbe2d939e4` |
| `procgen_tcpc_20261001_par26.c` | `1ac1c552a35b1d8c5435017bafbe1ffa4b84fb6726365735fe9731d737502f1e` |
| `procgen_tcpc_20261001_karp3.c` | `655b1911e607c45db17e719f5c4ea0d954fe8db62d1d17a03b6a7af15fa3c815` |
| `procgen_tcpc_20261001_karp4.c` | `2d004943351ef91ef126e650edb8b108ebaa696669f46f5358758eda8b8ad45c` |

*Process note.* The runner caught one false draft statement: the one-term cross-arc parity formula for lexicographic products (§2.2), now corrected and recorded as REFUTED. One runner bug was also fixed: the `r = 13` interval control had used `D_13`'s arc representatives.
