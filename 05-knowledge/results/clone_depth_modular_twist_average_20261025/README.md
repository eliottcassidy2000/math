# Clone-fibre cancellation depth modulo primes, and the twist average of the hourglass root surface

2026-10-25, klein (session worktree `codex/session-tournament-twist-unify`).
**PROVED** elementary theorems with complete proofs below; **FINITE-EXACT** receipts
(`clone-normal.json`, `twist-normal.json`, `census_*.out`) in both normal and `-O` modes;
one **OPEN** hypothesis ([HYP-9374](../../hypotheses/HYP-9374-mixed-order-clone-depth-equals-the-candidate-minimum.md)).
Canon entries: [THM-4764](../../../01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md)
and [THM-4765](../../../01-canon/theorems/THM-4765-quadratic-twists-of-a-square-root-cover-average-to-the-branch-complement.md).
No historical-novelty claim, no faster enumerator, no LRC/JC/HN transfer.

## 0. Inheritance, board, and the owner's claims

The owner's prompt described "one 50-vertex tournament" with one independent symmetry-respecting
toggle, cancellation depth 5 over the rationals, 21 modulo 5 and 25 modulo 2, and, separately,
a twist calculation averaging the 64 quadratic twists to `q^2-6q+12 = (q-1)^2-4(q-1)+7`.
Every one of those numbers is now verified and placed inside a proved law. The objects are:

- **Clone fibre.** `ports_regular_actions_hn_20261025/counting/CLONE_INVOLUTION.md` identifies the
  weighted-inflation fibre of a tournament `T` with the orbits of cut subsets of its canonical module
  gaps `G` under `A = Aut(T)` (input THM-4463 §2). `session62_agent_theme_bridges_20261025/clone/README.md`
  proves that the subset-orbit polynomial `P_A(t) = sum_[C] t^|C|` vanishes at `t=-1` to order
  `d = min_g c(g)` (minimum cycle count of one element), that the rank of legal invariant toggles is the
  number `r` of `A`-orbits on `G`, and realizes the regular `C3^2` action by a 24-vertex tournament.
- **Twists.** Session59's `toric-shared-roots/TORIFICATION_OBSTRUCTION.md` (Theorem T1) and Session60/the
  Gale-cube packet treat the hourglass root surface `H` (seven square roots of `E, E±x, E±z, E±(x+z)`),
  `p_g = 7`, `K^2 = 16`, 24 nodes, `e(Y) = 80`; Session59's `arithmetic-configurations/RESULT.md` §4 and
  `graph-magic-degeneration/branch-addendum.md` §4 use negative `(q-1)`-coefficients to exclude
  split-torus stratifications.

Closest proved mechanism: Session62's exact root order (2)–(3). Canonical hostile: the transposition
(even action: `1+t+t^2`, no zero). Corrected near miss: the single labelled gap toggle (MISTAKE-free but
refuted in CLONE_INVOLUTION §5). Least-used sidecar: the **prime-by-prime reduction** of the integer
polynomial `P_A`, never examined before in the repository (bounded search: no mod-`p` vanishing-order
statement for orbit or Burnside polynomials was found; 2-adic valuations of A000568 are THM-305).

| Live concept | Object / operation | What this packet changes | Lost coordinate / test |
|---|---|---|---|
| Clone fibre polynomial | `P_A(t)` of the gap action; expand at `t=-1` | depth mod `p` is a new invariant: `(5,21,25)` versus `(1,5,25)` for the two regular groups of order 25 | depends only on the gap action; two different Cayley sets give the same `P` |
| Cyclic-subgroup counts | `c_j` = number of cyclic subgroups of order `p^j` | `J = max{j : p ∤ c_j}` controls depth mod `p` | abelian non-cyclic forces `J=1`; non-abelian uses a classical count |
| Mixed-order groups | `p'`-part census `E_{o'}` and `e_{j,o'}` | candidate-order lower bound, equality in all 22 tested groups | collisions of candidate orders are the open case |
| Quadratic twists | `(k^×/k^{×2})^n / diag`, `2^(n-1)` classes | twist average = branch complement; Walsh spectrum = even character sums | the average forgets all `S_J`, `J ≠ ∅`; retain the additive relations |
| Critical expansion | `t=-1` or `q=1` | one formal lemma (jump iff `ℓ` divides the leading coefficient) serves both | no map between tournament and surface; only the method transfers |

Anchor: verify and prove the owner's numbers. Niche: the modular depth law. Wildcard: the Walsh
spectrum of the twisted counts.

## 1. The two 50-vertex tournaments

Let `S = {(0,1),(0,2),(1,0),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(3,1)} ⊂ Z5^2`, a half-system
(`S ∩ -S = ∅`, `S ∪ -S = Z5^2 \ 0`), and `R = Cay(Z5^2, S)` the Cayley tournament `v -> v+s`.
Inflate every vertex to the transitive pair `TT2`: `T = R[TT2]`, 50 vertices. The cyclic twin uses
`R' = Circ(Z25, {1,...,12})` and `T' = R'[TT2]`.

| | `T = Cay(Z5^2,S)[TT2]` | `T' = Circ(Z25,{1..12})[TT2]` |
|---|---|---|
| vertices / `|Aut|` (enumerated, no group assumed) | 50 / 25 | 50 / 25 |
| homogeneous pairs = canonical modules | exactly the 25 clone blocks | exactly the 25 clone blocks |
| gap action | regular `Z5^2` (translations) | regular `Z25` (rotations) |
| invariant toggle rank `r` | 1 | 1 |
| cycle profiles | `1^25` once, `5^5` 24 times | `1^25` once, `25^1` 20 times, `5^5` 4 times |
| `|Ω|` subset orbits (Burnside and exact `2^25` census agree) | 1,342,208 | 1,342,184 |
| depth over `Q`, leading coefficient `κ` | 5, `κ = 3000 = 2^3·3·5^3` | 1, `κ = 20 = 2^2·5` |
| depth mod 2, 3, 5, 7, 11, 13 | 25, 25, 21, 5, 5, 5 | 25, 1, 5, 1, 1, 1 |

The polynomial of `T` is

```
P(t) = ((1+t)^25 + 24 (1+t^5)^5)/25
     = 1 + t + 12t^2 + 92t^3 + 506t^4 + 2130t^5 + 7084t^6 + 19228t^7 + 43263t^8 + 81719t^9
       + 130760t^10 + 178296t^11 + 208012t^12 + 208012t^13 + ... + t^25,
P(-1+x) = 3000x^5 - 30000x^6 + 150000x^7 - 495000x^8 + ... + x^25.
```

Both tournaments have the same order, automorphism order, number of gaps and toggle rank; the depth
triple separates them. By `center_transitivity_20261021c` §3 every vertex-transitive 25-tournament is a
Cayley tournament on one of these two groups, so for every rigid such tournament the inflated clone
fibre has one of these two polynomials (3840 of the 4096 half-systems of `Z5^2`, and 4080 of the 4096
of `Z25`, are pair-free with automorphism group exactly the regular group; the others are not used).

Everything in the table is computed from literal arcs: colour-refinement/individualisation
enumeration of all automorphisms, a literal homogeneous-pair search, the induced permutations of the
25 first gaps, Burnside from those permutations, and (independently) the `2^25` minimal-mask census
in `census_orbits.c`.

## 2. Depth modulo primes

Notation. `A` finite of odd order acting on a finite set `G`, `n = |G|`, `P_A(t) = sum_[C] t^|C| ∈ Z[t]`,
`s = 1+t`, `Q_j = [x^j] P_A(-1+x)`, `d = ord_{t=-1} P_A = min_g c(g)` and `κ = Q_d > 0` (Session62).
For a prime `ℓ`, the **depth modulo `ℓ`** is `d_ℓ = min{ j : ℓ ∤ Q_j }`, the order of vanishing at
`t=-1` of `P_A` reduced in `F_ℓ[t]` (finite, because `P_A` is monic).

**Lemma B0.** For every prime `ℓ`: `d_ℓ ≥ d`, with equality iff `ℓ ∤ κ`.
*Proof.* `P_A = s^d R(s)` with `R ∈ Z[s]`, `R(0) = κ`; reducing mod `ℓ` multiplies the order by nothing
and adds `ord_s(R mod ℓ)`, which is zero iff `ℓ ∤ κ`. ∎

For a **regular** action (`G = A`), `P_A(t) = (1/|A|) sum_o e_o (1+t^o)^{|A|/o}` where `e_o` is the number
of elements of order `o`.

**Theorem B′ (primes prime to the order).** If `ℓ ∤ |A|` and the action is regular, then
`d_ℓ = min{ |A|/o : ℓ ∤ e_o }`.
*Proof.* `1/|A|` is an `ℓ`-adic unit. For each `o`, `1+t^o = o·s + O(s^2)` with `ℓ ∤ o` (as `o` divides
`|A|`), so `(1+t^o)^{|A|/o}` has exact order `|A|/o` mod `ℓ` with unit leading coefficient `o^{|A|/o}`.
The `o`-term vanishes identically mod `ℓ` iff `ℓ | e_o`. Distinct `o` have distinct orders `|A|/o`, so no
cancellation occurs among surviving terms. ∎

**Theorem B (`p`-groups).** Let `|A| = p^e`, `p` odd, act regularly, and let `c_j` be the number of cyclic
subgroups of order `p^j`. Then `c_1 ≡ 1 (mod p)`, and with `J = max{ j ≥ 1 : p ∤ c_j }`,

```
d_p = p^e - p^J + p^(J-1).
```

In particular `A` cyclic gives `J = e` and `d_p = p^(e-1)`; `A` abelian non-cyclic gives `J = 1` and
`d_p = p^e - p + 1`.

*Proof.* Write `P = (1/p^e) sum_{j=0}^{e} c_j φ(p^j) (1+t^{p^j})^{p^{e-j}}`, `c_0 = 1`. For `j ≥ 1`,
`1 + t^{p^j} = sum_{i=1}^{p^j} C(p^j,i)(-1)^{p^j-i} s^i = s^{p^j} + p·s·g_j(s)` with `g_j ∈ Z[s]`; by
Kummer `v_p(C(p^j,i)) = j - v_p(i)`, so modulo `p` the first surviving term of `g_j` is at `i = p^(j-1)`
with unit coefficient `C(p^j,p^(j-1))/p`: `g_j ≡ u_j s^{p^(j-1)-1} + O(s^{p^(j-1)})`, `u_j` a unit, and
`g_1(0) = 1`. Expand `(s^{p^j} + p s g_j)^{p^{e-j}} = sum_m C(p^{e-j},m) p^m s^m g_j^m s^{p^j(p^{e-j}-m)}`.
The coefficient `c_j φ(p^j) C(p^{e-j},m) p^m / p^e` has `p`-adic valuation
`v_p(c_j) + (j-1) + (e-j-v_p(m)) + m - e = v_p(c_j) + m - 1 - v_p(m)`, which is at least 1 for `m ≥ 2`
(`p` odd). Hence modulo `p`

```
P ≡ s^{p^e} + sum_{j≥1} c_j (p-1) s^{p^e - p^j + 1} g_j(s)      (mod p),
```

the `m=0` terms summing to `s^{p^e}` because `sum_j c_j φ(p^j) = p^e`. The `j`-summand vanishes iff
`p | c_j`, and otherwise has exact `s`-order `p^e - p^j + p^(j-1)`, which is strictly decreasing in `j`
and smaller than `p^e`; distinct `j` cannot cancel. So `d_p = p^e - p^J + p^(J-1)` with `J` the largest
admissible index. `J` exists: reducing `sum_{j≥1} c_j φ(p^j) = p^e - 1` modulo `p` gives `-c_1 ≡ -1`.
For cyclic `A`, `c_j = 1` for all `j`. For abelian non-cyclic `A = ⊕ C_{p^{a_i}}` with `r ≥ 2` summands,
the elements of order exactly `p^j` (`j ≥ 2`, `j ≤ a_1`) number `|A[p^{j-1}]|(p^{m_j}-1)` with
`m_j = #{i : a_i ≥ j} ≥ 1`, so `v_p(c_j) = sum_i min(a_i, j-1) - (j-1) = sum_{i≥2} min(a_i,j-1) ≥ 1`;
thus `J = 1`. ∎

For non-abelian non-cyclic `p`-groups the same conclusion `J = 1`, hence `d_p = p^e - p + 1`, follows
from the classical fact that the number of cyclic subgroups of order `p^j ≥ p^2` in a non-cyclic
`p`-group of odd order is divisible by `p` (G. A. Miller; P. Hall 1933; stated in Berkovich, *Groups of
prime power order I*, §1 — **CITED from a secondary source, not re-derived here**). It is verified below
on all five groups of order 27 and on `C27⋊C3`, `C9⋊C9`, `Heis(5)`, `C25⋊C5`.

**Theorem B″ (general regular actions, lower bound).** Let `|A| = p^e m`, `p ∤ m`, act regularly. For each
divisor `o'` of `m` let `E_{o'}` be the number of elements whose `p'`-part has order `o'`, and for `j ≥ 1`
let `e_{j,o'}` be the number of elements of order `p^j o'`. Then `p^e | E_{o'}` and `p^(j-1) | e_{j,o'}`, and

```
d_p ≥ min( { |A|/o' : p ∤ E_{o'}/p^e }  ∪  { |A|/o' - p^j + p^(j-1) : p ∤ e_{j,o'}/p^(j-1) } ).
```

*Proof.* Group the terms of `P` by the `p'`-part `o'` of the element order and put `s' = 1 + t^{o'}`,
which equals `o'·s·(unit power series)` since `p ∤ o'`. The term of order `p^j o'` is
`(e_{j,o'}/|A|) (s'^{p^j} + p s' g_j(s'))^{M}`, `M = p^{e-j} m / o'`, and the binomial expansion in `μ`
has coefficients of valuation `v_p(e_{j,o'}) - (j-1) + (μ - v_p(μ))`, using `v_p(C(M,μ)) = e - j - v_p(μ)`
for `p ∤` the cofactor. Divisibility `p^(j-1) | e_{j,o'}` holds because elements of order `p^j o'` come in
generator sets of cyclic subgroups, of size `φ(p^j)φ(o')`. Hence all `μ ≥ 2` terms vanish mod `p`, the
`μ = 1` terms survive exactly when `v_p(e_{j,o'}) = j - 1`, at `s`-order `|A|/o' - p^j + p^(j-1)`, and the
`μ = 0` terms of all `j` with the same `o'` add up to `(E_{o'}/|A|) s'^{|A|/o'}` times a unit series.
Divisibility `p^e | E_{o'}`: the elements with `p'`-part conjugate to a fixed `h` of order `o'` number
`[A : C_A(h)]` times the number of `p`-elements of `C_A(h)`, and by Frobenius' theorem the latter is a
multiple of `|C_A(h)|_p`. Every surviving term therefore has `s`-order in the displayed candidate set, so
the order of `P mod p` is at least its minimum. ∎

Equality fails only if candidate orders collide and cancel; this never happens in the 22 groups tested
(five mixed-order groups and seventeen `p`-groups, where Theorem B makes it a theorem). The equality
statement is [HYP-9374](../../hypotheses/HYP-9374-mixed-order-clone-depth-equals-the-candidate-minimum.md).
The **direct-product** case `A = S × H` (`S` the Sylow `p`-subgroup) reads
`d_p = p^e m / o* - p^J + p^(J-1)` with `o* = max{ o' : p ∤ e^H_{o'} }` and `J = J(S)`, provided no
collision; for `C7⋊C3` the candidate formula is needed in its general form (`E_3 = 14`, `v_7 = 1 = e`).

### Table (all regular; from `clone-normal.json`)

| group | `|A|` | `d_Q` | `d_2` | `d_3` | `d_5` | `d_7` | `d_11` | `d_13` | `J` (Thm B) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `C3^2` | 9 | 3 | 9 | **7** | 3 | 3 | 3 | 3 | 1 |
| `C9` | 9 | 1 | 9 | **3** | 1 | 1 | 1 | 1 | 2 |
| `C5^2` | 25 | 5 | 25 | 25 | **21** | 5 | 5 | 5 | 1 |
| `C25` | 25 | 1 | 25 | 1 | **5** | 1 | 1 | 1 | 2 |
| `C7^2` | 49 | 7 | 49 | 49 | 7 | **43** | 7 | 7 | 1 |
| `C27` | 27 | 1 | 27 | **9** | 1 | 1 | 1 | 1 | 3 |
| `C9×C3`, `M27 = C9⋊C3` | 27 | 3 | 27 | **25** | 3 | 3 | 3 | 3 | 1 |
| `C3^3`, `Heis(3)` | 27 | 9 | 27 | **25** | 9 | 9 | 9 | 27 | 1 |
| `C27⋊C3` | 81 | 3 | 81 | **79** | 3 | 3 | 3 | 3 | 1 |
| `C9⋊C9` | 81 | 9 | 81 | **79** | 9 | 9 | 9 | 9 | 1 |
| `C125` | 125 | 1 | 125 | 1 | **25** | 1 | 1 | 1 | 3 |
| `Heis(5)` | 125 | 25 | 125 | 25 | **121** | 25 | 25 | 25 | 1 |
| `C25⋊C5` | 125 | 5 | 125 | 5 | **121** | 5 | 5 | 5 | 1 |
| `C15` | 15 | 1 | 15 | 1 | 1 | 1 | 1 | 1 | mixed |
| `C21` | 21 | 1 | 21 | **19** | 1 | 1 | 1 | 1 | mixed |
| `C45` | 45 | 1 | 45 | **3** | 1 | 1 | 1 | 1 | mixed |
| `C7⋊C3` | 21 | 3 | 21 | 3 | 3 | **7** | 3 | 3 | mixed |
| `C3×C5^2` | 75 | 5 | 75 | **73** | **21** | 5 | 5 | 5 | mixed |

Bold entries are at primes dividing `|A|`. Primes not dividing `|A|` follow Theorem B′ exactly (for
example `C3^3` mod 13 jumps to 27 because `e_3 = 26 ≡ 0 (mod 13)`; `C5^2` mod 2 and mod 3 jump to 25
because `e_5 = 24`). The `C21` row is the witness that mixed order is genuinely different from the
`p`-group case: modulo 3 the `s^3` coefficients `760` (order-21 elements) and `98` (order-7 elements)
cancel, `858 ≡ 0`, and the depth climbs to `19 = 21 - 3 + 1`.

Numerology flag: for `C_p^2` the depth mod `p` is `p^2 - p + 1 = Φ_6(p)` (`7, 21, 43`). For `e ≥ 3` it is
`p^e - p + 1`, not cyclotomic; the equality with `Φ_6` is an artefact of `e = 2`.

### What the modular depth means in the fibre

By Session62 (7), the signed binomial moments of the cut-orbit sizes,
`sum_[C] (-1)^(n-|C|) C(|C|,j)`, vanish for `j < d` and equal `κ` at `j = d`. Lemma B0 and the theorems
say that these integers are divisible by `ℓ` for all `j < d_ℓ`: for `T`, every signed moment through
`j = 20` is a multiple of 5 and every one through `j = 24` is a multiple of 2 and of 3. The `ℓ`-adic
valuations of the whole expansion (receipt field `valuations`) are a finer invariant; for `T` at `ℓ = 5`
they read `3,4,5,4,3,4,3,3,2,2,2,3,1,1,1,3,0,...` from `j = 5`, reaching zero exactly at `j = 21`.

## 3. The twist average

Let `k = F_q`, `q` odd, and `L_0,…,L_{n-1}` linear forms on `k^{d+1}` with no common zero in `P^d`. The
root variety `X ⊂ P^{n-1}` is the inverse image of the linear `P^d` under coordinate squaring. For
`t ∈ (k^×)^n` the twist `X_t` replaces `x_j^2 = L_j` by `x_j^2 = t_j L_j`. Twists by `t` and `t'` are the
same cover of `P^d` when `t'_j/t_j` are squares, and a common scalar does nothing projectively; so the
twist classes form `(k^×/k^{×2})^n/diag ≅ F_2^n/F_2`, of size `2^(n-1)` (`64` for the hourglass).
Let `U ⊂ P^d` be the complement of the branch hyperplanes and `X_t^{ét}` the preimage of `U`.

**Theorem C.** `sum_t #X_t(k) = 2^(n-1) · #P^d(k)` and `sum_t #X_t^{ét}(k) = 2^(n-1) · #U(k)`, the sums over
the `2^(n-1)` twist classes. Moreover, on `U` the twisted counts are a Walsh transform of the even
character sums: with `χ` the quadratic character and `S_J = sum_{y∈U} prod_{j∈J} χ(L_j(y))` for
`J ⊂ {0,…,n-1}` of even size,

```
#X_t^{ét}(k) = sum_{|J| even} ( prod_{j∈J} χ(t_j) ) S_J,      S_∅ = #U(k).
```

*Proof.* Fix `y ∈ P^d(k)` with a lift `ỹ`, let `J(y) = {j : L_j(ỹ) ≠ 0}` and `m = n - |J(y)| < n`. The
`k`-points of `X_t` over `y` are the classes of `(λ, x)` with `λ ∈ k^×`, `x_j^2 = λ t_j L_j(ỹ)`, modulo
`(λ,x) ~ (μ^2 λ, μ x)`, a free action of `k^×`. For fixed `λ` the number of `x` is `prod_{j∈J(y)}
(1 + χ(λ t_j L_j(ỹ)))`, which is `2^{|J(y)|}` when all `χ(λ t_j L_j(ỹ))`, `j ∈ J(y)`, equal `+1` and `0`
otherwise. As `λ` runs over `k^×` exactly half the values have `χ(λ) = ε` for either sign. Hence the fibre
has `2^{n-m-1}` points if the vector `(χ(t_j L_j(ỹ)))_{j∈J(y)}` is constant and `0` otherwise; the
condition does not depend on the lift. Among the `2^n` sign vectors `(χ(t_j))_j`, those whose restriction
to `J(y)` equals `±(χ(L_j(ỹ)))_{J(y)}` number `2·2^m`; modulo the diagonal, `2^m` twist classes. So the
fibre counts of `y` summed over the twist classes give `2^m · 2^{n-m-1} = 2^{n-1}`, for every `y`. Summing
over `P^d(k)`, respectively over `U(k)` (where `m = 0`), gives the two identities. For `y ∈ U` the
indicator of "constant vector" is `2^{-(n-1)} sum_{|J| even} prod_{j∈J} χ(t_j)χ(L_j(ỹ))`, which is
well defined because `|J|` is even; multiplying by `2^{n-1}` and summing over `U` gives the Walsh
formula. ∎

The theorem is checked on the hourglass (`n = 7`), the six-root K3 subset `E, E±x, E±z, E+x+z` and all
nine magic forms, over `F_q` for `q = 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 25, 49` and the
characteristic-three controls `q = 3, 9`; the Walsh identity is checked for all 64 hourglass twists.

**Hourglass complement.** The seven lines have `t_2 = 12` double and `t_3 = 3` triple points away from
characteristics 2 and 3, so

```
#U(F_q) = q^2 + q + 1 - [7(q+1) - (t_2 + 2 t_3)] = q^2 - 6q + 12 = (q-1)^2 - 4(q-1) + 7.
```

In characteristic 3 two further triple points appear (`(1,1)` and `(-1,-1)` in the chart `E = 1`) and
`#U = q^2 - 6q + 10`; in characteristic 2 the seven lines collapse to four and `#U = q^2 - 3q + 3` (the
cover is then inseparable and the twist statement is void). Over `C`, `e(U) = 3 - 14 + 12 + 6 = 7`,
the value of the count polynomial at `q = 1`; in Session59's Euler stratification of the resolved
hourglass, `e(Y) = 2^6·7 + 2^5·(14-24-9) + 2^4·12 + 2^3·3 + 24 = 80`, the term `2^6·7 = 448` is exactly
the étale-sheet contribution whose finite-field shadow is the twist average. The coefficient `-4`
excludes a finite stratification of `U` by split tori (a disjoint union of `G_m^r`'s has count
`sum a_r (q-1)^r`, `a_r ≥ 0`; equality for infinitely many `q` is a polynomial identity), by the same
argument as Session59 §4. The surface itself is not torifiable for a stronger reason (`p_g = 7`,
Corollary T2); the complement statement is the one the twist calculation sees.

**The two sevens are different invariants.** For every 7-subset of the nine magic lines, `χ(U) = 10 - t_3`
while `p_g = 7` for all of them (Theorem T1). The labelled census is `t_3 = 2: 4, 3: 18, 4: 12, 5: 2`
subsets, so `χ(U) = 7 = p_g` holds for the hourglass only because it has exactly three triple points, and
fails for 18 of the 36 subsets. No map relates the two sevens; "another exact connection to seven" is
a coincidence of the incidence count, recorded here so that it is not re-promoted.

**What the twist average forgets, and what the additive equations retain.** The average is the `J = ∅`
Walsh coefficient; it depends only on the Kummer (multiplicative) construction and the branch locus. The
other 63 coefficients `S_J` are character sums of products of the forms, and they feel the linear
(additive) relations among the forms. Two exact facts, one proved and one finite:

- The symmetry group of the configuration (generated by `x ↔ z`, `(x,z) ↦ (-x,-z)` and `(x,z) ↦ (z,-x-z)`,
  order 12) fixes `E`, permutes the six outer forms transitively and the three parallel pairs
  transitively, and preserves `U` and every `S_J` up to relabelling. Hence the six "outer-deletion"
  sums `S_{J}`, `J = all roots except one outer root`, coincide, and the three parallel-pair sums
  `S_{{E+x,E-x}}, S_{{E+z,E-z}}, S_{{E+x+z,E-x-z}}` coincide. **PROVED.**
- For every tested `q ≥ 5` the central-deletion sum (`J = all six outer roots`) differs from the common
  outer value: `(3,7), (5,7), (9,31), (5,19), (37,19), (39,43), …` for `q = 5, 7, 11, 13, 17, 19, …`.
  **FINITE-EXACT.** This is the finite-field shadow of Session60's `1+6` split of the seven K3
  deletions (Picard ranks 20 central, 19 outer); no identification of the two statements is claimed.
- The parallel-pair sums are a period-12 quasi-polynomial in the data (`q`, `-(q-4)`, `-(q-8)`, `q-4` by
  residue of `q` mod 12 for `q ≥ 5`), as in Session59's period-12 counts. **DATA**, not a theorem.

## 4. One frame for both calculations

Neither object maps to the other. What they share is a formal structure, and one lemma applies
verbatim to both:

| | clone fibre | twist average |
|---|---|---|
| finite group averaged | `Aut(T)` on gaps (odd order) | twist classes `F_2^n/F_2` (a 2-group of characters) |
| projection onto invariants | Burnside: `(1/|A|) sum_g prod_D (1+t^{|D|})` | character orthogonality: the `J = ∅` Walsh coefficient |
| critical point | `t = -1`, where every odd-cycle factor vanishes | `q = 1`, the point where every `(q-1)^r` vanishes |
| order of vanishing = minimum rank | `d = min_g c(g)` (Session62) | minimum torus rank (Session59 §2) |
| positivity of the quotient = global free factor | exterior factor `(1+t)^{dim V^A}` | split-torus stratification |
| negative quotient coefficient = obstruction | `P/(1+t)^3 = 1 - 2t + 7t^2 - …` (Session62 (14)) | `(q-1)^2 - 4(q-1) + 7` |
| value / modular jump | `κ = 3000`; `d_ℓ > d` iff `ℓ | κ` (Lemma B0) | `7 = e(U)`; `(q-1) | #U mod 7` iff `7 | 7` |

The same lemma also covers THM-305 (`v_2 A000568(n) = (n-1)/2` is a 2-adic depth of a Burnside
average) and THM-587 (`P_n(-1) = SC(n)`: a critical value that is a count, like `e(U) = 7`). The
method card "Expand at the critical point" is proposed in the session reflection
(`07-reflections/critical-expansion-depth-fifty-vertex-tournament-and-hourglass-twists-klein-20261025.md`);
it is **not** in `META-PATTERNS.md`, because the bounded startup surface is at its byte budget (13 bytes of
slack at the inherited HEAD) and promotion requires a maintainer to compress an existing entry. Its
counterindications are the even action (`1+t+t^2`), complex scalar actions (`1+t^3`), and mixed-order
collisions (`C21`): the method is a reading of an expansion, not a bridge between objects.

The owner's closing sentence is confirmed in this precise sense: the multiplicative root
construction alone determines the twist average and nothing finer; every finer finite-field invariant
of the hourglass (`S_J`, the `1+6` split, the period-12 pair sums) requires the additive branch
equations, i.e. the linear relations among the seven forms.

## 5. Boundaries and non-claims

- Theorem B needs a **regular** `p`-group action; for non-regular actions the cycle types no longer
  depend only on the element order and surviving cycle types of equal cycle count may cancel mod `ℓ`
  (Theorem B′ becomes a lower bound `≥ min{c(τ) : ℓ ∤ e_τ}`).
- The non-abelian corollary `J = 1` is CITED, verified on seven non-abelian groups, not re-proved.
- HYP-9374 (equality in Theorem B″) is open; its failure mode would be a collision of candidate orders
  from different `p'`-parts with cancelling leading coefficients.
- The twist theorem needs odd characteristic and forms with no common zero; it says nothing about
  rational points of any single twist, Picard ranks, or `p_g`.
- No claim about integer magic squares of squares, Hadwiger–Nelson, LRC(14), or A000568 growth.

## 6. Reproduction

```
python3 check_clone_depth.py  > clone-normal.json ; python3 -O check_clone_depth.py  > clone-optimized.json
python3 check_twist_average.py > twist-normal.json ; python3 -O check_twist_average.py > twist-optimized.json
cc -O2 -o census_orbits census_orbits.c ; ./census_orbits z5sq > census_z5sq.out ; ./census_orbits z25 > census_z25.out
sha256sum -c SHA256SUMS
```

Both Python checkers use only the standard library, raise explicit exceptions at every gate (active under
`-O`), and produced byte-identical outputs in both modes. `check_clone_depth.py` imports no automorphism
package and assumes no group: it enumerates automorphisms from arcs. The C census is an independent path
for the two 25-point regular actions. Finite universes: the two tournaments, 22 regular actions, 16
finite fields times three arrangements, 36 seven-subsets.
