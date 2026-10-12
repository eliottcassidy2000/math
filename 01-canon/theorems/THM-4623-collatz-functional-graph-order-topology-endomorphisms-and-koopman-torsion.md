---
id: THM-4623
title: "The Collatz functional graph under order and module lenses: for any self-map T of a set, the T_0 quotient P of the reachability preorder is a dcpo iff its Alexandrov space is sober iff every nonempty directed subset has a greatest element iff T has no divergent orbit (and then the Scott topology is the Alexandrov topology), so the Collatz conjecture says exactly that P has a greatest element (no divergence is the sober half, a unique maximal element the cycle half); the endomorphism monoid of (N, T) has cardinality 2^(aleph_0) although its automorphism group is trivial (THM-4523), by collapsing the ancestor tree of the even preimage of a non-periodic merge vertex m = 4, 16 (mod 18) onto the ancestor tree of its odd preimage; and the Koopman operator f -> f o T on locally constant functions on Z_2 has (t - 1)-primary torsion of height one only (Ker rho(theta) = K.1 if rho(1) = 0, else 0), while its transpose at level L, the Terras/de Bruijn transition matrix on Z/2^L, is pure torsion with P_L^L = J/2^L"
status: >
  PROVED (elementary), WITH CORRECTIONS applied after independent audit B (2026-10-12; MISTAKE-594):
  statement 2 now requires m = 4 or 16 (mod 18) and m non-periodic (for m = 10 (mod 18) the odd
  preimage is divisible by 3 and the collapsing map is unique, so the continuum count was FALSE as
  first stated); the transpose of the Koopman operator was misidentified with THM-4520's 3-adic
  level operator (a different object) and is now the Terras/de Bruijn matrix P_L with P_L^L = J/2^L
  (Terras's parity-vector bijection; minimal polynomial divides t^L (t - 1), equal to it for
  L <= 8 by the audit's computation); directed sets are nonempty; 'dcpo with a top' is simplified to
  'has a greatest element'. The Koopman statement uses the KNOWN topological conjugacy of T on Z_2
  with the one-sided shift (Lagarias 1985; Bernstein 1994) only through two consequences: periodic
  points of every sufficiently large period are dense in every residue class (checked exactly by the
  audit for L <= 6), and some orbit is dense. The endomorphism construction was checked by the audit
  at all n <= 2*10^5 for m = 16, 22. Nothing here bears on the Collatz conjecture itself; the
  dictionary is a reformulation. Inspired by Xu-Ji (arXiv:2609.19612), Pinsky (arXiv:2305.11382)
  and Chini (arXiv:2609.20184).
session: mac-mini-2026-10-12-logic4
source: 05-knowledge/results/logic_lenses_jc2_collatz_hn_20261012.md
related:
  - THM-4523 (the Collatz functional graph is rigid: Aut = 1)
  - THM-4520 (level operators: Gauss-twisted circulants, spectrum on the half circle)
  - THM-4581 (Haar coalescence of affinely related orbits)
  - THM-4620 (the same t-primary kernel structure for Keller dual fields)
---

# THM-4623 — order topology, endomorphisms and Koopman torsion of the Collatz graph

## Setting

`T` is a self-map of a set `X`; `x ≼ y` iff `T^k x = y` for some `k ≥ 0` (reachability, a preorder);
`x ~ y` iff `x ≼ y ≼ x` (same cycle, or equal); `P = X/~` with the induced partial order. The
Alexandrov topology `Alex(P)` has the up-sets as open sets, so closed sets are down-sets and the
closure of a point is `↓[x]`. A *divergent* orbit is one that is not eventually periodic. For the
Collatz map `T(n) = n/2` (even), `3n + 1` (odd) on `N = {1, 2, ...}`, `D(n)` is the number of steps to
reach the cycle `{1, 2, 4}` (when it is reached).

## Statement

1. **Order dictionary.** The following are equivalent: (i) `T` has no divergent orbit; (ii) `P` is a
   dcpo; (iii) `Alex(P)` is sober; (iv) every nonempty directed subset of `P` has a greatest element.
   When they hold the Scott topology of `P` equals `Alex(P)`. Consequently, for the Collatz map on `N`:

       Collatz conjecture  ⟺  P has a greatest element  ⟺  Alex(P) is sober and irreducible,

   and the honest split of the conjecture is: no divergence ⟺ (i)–(iv); given that, no other cycle ⟺
   `P` has a unique maximal element. (A greatest element alone already excludes divergence, so the
   dcpo clause is redundant in the conjunction.)

2. **Endomorphisms.** Let `m ∈ N` be a non-periodic merge vertex with `m ≡ 4` or `16 (mod 18)` (so
   `m ≡ 4 (mod 6)` and the odd preimage `v = (m-1)/3` is prime to `3`; e.g. `m = 16, 22, 34, 40`), with
   preimages `u = 2m` and `v`, and let `anc(u)` be the set of ancestors of `u` (vertices `n` with `T^d n = u`, `d ≥ 1`), with
   `depth(n) = d`. For every map `h : anc(u) ∪ {u} → N` with `h(u) = v` and `h(T n) = T(h n)` on `anc(u)`,
   the map `φ_h` equal to `h` on `anc(u) ∪ {u}` and to the identity elsewhere is an endomorphism of
   `(N, T)` with `φ_h(u) = φ_h(v)`. There are `2^(ℵ_0)` such `h` (for `m ≡ 10 (mod 18)`, where `3 | v`,
   the ancestor tree of `v` is the bare doubling chain and `h` is unique — audit B); so
   `|End(N, T)| = 2^(ℵ_0)` while `Aut(N, T) = 1` (THM-4523). Every endomorphism `φ` maps the cycle
   `{1, 2, 4}` to itself by a rotation (the only solutions of `T^3 x = x` in `N` are `1, 2, 4`) and
   satisfies `D(φ(n)) ≤ D(n)` wherever `D(n)` is defined.
3. **Koopman torsion.** Let `K` be a field of characteristic zero, `LC(Z_2, K)` the locally constant
   functions, and `θ f = f ∘ T` (with `T` the 2-adic Collatz map `T(x) = x/2`, `(3x+1)/2` on `Z_2`). For
   every nonzero `ρ ∈ K[t]`:

       Ker ρ(θ) = K · 1  if ρ(1) = 0,      Ker ρ(θ) = 0  otherwise.

   In particular `θ - λ` is injective for every `λ ≠ 1`, `Ker(θ - 1) = K·1`, and `Ker(θ - 1)^2 = Ker(θ - 1)`.
   The transpose statement: restricted to level `L` (functions on `Z/2^L`, averaged over the two
   lifts to `Z/2^(L+1)`), the transpose of `θ` is the Terras/de Bruijn transition matrix
   `P_L(a → b) = ½ #{x ∈ {a, a + 2^L} : T x ≡ b (mod 2^L)}`, and `P_L^L = J/2^L` (every entry `2^(-L)`),
   so `P_L` is pure torsion with minimal polynomial dividing `t^L (t - 1)` (equal to it for `L ≤ 8`):
   eigenvalue `0` besides the stationary `1`. (THM-4520's level operator lives on the units modulo
   `3^n` and is a different object.)

## Proofs

**(1).** *(i)⇒(iv).* Let `S ⊆ P` be directed (nonempty). Each class `[x] ∈ S` reaches a unique cycle `C_x`; a common
upper bound of `[x], [y]` in `S` is reachable from both, so `C_x = C_y`: all of `S` reaches one cycle
`C`. If `[C] ∈ S` it is the greatest element (only `C` is reachable from `C`). Otherwise every class in
`S` is a single pre-periodic point, and for `x ∈ S` the set `{k : T^k x ∈ S}` is finite (the orbit enters
`C ∉ S`); let `z_x := T^(k(x)) x` with `k(x)` maximal. For `x, y ∈ S` an upper bound `w ∈ S` of `z_x, z_y`
has `w = T^i z_x` with `T^(k(x)+i) x ∈ S`, so `i = 0` and `w = z_x`; likewise `w = z_y`. So all `z_x` are one
point `z ∈ S`, which is the greatest element.
*(iv)⇒(ii):* a greatest element is a supremum.
*(ii)⇒(i):* a divergent orbit `x, Tx, T^2 x, ...` has pairwise distinct points and is a chain in `P`;
an upper bound `[w]` would satisfy `w = T^(k + m_k) x` for every `k`, with `k + m_k` constant by
injectivity of the orbit — impossible. So the chain has no supremum.
*(i)⇒(iii):* in an Alexandrov space a closed set `D` is irreducible iff it is a directed down-set
(if directed, two proper closed subsets miss points `x, y ∈ D` and then miss a common upper bound; if
`x, y ∈ D` have no common upper bound in `D`, then `D = (D ∖ ↑x) ∪ (D ∖ ↑y)` is reducible). By (iv) a
directed down-set `D` has a greatest element `z`, so `D = ↓z` is the closure of `[z]`, unique by `T_0`.
*(iii)⇒(i):* for a divergent orbit, `D = ∪_k ↓[T^k x]` is a directed down-set; if `D = ↓[w]` then
`w ≼ T^k x` for some `k` and `T^(k+1) x ≼ w`, so `T^k x ~ T^(k+1) x`, making the orbit periodic.
*Scott = Alexandrov under (i):* for an up-set `U` and a directed `S` with `sup S ∈ U`, the greatest
element of `S` is `sup S ∈ U ∩ S`; so every up-set is Scott-open.
*Collatz:* the conjecture says every orbit reaches `{1, 2, 4}`, i.e. `[{1, 2, 4}]` is above every
class. Conversely a greatest element `[c]` is reachable from every point, so no orbit diverges and
`c` lies on a cycle reached by everything; since `1` reaches no other cycle, it is `{1, 2, 4}`.
A sober irreducible Alexandrov space is one whose whole space is a directed down-set with a greatest
element. The split: (i)–(iv) are 'no divergence'; given no divergence every orbit ends in a cycle,
the cycles are the maximal elements of `P`, and 'one cycle' is 'one maximal element'.

**(2).** For `n ∈ anc(u)`, `T n ∈ anc(u) ∪ {u}` and `φ_h(T n) = h(T n) = T(h n) = T(φ_h n)`. For `n = u`:
`T u = m ∉ anc(u) ∪ {u}` because `m` is not periodic, so `φ_h(T u) = m = T(v) = T(φ_h u)`. For
`n ∉ anc(u) ∪ {u}`, also `T n ∉ anc(u) ∪ {u}`, so `φ_h` is the identity on both. Thus `φ_h ∘ T = T ∘ φ_h`.
(`anc(u)` and `anc(v) ∪ {v}` are disjoint: a common ancestor would visit `m` twice.) Existence and
cardinality: define `h` along the doubling chain `u, 2u, 4u, ...` by a target sequence `t_0 = v`,
`t_(k+1) ∈ T^(-1)(t_k)`, and on the rest of `anc(u)` by `h(n) = 2^(depth(n) - k) t_k` when `n` hangs
off `2^k u` (then `h(T n) = T(h n)` holds at every vertex: along the chain by construction, elsewhere
because `T(2^j t_k) = 2^(j-1) t_k`). At every `k` with `t_k ≡ 4 (mod 6)` there are two admissible
`t_(k+1)` (`2t_k` and `(t_k - 1)/3`); choosing the odd one only when `(t_k - 1)/3 ≢ 0 (mod 3)` keeps
infinitely many future options, provided the start `t_0 = v` is itself prime to `3` (this is where
`m ≢ 10 (mod 18)` is needed: along the doubling chain of an odd `t ≢ 0 (mod 3)`, the odd preimages
`(2^k t - 1)/3` run through all residues modulo `3`, since `(2^(k+2) t - 1)/3 = 4 (2^k t - 1)/3 + 1`,
whereas the doubling chain of a multiple of `3` has no odd preimages at all).
So every infinite binary sequence of choices gives a distinct `h` (they differ at the first
divergence), hence `2^(ℵ_0)` endomorphisms; the upper bound `|N^N| = 2^(ℵ_0)` is trivial. For a
general endomorphism `φ`: `T^3 φ(c) = φ(T^3 c) = φ(c)` for `c` in the cycle, so `φ(c)` solves `T^3 x = x`,
whose only solutions in `N` are `1, 2, 4` (the parity patterns `(odd, even, even)`, `(even, odd, even)`,
`(even, even, odd)` give `x = 1, 2, 4`; two consecutive odd steps are impossible and three even steps
force `x = 0`); so `φ` restricted to the cycle commutes with the rotation;
`T^(D(n)) φ(n) = φ(T^(D(n)) n)` lies in the cycle, so `D(φ(n)) ≤ D(n)`.

**(3).** Kernels commute with the base change `K ⊆ K̄`, so assume `K = K̄` and, by the Chinese remainder
theorem, `ρ = (t - λ)^k`. *Eigenvalue `λ ≠ 1`.* Let `f ∘ T = λ f`. At a point `p` with `T^k p = p`,
`f(p) = λ^k f(p)`. If `λ` is not a root of unity, `f` vanishes at every periodic point; if `λ` is a
primitive `m`-th root of unity (`m ≥ 2`), `f` vanishes at every point of period `k` with `m ∤ k`. Under
the conjugacy with the one-sided shift, the residue classes modulo `2^L` are the cylinders of length
`L`, and every cylinder contains points with `T^k p = p` for every `k ≥ L` (periodically repeated
words), in particular for `k` prime to `m`. So `f` vanishes on a dense set; a locally constant function
has clopen fibres, so `f = 0`. Hence `Ker (θ - λ)^k = 0` by induction on `k`. *Eigenvalue `1`.* If
`f ∘ T = f` and `x_0` has a dense orbit, then `f^(-1)(f(x_0))` is a clopen set containing a dense set, so
`f` is constant. If `(θ - 1) f = κ·1` then `f(T^n x) = f(x) + nκ`, and a locally constant function on the
compact space `Z_2` takes finitely many values, so `κ = 0`: `Ker (θ - 1)^2 = Ker (θ - 1) = K·1`, and by
induction `Ker (θ - 1)^k = K·1`. The transpose remark: `P_L` is the matrix of `θ` restricted to level
`L` and averaged back; by Terras's theorem the parity vector of length `2L` is a bijection
`Z/2^(2L) → {0,1}^(2L)`, so for each class `a mod 2^L` the `2^L` lifts `x ≡ a (mod 2^L)`, `x mod 2^(2L)`,
realize every length-`L` parity vector of `T^L x` exactly once, i.e. every class `T^L x mod 2^L` exactly
once: `P_L^L(a → b) = 2^(-L)` for all `a, b`. Hence `P_L^L (P_L - I) = 0`.  ∎

## Reading

* **Sobriety is the divergence question.** Xu–Ji prove Scott sobriety for countable meet-continuous
  lattices by successive finite choices on a countable set; for the orbit space of a map the sober /
  dcpo property is exactly "no divergent orbit", and the Collatz conjecture adds "one cycle". No
  lattice theorem can decide it, but the dictionary places the two halves of the conjecture as
  (sober) and (unique maximal element); the single predicate "`P` has a greatest element" is the
  whole conjecture.
* **Frucht in reverse.** THM-4523's rigidity (a trivial automorphism group) coexists with a huge
  endomorphism monoid: the Collatz graph realizes the trivial group, yet it is far from "cored".
  Endomorphisms cannot raise the hitting time `D`, so they are not a route to divergence.
* **Chini's identities for the Koopman operator.** The kernel lattice of `θ` on `LC(Z_2, K)` is as
  small as possible: `(t - 1)`-primary of height one; its finite-level transposes are nilpotent plus
  rank one. THM-4620 finds the analogous `t`-primary
  structure for Keller dual fields. Both are "almost generic" endomorphisms in the sense that no
  finite kernel identity beyond the forced one holds.
