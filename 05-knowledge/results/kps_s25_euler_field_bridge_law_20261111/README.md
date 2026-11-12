# Euler fields in twisted de Rham classes; block form of the root-change grammar (kind-pasteur-2026-11-11-S25)

Owner directive (2026-11-11): continue the session on Stein modules of Keller pairs (S24, THM-5138, HYP-9668) and
extend it into a session that aims to prove the Chan-Pak conjecture.

Canon: THM-5147 (part A), THM-5148 (part B), HYP-9677 (bounded-run root-change words); HYP-9668 updated ((E')
REFUTED). Pre-registration: `PREREGISTRATION.md` (commits 7479352f36, 870ad0c5d7, d67741ce21, each pushed before
its runs; eleven predictions: C1, C2, D1, D3, E1 CONFIRMED; C5 CONFIRMED for tuples only; C3, C4, D2, E2, E3
FAILED; D2 and D3 were registered against a fact already in canon, see B4).
**JC(2) and the Chan-Pak conjecture remain OPEN. No proof of Chan-Pak was obtained.**
Audit: `audit/AUDIT.md`, independent, own code: **PASS WITH CORRECTIONS** (five major findings, sixteen minor,
eight upgrades; all applied in this version; MISTAKE-768). The largest correction: the tournament family of part B
is the O/Z/J root-change grammar already in canon, and its list of 29 exceptions is KNOWN to `10^9` (HYP-9494 (D),
a packet of this same agent). Part B is therefore a credited re-derivation plus two new finite measurements.

## 0. What comes out

- **Part A (Stein thread; PROVED).** An elementary identity for twisted de Rham classes: if a vector field `xi`
  satisfies `xi(W) = W` and has constant divergence `delta`, then `W h = (xi + delta) h` modulo the image. Three
  consequences: (1) **HYP-9668 (E') is false**: for Broughton's polynomial `W = x + x^2 y` the image has
  codimension one, does not contain 1, and is not a Mathieu-Zhao subspace; for `W = x + x^(l+1) y`: `W^m` lies in
  the image for `m >= l` and `x^l W^n = n! x^l` never does. The audit extended this, without any Euler field, to
  every `W = x^(l+1) y + r(x)` (section A2, U1-U2). (2) A radial identity for the Stein module of a Keller pair:
  `[(UV)^m h] = [(E + 1)(E + 2)...(E + m) h]`, `E = (U d_u + V d_v)/2`. (3) A one-line proof of the charge grading
  used in S24 (Proposition 6.3 there).
- **Part B (Chan-Pak; no proof).**
  - *Re-derived here in block form (KNOWN: SEAM.md (1)-(7), THM-4875 (B), THM-4904 Setting, THM-4888 (D) 1).* With
    the root moved to the new apex, joining a one-apex block `w` is a 3 x 3 integer matrix `M_w`; the path counts of
    iterated joins of one-apex tournaments at their apexes are the first coordinates of the orbit of `(1, 0, 0)`.
    This family is the O/Z/J root-change grammar of Sessions 197-200. New: only the closed entries of a block.
  - *Re-computed (KNOWN: Session197; THM-4884 (2b); HYP-9494 (D), to `10^9`).* The family misses exactly 29 odd
    numbers, the largest 2647 (here to `2^18`).
  - *New, FINITE-EXACT, pre-registered.* (i) The sub-language with a transitive chain inserted at every join (every
    root change preceded by a sink) is far from complete: 1214 odd numbers missed below `2^14`, 5.26 per cent of
    `[2^17, 2^18)` at depth 3, with 14 distinct chains per odd number. (ii) The sub-language with at most `L` spine
    letters per root change: `L = 3, 4, 5` miss 1355, 127, 49 odd numbers below `2^17`; for `L = 5` (63 blocks)
    there is no miss in `(10391, 2^19]`. HYP-9677 states (ii) as a conjecture; it is a variant, in the uniquely
    Hamiltonian grammar, of HYP-9640 (B) (56 two-port block matrices of order at most 4, every odd `v >= 493`).
  - Five of the eleven predictions failed; each time I predicted more coverage than a restricted family has.

## Part A. Euler fields

### A1. The lemma

Setting: `B = C[x, y]`, `W` in `B`, `M_W = Im(d/dx - W_x) + Im(d/dy - W_y)`, classes `[h]` in `B / M_W`
(the top twisted de Rham cohomology `H^2(C^2, d - dW ^)`). The same holds on a smooth affine surface with a
nowhere vanishing 2-form, in particular in the coordinates `(u, v)` of an etale map (S24 note, section 2).

**Lemma A1 (PROVED; elementary).** For every polynomial vector field `v`, `div(v) - v(W)` lies in `M_W`, and
`M_W` is the set of these elements. Hence for a vector field `xi` and `h` in `B`:
`xi(h) + (div xi - xi(W)) h` lies in `M_W`. In particular:

- (Euler) if `xi(W) = W` and `div xi = delta` is constant, then `[W h] = [(xi + delta) h]`, and for all `m >= 1`
  `[W^m h] = [(xi + delta)(xi + delta + 1) ... (xi + delta + m - 1) h]`; for an eigen-polynomial (`xi h = e h`):
  `[W^m h] = (e + delta)(e + delta + 1)...(e + delta + m - 1) [h]`;
- (symmetry) if `eta(W) = 0` and `div eta = 0`, then `eta(B)` lies in `M_W`.

*Proof.* For `v = p d/dx + q d/dy`: `(d/dx - W_x) p + (d/dy - W_y) q = div(v) - v(W)`. Take `v = h xi`:
`div(h xi) = xi(h) + h div xi`. For the powers: `xi(W^(m-1) h) = W^(m-1) (xi + m - 1) h`, and induction. QED

Known cases: for `W = UV` and the rotation field `eta = U d_u - V d_v` the symmetry statement is S24's
Proposition 2.2 (e); for `W = xy` and `xi = (x d/dx + y d/dy)/2`, `delta = 1`, the Euler statement gives
`[(xy)^m] = m! [1]`. The lemma itself is the standard description of the image as "divergence minus derivative
along `W`"; no novelty is claimed for it.

### A2. Broughton's polynomial: HYP-9668 (E') is false

**Theorem A2 (PROVED).** Let `l >= 1` and `W = x + x^(l+1) y`. Then:

1. the linear form `lambda(x^k y^j) = (-1)^n j! / n!` for `n = (l+1)(j+1) - 1 - k >= 0`, and 0 for `n < 0`,
   vanishes on `M_W`; `lambda(1) = (-1)^l / l!` and `lambda(x^l) = 1`. So 1 is not in `M_W`;
2. `W^m` lies in `M_W` for every `m >= l`, and `x^l W^n = n! x^l` modulo `M_W` for every `n >= 0`;
3. so `P = W^l` has all its powers in `M_W` while `x^l P^m` never lies in `M_W`: **`M_W` is not a Mathieu-Zhao
   subspace of `C[x, y]`**, and `M_W = ker lambda` has codimension one.

For `l = 1` this is Broughton's polynomial `x + x^2 y` (generic fibre `C^*`, no critical point, `[1] != 0`): the
test case of HYP-9668 (E').

*Proof.* (1) `lambda(h) = Res_(x=0) int_0^infinity h e^(-W) dy dx`: for `h = x^k y^j` the inner integral along
the ray where `x^(l+1) y > 0` is `j! x^(k - (l+1)(j+1)) e^(-x)`, and the residue is the coefficient of
`x^(-1)`. It vanishes on `(d/dy - W_y) B` because the inner integral of a `y`-derivative is
`-h(x, 0) e^(-x)`, which has no residue, and on `(d/dx - W_x) B` because the residue of an `x`-derivative is zero.
Term by term, for every `l` (audit, section 4.1; 672 generators checked): with `N = (l+1)(j+1) - k`,
`lambda((d/dx - W_x) x^k y^j) = (-1)^N (j!/N!) (k + N - (l+1)(j+1)) = 0`, and in `lambda((d/dy - W_y) x^k y^j)`
the two terms have the same index and cancel.
(2) `xi = x d/dx - l y d/dy` satisfies `xi(x) = x`, `xi(x^(l+1) y) = x^(l+1) y`, so `xi(W) = W`, and
`div xi = 1 - l`. By Lemma A1, `[W^m] = (1 - l)(2 - l)...(m - l) [1]`, which is zero for `m >= l`; and
`xi(x^l) = l x^l` gives `[W^n x^l] = (l + 1 - l)(l + 2 - l)...(n) [x^l] = n! [x^l]`. (3) `lambda(x^l) = 1`;
codimension one is Lemma U1 below. QED

**Lemma U1 (audit; PROVED, elementary; re-checked by hand).** Let `W = x^(l+1) y + r(x)`, `l >= 1`, `r` any
polynomial. Then `B = C x^l + M_W`. *Proof.* `(d/dy - W_y)(x^k y^j)` gives `x^(k+l+1) y^j = j x^k y^(j-1)` modulo
`M_W`, so `B/M_W` is spanned by the `x^k y^j` with `k <= l`; let `V_j` be the span of those of `y`-degree at most
`j`. For `1 <= k <= l`, `(d/dx - W_x)(x^k y^j)` and the previous relation give
`(k - (l+1)(j+1)) x^(k-1) y^j = r'(x) x^k y^j`, with a nonzero coefficient; by descending induction on `k`,
`x^(k-1) y^j` lies in `C x^l y^j + V_(j-1)`. For `k = 0`: `x^l y^(j+1) = -r'(x) y^j / (l+1)` lies in `V_j`. So
`V_j = V_0 = C x^l` modulo `M_W`. QED (The first draft quoted a comparison theorem and S24's machine control for
the codimension; neither is needed.)

**Theorem U2 (audit; PROVED; re-checked by hand).** Let `W = x^(l+1) y + r(x)`, `l >= 1`, `r(0) = 0`. The form
`lambda_r(x^k y^j) = j! [x^((l+1)(j+1)-1-k)] e^(-r(x))` vanishes on `M_W` (same proof as A2 (1)),
`lambda_r(x^l) = 1`, so `M_W = ker lambda_r`; and `lambda_r(x^k W^m) = [x^(l-k)] (m! e^(-r) e_m(r))` with `e_m` the
truncated exponential series. Since `m! e^(-s) e_m(s) = m! + O(s^(m+1))`: `W^m` lies in `M_W` for `m >= l` and
`x^l W^n = n! x^l` modulo `M_W`. **So `M_W` is not a Mathieu-Zhao subspace for any `W = x^(l+1) y + r(x)`**, and
1 lies in `M_W` iff `[x^l] e^(-r) = 0`. No Euler field is needed. (`M_W` depends on `dW` only, so `r(0) = 0` is no
restriction.) Exact certificates for five potentials: `audit/a3`, `audit/a5`.

Machine (`scripts/euler_field.py`, `runs/euler_field.txt`; the S24 engine, exact modulo a prime, membership
relative to the truncation): for `l = 1`: codimension 1, `[1] != 0`, `W h = (k - j) h` on 20 monomials,
`W^m` in `M` and `x W^m - m! x` in `M` for `m <= 4`, `x W^m` not in `M`; for `l = 2`: codimension 1, `W + 1`,
`W^2`, `W^4` in `M`, `x^2 W^n - n! x^2` in `M` for `n <= 4`, `[x^2] != 0`. The audit found 23 explicit vector
fields over `Q` certifying the memberships for `l = 1, 2, 3` (`audit/a2`).

Literature (audit, section 4.2: one search, three papers read through summaries; a search, not a clearance).
Theorem A2 is an instance of the statement of Zhao's Conjecture 1.3 (arXiv:0902.0210: images of commuting
operators `h(z) + sum c_i d/dz_i`), not of the Image Conjecture as restated by van den Essen, Wright and Zhao
(arXiv:1008.3962, with a regular sequence), whose special case SIC is what THM-2801 refutes. The general statement
is therefore already REFUTED in canon; A2 is not a first counterexample. What was not found in the sources read: a
failure with 1 outside the image for a pair `(d/dx - W_x, d/dy - W_y)` over a field of characteristic zero (the
nearest: arXiv:1004.0521, Example 2.4, a derivation, by the 1-property). No priority is claimed.

**Corollary A3 (a second obstruction to the Mathieu property; PROVED).** Let `W` have a polynomial vector field
`xi` with `xi(W) = W` and `div xi = delta` a non-positive integer. If some eigen-polynomial `h` of `xi`
(`xi h = e h`) has `[h] != 0` and `e + delta` not a non-positive integer, then `M_W` is not a Mathieu-Zhao subspace:
`P = W^(1 - delta)` has all powers in `M_W` and `h P^m` never lies in it. (The eigen-polynomial is not essential:
any `h` with `[(xi + delta)(xi + delta + 1)...(xi + delta + N - 1) h] != 0` for infinitely many `N` works; in
applications `[h] != 0` must come from a linear form, as `lambda` above, not from a truncated normal form.)

So for the question (Q) of HYP-9668 ("for which potentials is the image a Mathieu-Zhao subspace") there are two
necessary conditions: the 1-property (1 is not in the image unless the image is everything; S24, `x^p y^q`,
`p != q`) and the Euler obstruction, which can fail with 1 outside the image. For the Gaussian potential
`delta = 1`, and for `W = (f - a)(g - b)` of a Keller pair the field `E = (U d_u + V d_v)/2` has `delta = 1` as
well: the obstruction is absent there. (It is absent for every potential without an Euler field of non-positive
integer divergence, for instance `x^3 + y^3`, where `delta = 2/3`.)

What this does not settle. Every polynomial without critical point has fields with `xi(W) = W`; what can be
missing is one of constant divergence. Example (audit, U3): `W = x (x^2 y + x + 1)` has none (it would give
`[W] = delta [1]`, `[W^2] = delta (delta + 1) [1]`, but `W - 1` and `W^2` lie in `M_W` and `[1] != 0`), and it is
nevertheless settled by Theorem U2. Not settled: rank-one potentials outside the monomials `x^p y^q` (1-property)
and the family `x^(l+1) y + r(x)` (U2). The statement "among rank-one potentials only `xy`, up to automorphisms,
has a Mathieu-Zhao image" is consistent with everything computed and is not proved.

### A3. The Stein module of a Keller pair

**Proposition A4 (radial identity; PROVED).** For an etale map `(u, v) : X -> C^2`, a centre `z`, `U = u - a`,
`V = v - b`, `E = (U d_u + V d_v)/2`, and every `h` in `O(X)`:

    [(UV)^m h]_z = [(E + 1)(E + 2) ... (E + m) h]_z,      and      [(U d_u - V d_v) h]_z = 0.

*Proof.* Lemma A1 in the frame `du ^ dv`, where `div E = 1`, `E(UV) = UV`, and the rotation field has divergence
0 and kills `UV`. QED

Checks: 7 of 7 on the Keller pair `(x + y^2, y)` at the centre (2, 3); 6 of 6 on the slice map of S24, section 6.2;
audit: 64 + 126 identities with explicit functionals.

Uses. (i) The grading of the Stein kernel by charge at the special centre of S24's Proposition 6.3 (the plane minus
a hyperbola, `W = xy`): `eta = x d/dx - y d/dy` preserves `1 + xy`, kills `W` and has divergence 0, so every
element of nonzero charge is `eta` of something and lies in the kernel. (ii) `[(UV - c)^m]_z = m! e_m(-c) [1]_z`
with `e_m` the truncated exponential series. Where `[1]_z != 0` (every centre for an automorphism; every centre
algebraically independent over a field of definition, S24 Theorem 5.4), `UV - c` is therefore not a null element
for any `c` (`e_1(-c) = 1 - c`, `e_2(-1) = 1/2`). Where `[1]_z = 0` (the double cover of S24 at the centres
`(0, b)`) every `UV - c` is null. (The first draft said "never a null element" without the hypothesis: audit M5.)
(iii) It is the radial half of the structure whose angular half is the rotation field; in the proof of GMC(2)
these are the two layers (THM-1645).

## Part B. Block form of the root-change grammar

### B1. What is known (credits)

- Odd-cycle formula `H = I(Omega, 2)` (Grinberg-Stanley; canon) and its increments: THM-4876
  (`H(T) - H(T - v) = 2 sum over odd cycles C through v of H(T - V(C))`, KNOWN).
- The join of two rooted tournaments (the ordinal sum `X => Y` with the arc between the two roots reversed):
  `H = x_h y_h + x_R y_R - x_e y_e` in rooted triples (THM-4866 (B)); Session263, Theorem RJ2 (return-join
  recursion). In the notation below `x_h = H(X)`, `x_R = W`, `x_e = o - e` (= the number of Hamiltonian paths
  ending at the root; audit, 202,013 rooted tournaments).
- **The universal seam transfer** (`05-knowledge/results/seven_return_native_recursions_20261029/SEAM.md`, (1)-(7);
  THM-4875 (B)): every block `Y` with entry and exit ports, coincident ports included, has an integer 3 x 3 matrix
  `M_Y` with `v(X.Y) = M_Y v(X)`, `M_(X.Y) = M_Y M_X`. Gadget transport by explicit integer matrices and rays from
  inserted transitive chains: THM-5126 (A), THM-5130 (A) (this agent, S21, S22); port chains: THM-5070 (S18).
- **The O/Z/J root-change grammar** (Sessions 197-200; THM-4904, Setting): letters act on `(D, E, U, d)` from the
  rooted edge `(1, 0, 2, 1)`; `O`, `Z` add a vertex above / below the root, `J` adds a vertex and moves the root
  to it. `JZ = OJ` (Session205; THM-4888 (D) 1). Strong members = the uniquely Hamiltonian strong tournaments
  (Session197). **Its value set misses exactly 29 odd numbers**: Session197 `C1_COMPLETE.md` (to 5,001), THM-4884
  (2b) (to 50,001), HYP-9494 (D) (to `10^9`, with the root-change depth of every value; this agent, S14). Two
  root changes miss exactly 140 odd values below `10^9`, the largest 3,780,967 (HYP-9494 (D)). Multiplicity by
  number of root changes: HYP-9480 (D).
- Finitely generated versions: HYP-9640 (B) (this agent, S18): the 56 two-port block matrices of order at most 4
  reach every odd `v >= 493` (to `2 * 10^5`), a 13-matrix subset every odd `v >= 961`.
- One-apex tournaments (one vertex whose removal leaves a transitive tournament): `H = g(v)`, `g(0) = 1`,
  `g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`, rooted triple `(g(v), e(v), 1 + v)` (THM-4875, THM-4885,
  THM-4890). One-apex x one-apex joins `(1+u)(1+v) + g(u) g(v) - e(u) e(v)`: mass `(pi^2/24) X ln X` (THM-4885
  (B)). Joins of a bank of 2,732 rooted triples of order at most 7 with one-apex tournaments reach every odd number
  up to `2^35` except 7 and 21 (THM-4885 (D), FINITE-EXACT; `2^30` audited, `2^35` one run not independently
  recomputed); the mean multiplicity of a fixed bank is a constant (THM-4885 (C), `M_Y`).
- CONCURRENT today: opus-2026-11-11-S93 (THM-5141: counting exponent 3/4 for one-apex counts; its registered W1:
  two-apex words of length at most 16 reach every odd number below `2^23` except 7 and 21), opus-2026-11-11-S94
  (every odd `N < 226,079,703` except 7, 21 is a two-apex value of length at most 18), Sessions 275 and 276 (digit
  banks, carries, corridors). Positive density of the spectrum is OPEN.

I re-derived the join law, the block transfer and the family below from the odd-cycle formula before my sweep
returned, and the sweep missed my own S13, S14 and S18 packets (a grep for "7, 13, 21, 41, 43, 55" finds THM-4884
at once). None of the laws in B2 and B4 is a result of this session. See MISTAKE-768.

### B2. The chain law (KNOWN, in block form)

For a tournament `X` and a vertex `u` put `W^+_X(u) = sum over directed paths P starting at u of H(X - P)`,
`W^-_X(u)` the same with paths ending at `u`, and `tau(Z) = sum over all directed paths P of Z, the empty one
included, of H(Z - P)` (`H` of the empty tournament is 1).

**Proposition B1 (KNOWN, re-derived: THM-4866 (B) applied to `(X => C, u)`, whose triple is
`(H(X), 2^i W^+_X(u), 0)`, plus the third row of the seam matrix of `C => Y`).** Let `(X, u)`, `(Y, v)` be rooted
tournaments and `T = X => C => Y` the ordinal sum with a transitive chain `C` of `i >= 1` vertices in the middle,
with the arc between `u` and `v` reversed (`v -> u`). Then

    H(T) = H(X) H(Y) + 2^i W^+_X(u) W^-_Y(v),          W^+_T(v) = H(X) W^+_Y(v) + 2^i W^+_X(u) tau(Y - v).

*Proof (audit, U4; no cycle formula).* A Hamiltonian path of `T` either avoids `v -> u`, and is a path of
`X => C => Y` (`H(X) H(Y)` of them), or is `(A_X)(A_C)(A_Y -> v) -> (u -> B_X)(B_C)(B_Y)` with `B_X` a path of `X`
from `u`, `A_X` a Hamiltonian path of the rest of `X`, `A_C + B_C = C` (`2^i` splittings), `A_Y` a path of `Y`
ending at `v`, `B_Y` a Hamiltonian path of the rest. A path of `T` starting at `v` either stays in `Y`, with
complement `X => C => (Y - P)`, or is `v -> u`, a path of `X` from `u`, a subset of `C`, a path of `Y - v` or
nothing. QED Brute force: 87,825 joins (audit).

For a one-apex tournament rooted at its apex, with apex word of value `p` **without leading letter 0** (letters:
1 if the apex beats the spine vertex; first spine vertex = top bit), `L(p)` the bit length and `R(p)` the reversal
of the complement of `p` on `L(p)` bits: `H = g(p)`, `W^+ = p + 1`, `W^- = R(p) + 1`, `tau = 2^(L(p))`. (With `k`
leading letters 0: `W^- = 2^k (R(p) + 1)`, `tau = 2^(k + L(p))`.) So:

**Corollary B2 (block form of the KNOWN letter moves; brute force 1470 of 1470 chains of depth 2 and 3,
`runs/c1_bruteforce.txt`; audit: 400 random chains).** With

    M(p, i) = [[ g(p), 2^i (R(p) + 1) ], [ p + 1, 2^(i + L(p)) ]]        (p >= 0, i >= 1),

the chain `Y_0 => C_1 => Y_1 => ... => C_k => Y_k` of one-apex tournaments with apex words `p_0, ..., p_k`, chains
of `i_1, ..., i_k >= 1` vertices and the arcs between consecutive apexes reversed has
`H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1)`. Leading letters 0 of a word are absorbed in the
exponent `i`, and trailing letters 1 of a block are chain vertices of the next join:
`M(p', i') M(2p + 1, i) = M(p', i' + 1) M(p, i)`; both are Session205's `JZ = OJ`. In O/Z/J terms these *pure
chains* are the words in which every `J` is immediately preceded by `O` (every root change is preceded by a sink):
the closure of that sub-language computed from the letter moves gives the same miss list (audit, `d4`).

### B3. The pure chains: pre-registered test (FINITE-EXACT; new as far as one search shows)

`scripts/chain_semigroup.py`; outputs in `runs/`. The audit recomputed every number with independent programs
without normal forms (`audit/b3`, `b4`, `b5`).

| prediction | outcome |
|---|---|
| C1: formula agrees with brute force | **CONFIRMED**: 1470 cases, 0 mismatches |
| C2: depth at most 2 misses more than 2 per cent of the odd numbers of `[2^17, 2^18)` | **CONFIRMED**: 11.06 per cent (7250 of 65536) |
| C3: depth at most 3 misses fewer than 0.5 per cent there, fewer than 600 below `2^18` | **FAILED**: 5.26 per cent (3446); 9190 below `2^18` |
| C4: full closure below `2^14` misses only 7 and 21 | **FAILED**: 1214 odd numbers missed (28,155 live states); 414 below `2^12` |
| C5: at least 40 tuples per odd number at depth at most 3 | **CONFIRMED for tuples only** (audit M3): 43.9 tuples of the registered matrices; 82.9 in the program's parametrisation, which counts a chain with chain lengths `c_1, c_2` exactly `c_1 c_2 (c_2 + 1)/2` times. Distinct chains per odd number: **14.3** (671,731 of depth 2 and 1,199,555 of depth 3 below `2^18`), below the registered 40 |

First misses of the full closure: 7, 13, 21, 29, 35, 41, 43, 49, 55, 59, 67, 71, 77, 83, 89, 101, ... Depth beyond 3
adds almost nothing below `2^14` (depth at most 3 misses 1225 odd numbers there, the full closure 1214). Full miss
lists: `audit/b3_misses_closure_12.txt`, `audit/b3_misses_closure_14.txt`.

Reading. My expectation (C3, C4) was that the chains with a sink before every root change would cover the odd
numbers once three blocks are allowed. They do not: a pure join gives `N = A_1 A_2 + 2^i W_1 W_2` with `W >= A` on
both sides (cut a Hamiltonian path at the root: audit U5), so `N >= 3 A_1 A_2`; the inserted chain makes the second
term dominate and ties the low bits of `N` to the product `A_1 A_2` of two smaller counts. The values are numerous
(14 distinct chains per odd number; a Poisson model with that mean would miss a fraction `6 * 10^-7`) and badly
distributed. That the full grammar, with root changes not preceded by a sink (the term `- x_e y_e` of THM-4866),
is complete above 2647 in range and needs no bank was already in canon (HYP-9494 (D)); this test adds the size of
the gap between the sub-language and the grammar. It drops the bank and the direct joins together and cannot
separate their contributions.

Limitation of the script (audit m10): `depth()` builds its depth-1 table for `p <= X` only; the 16,328 one-apex
values of `[2^17, 2^18]` attained only for larger `p` are all depth-2 values, so no count changes.

### B4. Direct joins = the O/Z/J grammar (addendum D; KNOWN, re-computed)

Coordinates of a rooted tournament: `(A, W, d)`, `A = H`, `W = o + e`, `d = o - e`, where `o` (`e`) is the sum of
`H(X - P)` over the paths `P` from the root with an odd (even) number of vertices (`x_h`, `x_R`, `x_e` of THM-4866;
`(D, E, U) = (A, d, W)` in the coordinates of THM-4904).

**Proposition B3 (KNOWN: the case "one-apex block, both ports at the apex" of the seam transfer; the product of
the letter moves of the block. New: the closed entries. Brute force 1470 of 1470, `runs/d1_bruteforce.txt`; audit:
the state at every stage of 600 random and 3,234 exhaustive chains, and 797,161 O/Z/J words against the
products).** Let `(X, u)` have coordinates `(A, W, d)`, and let `Y` be a one-apex tournament with apex `v` and apex
word `w` of length `n`; `P` = value of its letters 1 with the first spine vertex as top bit, `Q` = value of its
letters 0 with the first spine vertex as unit bit. The join `X => Y` with the arc between `u` and `v` reversed,
rooted at `v`, has coordinates

    A' = g(P) A + (Q + 1) W - [Q even] d,      W' = (P + 1) A + 2^n W,      d' = [P even] A - [n = 0] d.

So the path counts of the family are the first coordinates of the vectors `M_(w_k) ... M_(w_0) (1, 0, 0)^T`, where

    M_w = [[ g(P), Q + 1, -[Q even] ], [ P + 1, 2^n, 0 ], [ [P even], 0, -[n = 0] ]].

Dictionary with the O/Z/J grammar (audit M1, `d4`): `J` = a new apex (the block with the empty word), `O` = a spine
letter 1, `Z` = a spine letter 0, the rooted edge = the block with word `1`.

In the form `(A, o, e)`, which is the one checked by brute force: `A' = A g + 2 (o E + e O)`,
`o' = A (1 + floor(P/2)) + o t_o + e t_e`, `e' = A ceil(P/2) + o t_e + e t_o`, with `O = 1 + floor(Q/2)`,
`E = ceil(Q/2)`, `t_e = t_o = 2^(n-1)` for `n >= 1` and `t_e = 1`, `t_o = 0` for `n = 0`. A chain of `i` vertices
before a block is the same as `i` leading letters 0 of its word. After the first block `0 <= d <= A <= W`
(KNOWN: Session200's inequalities `D - E >= d >= 1`, `D <= U`, quoted in THM-4904), and `A' >= W` at every step.

**Result (KNOWN; re-computed here).** This family is the O/Z/J root-change grammar. That it misses exactly the 29
odd numbers

    7, 13, 21, 41, 43, 55, 67, 77, 83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 583, 637,
    821, 1067, 1339, 1861, 2647

is Session197's `C1_COMPLETE.md` (to 5,001), THM-4884 (2b) (to 50,001) and HYP-9494 (D) (to `10^9`).
`scripts/direct_joins.py` re-computes the list below `2^18` in block coordinates (every value `A <= X` of the family
is generated; states with `o` or `e` above `X/2` are not stored, which loses no value: audit 6.2; 32,860,481 live
states at `2^18`); the audit re-computed it below `2^16` with an independent program and below `2^11` with two
more, one without any pruning. All 29 except 7 and 21 are path counts of other tournaments (KNOWN; explicit
witnesses with literal recounts: `audit/d3_witnesses_for_the_29.out`). Predictions: D1 and D3 CONFIRMED, D2 (only
7 and 21 are missed) FAILED; D2 and D3 were registered without knowledge of the canon record, so they were not
out-of-sample tests (audit M1 (iii)).

### B5. At most `L` spine letters per root change (addendum E; FINITE-EXACT)

`scripts/finite_blocks.py`, `scripts/finite_blocks_trend.py`: closure of the semigroup generated by the `M_w` with
`w` of length at most `L`; in O/Z/J terms, the words in which every run of spine letters has length at most `L`.

| `L` | blocks | odd numbers missed below `2^17` | largest miss | misses in `[2^k, 2^(k+1))`, `k = 8, ..., 16` |
|---|---|---|---|---|
| 3 | 15 | 1355 | 130247 | 29, 48, 58, 86, 119, 170, 257, 285, 277 |
| 4 | 31 | 127 | 130247 | 12, 17, 20, 25, 15, 11, 5, 3, 1 |
| 5 | 63 | 49 | 10391 | 9, 9, 8, 6, 1, 1, 0, 0, 0 |

Below `2^16`: `L = 0` gives 20 live states and misses almost everything; `L = 1` misses 28,170 odd numbers, `L = 2`
10,078. For `L = 5` the run was extended to `2^19` (`runs/finite_blocks_L5_19.txt`, 57,008,668 live states): the
same 49 misses, none in `(10391, 2^19]`. Predictions: E1 CONFIRMED; E2 (no miss in `[2^13, 2^16)` for `L = 3`) and
E3 (at most 60 misses for `L = 4`) FAILED. The audit reproduced the numbers below `2^16` and `2^17` as a by-product
(`audit/e1_finite_blocks_bonus.out`) and did not audit the statements.

Prior art (my check after the audit; a search, not a clearance): the values 10391 and 130247 occur nowhere else in
the repository; HYP-9494 (D) bounds the *number* of root changes (at most 7 above 2647), not the run length;
HYP-9640 (B) is a finitely generated statement of the same kind with other generators (two-port blocks of order at
most 4) and a lower threshold (493, to `2 * 10^5`). So B5 is a new measurement inside a known programme, not a new
reduction. HYP-9677 records it as a conjecture.

### B6. Where this leaves a proof of Chan-Pak

No proof, and no reduction that is not already in canon.

1. PROVED (THM-4885 (C)): a family with a fixed finite bank has bounded mean multiplicity. HEURISTIC (the Poisson
   reading of HYP-9480 (A), (C), OPEN): completeness up to `2^35` with a fixed bank is then not evidence that the
   same family is complete for all `N`.
2. KNOWN: families with growing multiplicity exist (joins of two one-apex tournaments: `(pi^2/12) ln N`; two-apex
   words: polynomial, opus-2026-11-11-S93), and positive density is still OPEN for all of them.
3. This session's measurement: many representations are far from sufficient (B3: 14 distinct chains per odd
   number, 5 per cent missed). What a proof has to control is the distribution of the counts in residue classes
   and short intervals, not their number.
4. The finitely generated statements (HYP-9640 (B); HYP-9677) are local-to-global questions for an orbit of a
   semigroup of integer matrices, of the kind studied for continued fractions (Zaremba's conjecture, where density
   one is known and the full statement is not). Nothing of that technology was used or checked here; the remark
   only names where such a statement lives. The thin-orbit lane of opus-2026-11-11-S93 took this up for
   HYP-9640 (B) (parity lemma, Zariski density).

## Reproduce

    cd scripts
    python euler_field.py                 # part A   -> runs/euler_field.txt
    python chain_semigroup.py c1          # C1       -> runs/c1_bruteforce.txt
    python chain_semigroup.py depth 18    # C2, C3, C5 -> runs/depth_18.txt
    python chain_semigroup.py closure 14  # C4       -> runs/closure_14.txt   (closure 12 -> runs/closure_12.txt)
    python direct_joins.py d1             # D1       -> runs/d1_bruteforce.txt
    python direct_joins.py closure 13     # D2, D3   -> runs/direct_closure_13.txt (also 11, 16; 18 needs about 2 GB)
    python finite_blocks.py 16 4          # E1-E3    -> runs/finite_blocks_16.txt
    python finite_blocks_trend.py 17      # L = 3, 4, 5 by dyadic block -> runs/finite_blocks_trend_17.txt
    python finite_blocks_trend.py 19 5    # L = 5 to 2^19 -> runs/finite_blocks_L5_19.txt

`euler_field.py` imports the engine of `../kps_s24_stein_modules_keller_20261111/scripts/`; `chain_semigroup.py`
needs numpy and numba. Audit programs: `audit/` (see `audit/AUDIT.md`, section 9).
