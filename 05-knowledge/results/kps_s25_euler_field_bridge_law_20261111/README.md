# Euler fields in twisted de Rham classes, and the one-apex chain semigroup (kind-pasteur-2026-11-11-S25)

Owner directive (2026-11-11): continue the session on Stein modules of Keller pairs (S24, THM-5138, HYP-9668) and
extend it into a session that aims to prove the Chan-Pak conjecture.

Canon: THM-5147 (part A), THM-5148 (part B); HYP-9668 updated ((E') REFUTED). Pre-registration:
`PREREGISTRATION.md` (commit 7479352f36, pushed before the runs of part B; C1, C2, C5 CONFIRMED, C3, C4 FAILED).
**JC(2) and the Chan-Pak conjecture remain OPEN. No proof of Chan-Pak was obtained.** Audit: `audit/AUDIT.md`
(pending at this commit: PRE-AUDIT DRAFT).

## 0. What comes out

- **Part A (Stein thread; PROVED).** An elementary identity for twisted de Rham classes: if a vector field `xi`
  satisfies `xi(W) = W` and has constant divergence `delta`, then `W h = (xi + delta) h` modulo the image. Three
  consequences: (1) **HYP-9668 (E') is false**: for Broughton's polynomial `W = x + x^2 y` the image has
  codimension one, does not contain 1, and is not a Mathieu-Zhao subspace (`W^m` lies in it for all `m >= 1`,
  `x W^m = m! x` never does); the same for every `W = x + x^(l+1) y`. So the Gaussian potential `xy` is not just the
  simplest rank-one potential with nonzero class of 1; in every other case tested the Mathieu property fails.
  (2) A radial identity for the Stein module of a Keller pair: `[(UV)^m h] = [(E + 1)(E + 2)...(E + m) h]`,
  `E = (U d_u + V d_v)/2`. (3) A one-line proof of the charge grading used in S24 (Proposition 6.3 there).
- **Part B (Chan-Pak; a formulation, a negative result, no proof).** The join law of rooted tournaments is KNOWN
  in the repo (THM-4866, Session263), and so are the one-apex families built on it (THM-4875, THM-4885). What is
  added: when a transitive chain is inserted at each join, the join law becomes a product of 2 x 2 integer matrices
  (PROVED from the known law; brute force 1470 of 1470), so the path counts of one-apex chains are the top-left
  entries of a matrix semigroup. Pre-registered test: this "pure" semigroup is **not** complete in range, by a wide
  margin (1214 odd numbers missed below `2^14` at full closure; 5.3 per cent missed in `[2^17, 2^18)` at depth 3),
  although its mean multiplicity there is 83. Two of my five predictions failed. High mean multiplicity does not
  give coverage; the joins without an inserted chain and the bank of small tournaments carry the coverage.

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
   subspace of `C[x, y]`.**

For `l = 1` this is Broughton's polynomial `x + x^2 y` (generic fibre `C^*`, no critical point, `dim B/M_W = 1`,
`[1] != 0`): the test case of HYP-9668 (E').

*Proof.* (1) `lambda(h) = Res_(x=0) int_0^infinity h e^(-W) dy dx`: for `h = x^k y^j` the inner integral along
the ray where `x^(l+1) y > 0` is `j! x^(k - (l+1)(j+1)) e^(-x)`, and the residue is the coefficient of
`x^(-1)`. It vanishes on `(d/dy - W_y) B` because the inner integral of a `y`-derivative is
`-h(x, 0) e^(-x)`, which has no residue, and on `(d/dx - W_x) B` because the residue of an `x`-derivative is zero.
(Direct check on monomials for `l = 1`: both families of generators are killed, term by term.)
(2) `xi = x d/dx - l y d/dy` satisfies `xi(x) = x`, `xi(x^(l+1) y) = x^(l+1) y`, so `xi(W) = W`, and
`div xi = 1 - l`. By Lemma A1, `[W^m] = (1 - l)(2 - l)...(m - l) [1]`, which is zero for `m >= l`; and
`xi(x^l) = l x^l` gives `[W^n x^l] = (l + 1 - l)(l + 2 - l)...(n) [x^l] = n! [x^l]`. (3) `lambda(x^l) = 1`. QED

Machine (`scripts/euler_field.py`, `runs/euler_field.txt`; the S24 engine, exact modulo a prime, membership
relative to the truncation): for `l = 1`: codimension 1, `[1] != 0`, `W h = (k - j) h` on 20 monomials,
`W^m` in `M` and `x W^m - m! x` in `M` for `m <= 4`, `x W^m` not in `M`; for `l = 2`: codimension 1, `W + 1`,
`W^2`, `W^4` in `M`, `x^2 W^n - n! x^2` in `M` for `n <= 4`, `[x^2] != 0`.

**Corollary A3 (a second obstruction to the Mathieu property; PROVED).** Let `W` have a polynomial vector field
`xi` with `xi(W) = W` and `div xi = delta` a non-positive integer. If some eigen-polynomial `h` of `xi`
(`xi h = e h`) has `[h] != 0` and `e + delta` not a non-positive integer, then `M_W` is not a Mathieu-Zhao subspace:
`P = W^(1 - delta)` has all powers in `M_W` and `h P^m` never lies in it.

So for the question (Q) of HYP-9668 ("for which potentials is the image a Mathieu-Zhao subspace") there are now two
necessary conditions that do not involve 1 being in the image only: the 1-property (S24, `x^p y^q`, `p != q`) and
the Euler obstruction. For the Gaussian potential `delta = 1`, and for `W = (f - a)(g - b)` of a Keller pair the
field `E = (U d_u + V d_v)/2` has `delta = 1` as well: the obstruction is absent exactly there.

What this does not settle: rank-one potentials without an Euler field, such as `x (x^2 y + x + 1)` (in the
classification of polynomials with general fibre `C^*`, the ones with a non-constant inner polynomial). The
statement "among rank-one potentials only `xy`, up to automorphisms, has a Mathieu-Zhao image" is consistent with
everything computed and is not proved.

### A3. The Stein module of a Keller pair

**Proposition A4 (radial identity; PROVED).** For an etale map `(u, v) : X -> C^2`, a centre `z`, `U = u - a`,
`V = v - b`, `E = (U d_u + V d_v)/2`, and every `h` in `O(X)`:

    [(UV)^m h]_z = [(E + 1)(E + 2) ... (E + m) h]_z,      and      [(U d_u - V d_v) h]_z = 0.

*Proof.* Lemma A1 in the frame `du ^ dv`, where `div E = 1`, `E(UV) = UV`, and the rotation field has divergence
0 and kills `UV`. QED

Checks: 7 of 7 on the Keller pair `(x + y^2, y)` at the centre (2, 3); 6 of 6 on the slice map of S24, section 6.2.

Uses. (i) The grading of the Stein kernel by charge at the special centre of S24's Proposition 6.3 (the plane minus
a hyperbola, `W = xy`): `eta = x d/dx - y d/dy` preserves `1 + xy`, kills `W` and has divergence 0, so every
element of nonzero charge is `eta` of something and lies in the kernel. (ii) `UV - c` is never a null element:
`[(UV - c)^m] = m! e_m(-c) [1]` with `e_m` the truncated exponential series, and consecutive truncations have no
common root. (iii) It is the radial half of the structure whose angular half is the rotation field; in the proof
of GMC(2) these are the two layers (THM-1645).

## Part B. The one-apex chain semigroup

### B1. What is known (credits)

- Odd-cycle formula `H = I(Omega, 2)` (Grinberg-Stanley; canon) and its increments: THM-4876
  (`H(T) - H(T - v) = 2 sum over odd cycles C through v of H(T - V(C))`, KNOWN).
- The join of two rooted tournaments (the ordinal sum `X => Y` with the arc between the two roots reversed):
  `H = x_h y_h + x_R y_R - x_e y_e` in rooted triples (THM-4866 (B)); Session263, Theorem RJ2 (return-join
  recursion). In the notation below `x_h = H(X)`, `x_R = W`, `x_e = o - e`.
- One-apex tournaments (one vertex whose removal leaves a transitive tournament): `H = g(v)`, `g(0) = 1`,
  `g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`, rooted triple `(g(v), e(v), 1 + v)` (THM-4875, THM-4885,
  THM-4890). One-apex x one-apex joins `(1+u)(1+v) + g(u) g(v) - e(u) e(v)`: mass `(pi^2/24) X ln X` (THM-4885
  (B)). Joins of a bank of 2,732 rooted triples of order at most 7 with one-apex tournaments reach every odd number
  up to `2^35` except 7 and 21 (THM-4885 (D), FINITE-EXACT); the mean multiplicity of a fixed bank is a constant
  (THM-4885, `M_Y`).
- CONCURRENT today: opus-2026-11-11-S93 (THM-5141: counting exponent 3/4 for one-apex counts; two-apex words of
  length at most 16 reach every odd number below `2^23` except 7 and 21), opus-2026-11-11-S94 (proof-oriented),
  Sessions 275 and 276 (digit banks, carries, corridors). Positive density of the spectrum is OPEN.

I re-derived the first three items from the odd-cycle formula before the sweep returned; they are not results of
this session.

### B2. The chain law

For a tournament `X` and a vertex `u` put `W^+_X(u) = sum over directed paths P starting at u of H(X - P)`,
`W^-_X(u)` the same with paths ending at `u`, and `tau(Z) = sum over all directed paths P of Z, the empty one
included, of H(Z - P)` (`H` of the empty tournament is 1).

**Proposition B1 (PROVED from the KNOWN increment formula THM-4876).** Let `(X, u)`, `(Y, v)` be rooted
tournaments and `T = X => C => Y` the ordinal sum with a transitive chain `C` of `i >= 1` vertices in the middle,
with the arc between `u` and `v` reversed (`v -> u`). Then

    H(T) = H(X) H(Y) + 2^i W^+_X(u) W^-_Y(v),          W^+_T(v) = H(X) W^+_Y(v) + 2^i W^+_X(u) tau(Y - v).

*Proof.* In `X => C => Y` the arc `u -> v` lies on no cycle, so by THM-4876
`H(T) = H(X) H(Y) + 2 sum over odd cycles Z through v -> u of H(T - Z)`. Such a cycle is `v -> u`, a path `P_X` of
`X` starting at `u`, an increasing subset `S` of `C`, a path `P_Y` of `Y` ending at `v`; its complement is
`(X - P_X) => (C - S) => (Y - P_Y)`, with `H = H(X - P_X) H(Y - P_Y)`; for each pair `(P_X, P_Y)` exactly
`2^(i-1)` subsets `S` give an odd length. A path of `T` starting at `v` either stays in `Y`, with complement
`X => C => (Y - P)`, or is `v -> u`, a path of `X` from `u`, a subset of `C`, a path of `Y - v` or nothing. QED

For a one-apex tournament rooted at its apex, with apex word of value `p` (letters: 1 if the apex beats the spine
vertex; first spine vertex = top bit), `L(p)` the bit length and `R(p)` the reversal of the complement of `p` on
`L(p)` bits: `H = g(p)`, `W^+ = p + 1`, `W^- = R(p) + 1`, `tau = 2^(L(p))`. So:

**Corollary B2 (PROVED; brute force 1470 of 1470 chains of depth 2 and 3, `runs/c1_bruteforce.txt`).** With

    M(p, i) = [[ g(p), 2^i (R(p) + 1) ], [ p + 1, 2^(i + L(p)) ]]        (p >= 0, i >= 1),

the chain `Y_0 => C_1 => Y_1 => ... => C_k => Y_k` of one-apex tournaments with apex words `p_0, ..., p_k`, chains
of `i_1, ..., i_k >= 1` vertices and the arcs between consecutive apexes reversed has
`H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1)`. (Leading letters 0 of a word are absorbed in the
exponent `i`.) The path-count spectrum contains the set `Sigma` of top-left entries of all products of the
`M(p, i)`.

Without inserted chains (`i = 0`) the law needs the parity split of the paths, i.e. the three-dimensional rooted
triples of THM-4866; that case is the known join law and is not treated here.

### B3. The pre-registered test (FINITE-EXACT)

`scripts/chain_semigroup.py`; outputs in `runs/`.

| prediction | outcome |
|---|---|
| C1: formula agrees with brute force | **CONFIRMED**: 1470 cases, 0 mismatches |
| C2: depth at most 2 misses more than 2 per cent of the odd numbers of `[2^17, 2^18)` | **CONFIRMED**: 11.06 per cent (7250 of 65536) |
| C3: depth at most 3 misses fewer than 0.5 per cent there, fewer than 600 below `2^18` | **FAILED**: 5.26 per cent (3446); 9190 below `2^18` |
| C4: full closure below `2^14` misses only 7 and 21 | **FAILED**: 1214 odd numbers missed (28,155 live states); 414 below `2^12` |
| C5: mean multiplicity at depth at most 3 at least 40 | **CONFIRMED**: 82.9 (2,659,887 products of depth 2 and 8,206,973 of depth 3 below `2^18`) |

First misses of the full closure: 7, 13, 21, 29, 35, 41, 43, 49, 55, 59, 67, 71, 77, 83, 89, 101, ... Depth beyond 3
adds almost nothing below `2^14` (depth at most 3 misses 1225 odd numbers there, the full closure 1214).

Reading. My expectation (C3, C4) was that the chain semigroup alone, without the bank of small tournaments, would
cover the odd numbers once three blocks are allowed. It does not, and the reason is visible in Proposition B1: a
pure join gives `N = A_1 A_2 + 2^i W_1 W_2` with `W >= A` on both sides, so `N >= 3 A_1 A_2`; the inserted chain
makes the second term dominate and ties the low bits of `N` to the product `A_1 A_2` of two smaller counts. The
values are numerous (mean multiplicity 83) and badly distributed. The joins without an inserted chain (the term
`- x_e y_e` of THM-4866) and the bank are what the completeness results in canon use, and this test shows they
cannot be dropped at these sizes.

### B4. Where this leaves a proof of Chan-Pak

No proof, and no reduction that is not already in canon. Three remarks for whoever continues (the first two are
KNOWN, the third is this session's measurement):

1. A family with a fixed finite bank has bounded mean multiplicity (THM-4885), so completeness up to `2^35` with a
   fixed bank is not evidence that the same family is complete for all `N`.
2. Families with growing multiplicity exist (joins of two one-apex tournaments: `(pi^2/12) ln N`; two-apex words:
   polynomial, opus-2026-11-11-S93), and positive density is still OPEN for all of them.
3. A large mean multiplicity is far from sufficient (part B3: 83 on average, 5 per cent missed). What a proof has to
   control is the distribution of the counts in residue classes and short intervals, not their number.

## Reproduce

    cd scripts
    python euler_field.py                 # part A   -> runs/euler_field.txt
    python chain_semigroup.py c1          # C1       -> runs/c1_bruteforce.txt
    python chain_semigroup.py depth 18    # C2, C3, C5 -> runs/depth_18.txt
    python chain_semigroup.py closure 14  # C4       -> runs/closure_14.txt   (closure 12 -> runs/closure_12.txt)

`euler_field.py` imports the engine of `../kps_s24_stein_modules_keller_20261111/scripts/`; `chain_semigroup.py`
needs numpy and numba.
