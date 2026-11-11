---
id: THM-5130
title: "Sign-flip towers of the Hamiltonian-path spectrum. (A) Gadget law: a pair (B, O) has a state s = (r, M, h, E) (offset, slope, H(B), paths ending in O), the closed tournament B => TT_t => v with v -> O has H = r + 2^t M, and for every gadget (X, P) the pair (B => X => v with v -> O, P) has state G(X, P) s for an explicit 4x4 integer matrix. (B) For X transitive and P two vertices at distance delta the new offset is rho_delta r + 2^(a-delta) (2^delta - 1 - 3^(delta-1)) M; it is -r for every base exactly when 2^delta - 3^(delta-1) = 1, that is delta = 1 or 2. (C) Towers: iterating gives multi-apex tournaments with H = (+1 or -1) + 2^t M_k whose slopes obey M_(k+1) = (2^delta+1) 2^(m+a_k) M_(k-1) + 2 (3^delta+2^delta) 2^(i-1) 3^m M_k - [i=1] E_k. (D) Explicit two-closer families: 1 + q 3^j 2^s (s >= 3) and 1 + (2 q 3^j - 1) 2^s (s >= 1) for q in {15, 25, 39, 65}; 2^s (3 2^m + 5 2^K 3^m) - 1 and 2^s (5 2^m + 13 2^K 3^m) - 1; rays 2^s - 1 from one pair. (E) For every odd c <= 4097 there is a seed (ray 1 + c 2^s) and an anti-seed (ray c 2^s - 1), each recounted literally; this gives seeds to the 33 values of HYP-9658 (B). (F) Unit sets lie in one strong component; W_ef >= H(B) across components; U_B(o) >= H(B); twin lemma"
status: >
  (A) PROVED (path surgery). KNOWN parts: the closing law is the cut form of the tail law (Session216 star closure,
  THM-4960, THM-5122 (1), THM-5126 (A)); row h of the matrix, H(B => X => v) = r H(X) + M Pi_2(X), is Session216
  star/PROOF.md Theorem 3 (a special case of the substitution law THM-3121 (4); linearity in the bin polynomial of a
  module is Session197 UNIVERSAL_CONTEXTS section 1); the four coordinates are those of the Session216 interface
  (h, alpha = M, tau = h - r, e = E). New here, as far as a canon search shows: the whole state is transported
  linearly through a new closer, with the matrix (rows r, M, E).
  (B) PROVED. For a transitive base the offsets 1 - 2^(delta-1) and the two lines 3 2^j - 1, 5 2^j - 1 are the
  shared-endpoint case of THM-4748 (1), (2). New: the coefficient of M for an arbitrary base and the criterion
  2^delta - 3^(delta-1) = 1 (stated in the pre-registration of this session, pushed 2026-11-10 16:20). CONCURRENT:
  Session260 (coverage/PROOF.md, 2026-11-11) proves the same bracket, in the form 3^c - 2^(c+1) + 1 with
  c = delta - 1, for a selected pair inside a transitive module of one five-vertex carrier (unit iff c = 0 or 1).
  (C) PROVED from (A), (B).
  (D) PROVED from (C); every family is also recounted on a grid by a transfer-matrix count that shares no code with
  the subset count. The one-closer lines are KNOWN (THM-4748 (2)). CONCURRENT / KNOWN for D3: the unit slopes
  15, 25, 39, 65 times 3^b are in Session260 coverage/PROOF.md, and the lines 1 + (b 3^a - 1) 2^s for b = 30, 50, 78
  in endpoint_coverage_compiler_20261111; here they are corollaries of the tower recursion (b = 130 is added). That 2^k - 1 is a path count for every k >= 4 is
  KNOWN (THM-4748 (4), the Fibonacci fan, a one-apex tournament); new is only that a single pair has offset -1 and a
  power of two as slope, so that one closing family gives 2^s - 1 for all s >= 11.
  (E) FINITE-EXACT + the KNOWN closing law: 2,049 + 2,049 certificates (certs_plus.txt, certs_minus.txt), each state
  recounted by subset dynamic programming on the explicit tournament. The seed half is KNOWN: 2,016 seeds are in
  THM-5126, and unit-offset roots of order eleven for the 33 exceptional values of HYP-9658 (B) were published first
  by the session of 05-knowledge/results/odd_c_seed_symmetry_20261110/ (2026-11-10), which settled HYP-9658 (B); the
  33 tower seeds here are a second, independent set of witnesses (rooted orders 11 to 15). New in (E), as far as a
  canon search shows: the anti-seeds (offset -1) for every odd c <= 4097.
  (F) PROVED, elementary. KNOWN (pushed first by other sessions while this one was running): the strong-core
  structure F2 (c) (Session250, 2026-11-10 16:53; Session258 CORE_CLASSIFICATION), M >= |O| H(B) (Session258), the
  twin identity F3 (Session260). F1 in the exact form, F2 (a), (b) and F4 were not found elsewhere.
  (H) H1, H2 PROVED (H2 uses the domination inequality of Session261, EXTENSION4.md); H3 FINITE-EXACT; H4 PROVED
  (lower bound KNOWN, Session197).
  (G) (classification of neutral gadgets) FINITE-EXACT, orders <= 7.
  Pre-registered runs (PREREGISTRATION.md, pushed before each run): R2, R4, R5, R6c CONFIRMED; R1, R3, R6a, R6b
  FAILED (the towers reach more than predicted: a power of two as slope at depth 3, and c = 5 at depth 6).
  Independent audit: PENDING at this commit.
  The Chan-Pak conjecture and HYP-9658 (A) (a seed for every odd c) remain OPEN. Nothing here gives a positive
  density of path counts.
session: kind-pasteur-2026-11-10-S22
source: 05-knowledge/results/kps_s22_sign_flip_towers_20261110/README.md
scripts:
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/scripts/verify_s22.py (+ verify_s22.out)
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/scripts/towers.py, tower_build.py, mersenne_towers.py, neutral.py
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/scripts/pairstate.c, certs_verify.py, make_certs.py, portrait.c, tower_reach.c, ray_cover.py
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/certs_plus.txt, certs_minus.txt, runs/
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/scripts/sources_check.py (+ .out), singles.c, src10.c
  - 05-knowledge/results/kps_s22_sign_flip_towers_20261110/SOURCES_AND_DOMINATION.md
related:
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md (setting, lines, seeds)
  - 05-knowledge/hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md (a seed for every odd c)
  - 05-knowledge/hypotheses/HYP-9662-two-sided-seeds-and-the-reach-of-sign-flip-towers.md
  - 01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md
  - 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md (substitution law, eq. (4))
  - 05-knowledge/results/session216_overnight_constructions_20261107/ (star/PROOF.md Theorem 3; three-core/PROOF.md (3)-(4))
  - session197/path-polynomials/UNIVERSAL_CONTEXTS.md, UNIT_OFFSET_RIGIDITY.md; session196/spectral-moves/PROOF.md (9)
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md (the exclusions)
  - 05-knowledge/results/odd_c_seed_symmetry_20261110/ (order-eleven roots for the 33 lines; settled HYP-9658 (B) first)
  - 05-knowledge/results/session250_sixteen_20261110/, session258_endpoint_seed_20261111/ (strong cores of unit sets)
  - 05-knowledge/results/session260_quasiperfect_chanpak_20261111/ (coverage/PROOF.md; bridges/BIPARTITE_DEFECT.md)
  - 05-knowledge/results/session261_creative_chanpak_20261111/ (bridges/EXTENSION4.md; constructors/NINE_OFFSET_EXTENSION.md)
  - 05-knowledge/results/multiplicand_saturation_chanpak_20261111/ (chanpak/OFFSET_TRANSPORT.md; exchange/PROOF.md)
  - 05-knowledge/results/endpoint_coverage_compiler_20261111/ (lines b 3^a - 1)
---

# THM-5130 - sign-flip towers: gadget law, the equation 2^delta - 3^(delta-1) = 1, seeds of both signs

## Setting

`H(T)` is the number of Hamiltonian paths of a tournament `T`; `S` is the set of its values. `A => B` means that
every vertex of `A` beats every vertex of `B`; `TT_a` is the transitive tournament on `1 -> 2 -> ... -> a`.

A *pair* `(B, O)` is a nonempty tournament `B` with a vertex set `O` (possibly empty). As in THM-5126:

- `Q_B(e, f)` = number of ordered pairs `(P1, P2)` of vertex-disjoint nonempty paths covering `B`, `P1` ending at `e`,
  `P2` starting at `f`; `W_B(O) = sum over ordered pairs e != f in O of Q_B(e, f)`;
- `U_B(o)` = number of ordered pairs `(P1, P2)` covering `B` with `P2` starting at `o`, `P1` possibly empty;
  `U_B(O) = sum_{o in O} U_B(o)`; `V_B(o)` = the mirror count (`P1` ending at `o`, `P2` possibly empty);
- `Start_B(O)`, `End_B(O)` = number of Hamiltonian paths of `B` starting, ending in `O`;
- `Pi_2(B)` = number of ordered pairs of vertex-disjoint possibly empty paths covering `B`.

The **state** of the pair is `s(B, O) = (r, M, h, E)` with

    h = H(B),   r = h - W_B(O)  (offset),   M = U_B(O)  (slope),   E = End_B(O).

For a tournament `X`, the **closing** `C(B, O; X)` is the tournament on `V(B) + V(X) + {v}` with `B => X => v`,
`v -> O`, `B - O -> v`; the new vertex `v` is the *closer*. `T_t(B, O) = C(B, O; TT_t)`.

**Closing law (KNOWN).** `H(T_t(B, O)) = r + 2^t M` for `t >= 1`.

A pair with `r = 1` is a *seed* (it puts the ray `1 + M 2^t` into `S`), a pair with `r = -1` an *anti-seed* (ray
`M 2^t - 1`). Every odd `N > 1` is `1 + c 2^s` and `c' 2^s' - 1` in exactly one way (`c`, `c'` odd), and
`min(s, s') = 1`.

A **gadget** `(X, P)` is a nonempty tournament `X` with a vertex set `P`. Its constants:

- `H_X = H(X)`, `Pi = Pi_2(X)`, `w = W_X(P)`;
- `w3 = sum over ordered pairs e != f in P` of the number of ordered triples `(D1, D2, D3)` of disjoint nonempty paths
  covering `X`, `D1` ending at `e`, `D2` starting at `f`;
- `sigma = Start_X(P)`, `u = U_X(P)`, `eps = End_X(P)`, `nu = V_X(P)`;
- `u3 = sum_{o in P}` of the number of ordered triples `(D, D', D'')` of disjoint nonempty paths covering `X`, `D`
  starting at `o`;
- `r_X = H_X - w` and `kappa_X = Pi - 2 w - 2 w3`.

## (A) The gadget law

**Theorem A.** Let `(B, O)` be a pair with state `(r, M, h, E)`, `(X, P)` a gadget, and `B+ = C(B, O; X)`. Then the
state `(r+, M+, h+, E+)` of the pair `(B+, P)` is

    r+ = r_X r + kappa_X M
    M+ = u (r + h) + (4 u - 2 sigma + 2 u3) M - sigma E
    h+ = H_X r + Pi M
    E+ = eps (r - h) + nu M.

*Proof.* Contract `X` to one vertex `z`: `L' = B + z + v` with `B => z -> v`, `v -> O`, `B - O -> v`. In a family of
disjoint paths covering `B+`, the vertices of `X` occur in maximal blocks; a block is initial in its path or entered
from `B`, and terminal in its path or left towards `v` (because `B => X => v` and `v` beats no vertex of `X`).
Contracting blocks gives walks of `L'` that cover `V(B) + {v}` once and visit `z` some number `k` of times, and the
blocks form an ordered `k`-tuple of disjoint nonempty paths covering `X`; the correspondence is bijective. In `L'`:

- (o1) a visit of `z` is first in its walk or preceded by a vertex of `B`; it is last in its walk or followed by `v`;
- (o2) `v` occurs once, so at most one visit is followed by `v`; `v` is last in its walk or followed by a vertex of
  `O`, and first in its walk or preceded by `z` or by a vertex of `B - O`.

Write `(Y', Y)` for a pair of disjoint paths covering `B`, and recall `M = U_B(O)` (`Y'` possibly empty, `Y` nonempty
starting in `O`), `W = W_B(O) = h - r`.

*One walk.* One visit: `Y' z v Y` (`Y` nonempty: `M`; `Y` empty: `h`) or `Y' v Y z` (`Y` nonempty starting in `O`,
`Y'` empty or ending outside `O`: `M - W`). Total `r + 2M`. Two visits: `Y' z v Y z`, `Y` nonempty: `M`. Three visits
are impossible by (o2). Walks ending at `z`: one visit `M - W = M - h + r`, two visits `M`.

*Two walks, the first ending at a visit and the second starting at one.* Gluing the two end visits into one is a
bijection onto single walks with a marked visit. So there are `r + 2M` such pairs with two visits and `2M` with three.

*Two walks `(W1, W2)`, `W1` possibly empty, `W2` starting at a visit.* Case I: `W2` begins `z v`. Then every other
visit is last in its walk: `W2 = z v Y` or `z v Y z` (`Y` nonempty in the second form), `W1` is empty, `Y'`, `Y' z`
(`Y'` nonempty) or `(z)`. Case II: `W2 = (z)`; then `W1 = alpha v beta` with `alpha` empty, `Y'` (ending outside
`O`), `(z)` or `Y' z`, and `beta` empty, `Y` or `Y z`. Counting by the number of visits, with `S = Start_B(O)`:

    one visit:    (M + h) + (M - W + h - E)                              = 2M + h + r - E
    two visits:   (M + (M - S) + h + S) + (S + (M - S) + h + (M - W))    = 4M + h + r
    three visits: ((M - S) + S) + (S + (M - S))                          = 2M.

*The four coordinates.* `h+`: `(r + 2M) H_X + M (Pi - 2 H_X)`, because two blocks form an ordered pair of nonempty
paths. `W_{B+}(P)`: the block ending `W1` ends at `e`, the block starting `W2` starts at `f`, so
`W_{B+}(P) = (r + 2M) w + 2M w3`, and `r+ = h+ - W_{B+}(P)`. `M+`: the block at the start of `W2` starts at `o`; with
one, two, three visits the blocks are counted by `sigma`, `u - sigma`, `u3`, so
`M+ = (2M + h + r - E) sigma + (4M + h + r)(u - sigma) + 2M u3`. `E+`: walks ending at `z`, the last block ending at
`o`: `(M - h + r) eps + M (nu - eps)`. ∎

**Remarks.** (1) The same matrix applied to `(1, 0, 1, 0)` gives the state of `(X => pt, P)` (direct check), so the
empty base may be given this state. (2) Adding a sink to `B` changes the state to `(r, 2M, h, 0)`; adding a source
changes nothing (direct; verify_s22.py V4). (3) Row `h` with `X = TT_t` is the closing law.

## (B) Transitive pairs and the equation 2^delta - 3^(delta-1) = 1

**Theorem B.** Let `X = TT_a`, `P = {i, i + delta}`, `m = a - i - delta >= 0`. Then

    r_X = -1 (delta = 1),   r_X = 1 - 2^(delta-1) (delta >= 2),   kappa_X = 2^(a-delta) (2^delta - 1 - 3^(delta-1)),
    u = (2^delta + 1) 2^m,   4u - 2 sigma + 2 u3 = 2 (3^delta + 2^delta) 2^(i-1) 3^m,   sigma = [i = 1],
    eps = [m = 0],   nu = (2^delta + 1) 2^(i-1),   H_X = 1,   Pi = 2^a.

*Proof.* In `TT_a` a set of vertices carries exactly one path, so covers are set partitions. `Q(i, i+delta)`: the
first block lies in `{1..i}`, the second in `{i+delta..a}`: one cover if `delta = 1`, none otherwise.
`Q(i+delta, i)`: vertices below `i` go to the first block, above `i + delta` to the second, the `delta - 1` between
to either: `2^(delta-1)`. So `w = 2` for `delta = 1` and `2^(delta-1)` for `delta >= 2`. With a third nonempty block:
`2^(i-1+m) - [delta = 1]` and `2^(i-1) 3^(delta-1) 2^m - 2^(delta-1)`, so
`w3 = 2^(a-delta-1) (1 + 3^(delta-1)) - [delta = 1] - 2^(delta-1)` and
`kappa_X = 2^a - 2w - 2 w3 = 2^a - 2^(a-delta)(1 + 3^(delta-1))`. For a vertex `j`: `U(j) = 2^(a-j)`,
`u3(j) = 2^(j-1) 3^(a-j) - 2^(a-j+1) + [j = 1]`, `V(j) = 2^(j-1)`; summing over `j = i, i + delta` gives the rest. ∎

**Corollary B (sign flip).** For a transitive pair at distance `delta` put at the tail of an arbitrary pair `(B, O)`,

    r+ = r_X r + 2^(a-delta) (2^delta - 1 - 3^(delta-1)) M.

Since `M > 0` whenever `O` is nonempty, `r+ = -r` holds for every base if and only if `2^delta - 3^(delta-1) = 1`,
that is `delta = 1` or `delta = 2` (for `delta >= 3`, `3^(delta-1) > 2^delta - 1`). For these two distances the level
map is, with `alpha = 2^delta + 1`, `gamma = 3^delta + 2^delta` (so `(alpha, gamma) = (3, 5)` or `(5, 13)`):

    r+ = -r,   h+ = r + 2^a M,   M+ = alpha 2^m (h + r) + 2 gamma 2^(i-1) 3^m M - [i = 1] E,
    E+ = alpha 2^(i-1) M + [m = 0] (r - h).                                                      (1)

## (C) Towers

A **tower** over a pair `(B_0, O_0)` with levels `(delta_j, i_j, m_j)`, `j = 1..k`, `delta_j` in `{1, 2}`, is defined
by `B_j = C(B_(j-1), O_(j-1); TT_(a_j))`, `a_j = i_j + delta_j + m_j`, and `O_j` = the vertices `i_j`, `i_j + delta_j`
of the tail `TT_(a_j)`. When `B_0` is transitive, `T_t(B_k, O_k)` is a tournament with a transitive chain
`Z_0 => Z_1 => ... => Z_(k+1)` and `k + 1` closers `v_1, ..., v_(k+1)`: the closer `v_j` beats its out-set (inside
`Z_(j-1)`), all chain vertices after `Z_j`, and all later closers, and loses to everything else. For `j >= 2` the
out-set is two chain vertices at distance 1 or 2.

**Theorem C.** (i) The state of `(B_k, O_k)` is obtained from that of `(B_0, O_0)` by `k` applications of (1); its
offset is `(-1)^k r_0`. (ii) `h_(j+1) + r_(j+1) = 2^(a_(j+1)) M_j`, hence for `j >= 1`

    M_(j+1) = alpha 2^(m + a_j) M_(j-1) + 2 gamma 2^(i-1) 3^m M_j - [i = 1] E_j        (2)

with `(delta, i, m)` the parameters of level `j + 1`. (iii) `H(T_t(B_k, O_k)) = (-1)^k r_0 + 2^t M_k`, and the base
of level `j + 1` is the closed tournament of level `j`: `H(B_(j+1)) = r_j + 2^(a_(j+1)) M_j`.

*Proof.* (i) is Corollary B. (ii): `h_(j+1) = r_j + 2^a M_j` and `r_(j+1) = -r_j`; insert into (1). (iii) is the
closing law and row `h`. ∎

So the slope of one level (a count of two-path covers of that level's base) is, up to the offset and a power of two,
the Hamiltonian path count of the next base; and the slopes satisfy a three-term recurrence whose coefficients are
`3 2^x`, `5 2^x` and `10 2^x 3^y`, `26 2^x 3^y`.

## (D) Explicit families

Transitive starts: `s(TT_n, {p}) = (1, 2^(n-p), 1, [p = n])` and, for `delta` in `{1, 2}`,
`s(TT_a, {i, i+delta}) = (-1, (2^delta + 1) 2^m, 1, [m = 0])` (`m = a - i - delta`).

**D1 (one closer; KNOWN, THM-4748 (2)).** A transitive chain and one vertex that beats exactly two chain vertices at
distance `delta` in `{1, 2}`, with `n' >= 1` chain vertices after them: `H = (2^delta + 1) 2^n' - 1`.

**D2 (two closers over a fan; anti-seeds).** For `n >= 1`, `1 <= p <= n`, `delta` in `{1, 2}`, `i >= 1`, `m >= 0`,
`t >= 1`, the tower over `(TT_n, {p})` with one level `(delta, i, m)` has

    H = 2^t [ (2^delta + 1) 2^(m+1) + 2 (3^delta + 2^delta) 2^(i-1+n-p) 3^m - [i = 1][p = n] ] - 1.

In particular `2^s (3 2^m + 5 2^K 3^m) - 1` and `2^s (5 2^m + 13 2^K 3^m) - 1` are path counts for all `s >= 2`,
`K >= 1`, `m >= 0`.

**D3 (two closers over a transitive pair; seeds).** For `delta_1, delta_2` in `{1, 2}`, the tower over
`(TT_a, {i_1, i_1 + delta_1})` with one level `(delta_2, i_2, m_2)` has, with
`q = (2^delta_1 + 1)(3^delta_2 + 2^delta_2)` in `{15, 39, 25, 65}`,

    H = 1 + 2^t [ 2 q 2^(m_1 + i_2 - 1) 3^(m_2) - [i_2 = 1][m_1 = 0] ].

In particular, for every `j >= 0`: `1 + q 3^j 2^s` is a path count for all `s >= 3`, and `1 + (2 q 3^j - 1) 2^s` for
all `s >= 1`.

**D4 (rays 2^s - 1).** The pair `(B, O)` with `B` of order 8, out-lists `0:[1]; 1:[2,4]; 2:[0,3,4]; 3:[0,1,7];
4:[0,3,5,6]; 5:[0,1,2,3]; 6:[0,1,2,3,5]; 7:[0,1,2,4,5,6]` and `O = {1, 3, 4}` has state `(-1, 2^10, 217, 64)`. So the
closed tournaments, on `9 + t` vertices, have exactly `2^(10+t) - 1` Hamiltonian paths (`t >= 1`). The towers over
`(TT_n, {1})`, `1 <= n <= 7`, with levels `(1, 8 - n, 0), (1, 4, 0), (2, 2, 1)` have offset `-1` and slope `2^23`:
tournaments on `N >= 24` vertices with four closers and `2^N - 1` paths. (That every `2^k - 1`, `k >= 4`, is a path
count is KNOWN: THM-4748 (4).)

*Proofs.* D1-D3: insert the start states into (1). D4: FINITE-EXACT state (verify_s22.py V8a, pairstate.c) and
the closing law; the towers by (1) and by the transfer-matrix count (mersenne_towers.py). ∎

## (E) Seeds and anti-seeds for every odd c <= 4097

**Theorem E (FINITE-EXACT).** For every odd `c <= 4097` there are a pair with `r = +1` and slope `c 2^v` and a pair
with `r = -1` and slope `c 2^v'`. Hence there are `s_+(c)`, `s_-(c)` with `1 + c 2^s` in `S` for all `s >= s_+(c)`
and `c 2^s - 1` in `S` for all `s >= s_-(c)`.

Certificates: `certs_plus.txt`, `certs_minus.txt` (2,049 lines each; `v <= 7`, `v' <= 10`; pairs of order at most
14). 1,969 seeds and 1,959 anti-seeds are pairs of order at most 9; the others are towers of one or two levels over
such pairs. Every certificate is rebuilt as an explicit tournament and its state recounted (certs_verify.py with
pairstate.c: 0 mismatches). Among them are tower seeds for the 33 odd `c <= 4097` of HYP-9658 (B), with rooted orders
11 to 15; these 33 lines already have roots of order eleven (odd_c_seed_symmetry_20261110, first), so the seed half of
this theorem is a second proof by another construction. The anti-seed half is the new part.

## (F) Structure of unit sets

**F1.** `U_B(o) = H(B) + sum over out-neighbours x of o of Q_B(x, o)`. In particular `U_B(o) >= H(B)` and
`M >= |O| H(B)`.

*Proof.* A Hamiltonian path starts at `o` or is cut before `o` into `(P1, P2)` with `P1` ending at an in-neighbour of
`o`: `H(B) = Start_B(o) + sum_{x -> o} Q_B(x, o)`, while `U_B(o) = Start_B(o) + sum_{x != o} Q_B(x, o)`. ∎

**F2.** Let `C_1 => ... => C_s` be the strong components of `B`. (a) If `e` in `C_a`, `f` in `C_b`, `a < b`, then
`Q_B(f, e) = prod_{i<a} H(C_i) * U_{C_a}(e) * prod_{a<i<b} Pi_2(C_i) * V_{C_b}(f) * prod_{i>b} H(C_i) >= H(B)`.
So `r(B, O) <= 0` whenever `O` meets two components. (b) If `O` lies in `C_a`:
`r(B, O) = r(C_a, O) prod_{i != a} H(C_i)` and `U_B(O) = U_{C_a}(O) prod_{i<a} H(C_i) prod_{i>a} Pi_2(C_i)`.
(c) Hence a unit set (`r = 1`) with `|O| >= 2` lies in one strong component `Q`, every other component is a single
vertex, `(Q, O)` is a unit pair, and `U_B(O) = 2^(number of vertices after Q) U_Q(O)`.

*Proof.* (a) `P1` ends in `C_b` and `P2` starts in `C_a`, so components before `a` lie in `P1`, after `b` in `P2`,
and each `C_i`, `a <= i <= b`, is split into an ordered pair of possibly empty paths, with the end conditions in
`C_a` and `C_b`. The bound uses F1, its mirror `V >= H`, and `Pi_2 >= 2H`. (b) is the same bookkeeping. (c) By (a)
and (b), `1 = r(Q, O) prod H(C_i)` forces all other components to have `H = 1`. ∎

**F3 (twin lemma).** For a tournament `K` and a vertex `y`, let `K_2 = K[TT_2 at y]` with twins `e -> f`. Then
`W_ef = 2 H(K)` and `r(K_2, {e, f}) = H(K_2) - 2 H(K)`. So the twins are a unit pair iff `H(K_2) = 2 H(K) + 1`.

*Proof.* `Q(e, f)` counts the Hamiltonian paths of `K_2` in which `e` is followed by `f`, that is `H(K)`. The covers
counted by `Q(f, e)` do not use the arc `ef`, and exchanging the names of the twins maps them to those of
`Q(e, f)`. ∎

**F4.** `W_ef` does not depend on the orientation of the arc between `e` and `f`. So if `B'` is `B` with that arc
reversed, `r(B, {e,f}) - r(B', {e,f}) = H(B) - H(B')`, and for `e -> f` the pair is unit iff the number of
Hamiltonian paths of `B` avoiding the arc `ef` is one more than `Q_B(f, e)`.

## (G) Neutral gadgets (FINITE-EXACT)

Among all gadgets `(X, P)` with `X` of order at most 7 (one representative per isomorphism class, every nonempty
`P`), `kappa_X(P) = 0` occurs only for `|P| = 2`, and then `r_X(P) = -H(X)`; those with `H(X) = 1` are exactly the
transitive pairs at distance 1 or 2 (`2a - 3` of them at order `a`). So, up to order 7, these are the only gadgets
that transport a unit offset to a unit offset for every base (neutral.py, neutral.out; orders at most 5 over all
labelled tournaments in verify_s22.py V9).

## (H) Sources, the domination skeleton, and an abundancy index

Details, credit and discussion: SOURCES_AND_DOMINATION.md of the packet; checks: sources_check.py (S1-S8).

**H1 (one-vertex gadget).** For a pair with state `(r, M, h, E)` and the gadget `(TT_a, {p})`:
`r+ = h+ = r + 2^a M`, `M+ = 2^(a-p)(r + h) + 2^p 3^(a-p) M - [p = 1] E`, `E+ = [p = a](r - h) + 2^(p-1) M`.
(Theorem A with `w = w3 = 0`, `u = 2^(a-p)`, `u3 = 2^(p-1) 3^(a-p) - 2^(a-p+1) + [p = 1]`.) For a seed of slope
`2^k` this is a pair with offset `1 + 2^(a+k)`; the pairs (offset, slope) = (5, 9), (17, 18), (33, 36) used by
Session261 for `1 + 2^s + 9 * 2^(r+s)` arise this way.

**H2 (domination ray).** For a pair `(B, O)`, `O` nonempty, let `G_t` be the bipartite graph of the closed tournament
`T_t(B, O)` at its closer (left side `O`, right side the other vertices of `B` and the tail, `o ~ x` iff `o -> x`).
Its number of dominating sets is `d(G_t) = rho_dom + 2^t delta_dom`, `t >= 1`, with
`delta_dom = sum over nonempty S in O of 2^(number of vertices outside O beaten by some vertex of S)`.
With the inequality `H(T) - H(T - z) >= d(G_z) - 1` of Session261: `(r - h - rho_dom + 1) + 2^t (M - delta_dom) >= 0`
for all `t >= 1`; in particular `M >= delta_dom`, and `U_B(o) >= 2^(outdeg o)`.

*Proof.* A tail vertex is adjacent to all of `O`. Dominating sets containing the whole tail:
`D_1 = sum over all S_O in O of 2^|N(S_O)|` (the vertices outside `O` not beaten by `S_O` are forced, the others
free). Those with a proper nonempty tail part: `(2^t - 2)(D_1 - 1)`, since `S_O` must be nonempty. Those with empty
tail part do not depend on `t`. So `delta_dom = D_1 - 1`. ∎

**H3 (FINITE-EXACT).** (i) `M = delta_dom` exactly when `B` is transitive (all pairs of order at most 6).
(ii) Pairs with offset `1 + 2^s` and slope `9 * 2^v`, `v <= s`, exist for every `2 <= s <= 11` among selected sets
over tournaments of order at most 9; so `1 + 2^s + 9 * 2^(r+s)` is a path count for all `r >= 1`, `s <= 11`
(`s <= 6`: Session261). (iii) The unit root of order 9 with code 45174751748 and root 2 has a port (vertex 8,
in-degree 1, not selected) with slope family `16 + 172 * 2^a`; with pairs of offsets 3, 5, 9, 17 and slopes 172, 172,
86, 86 and the count 175 of an order-7 tournament this gives: `1 + 2^s + 43 * 2^(r+s)` is a path count for all
`r, s >= 1`. (A module at an unselected port of in-degree 1 is visited at most twice by the paths counted in `H`,
`W` and `U`, so these are affine in `2^a` and a unit offset at `a = 1, 2` persists.)

**H4 (abundancy).** `H(X) >= H(X[S]) H(X - S)` for every vertex set `S` (two disjoint paths merge, keeping both
orders, into one path, injectively). Hence `a(X) = Pi_2(X)/H(X) = sum over S of H(X[S]) H(X - S)/H(X) <= 2^n`, with
equality iff `X` is transitive; `a(X) >= n + 1` (Session197); `a(X => Y) = a(X) a(Y)`; and the closing law reads
`H(C(B, O; X)) = H(X) (r + M a(X))`.

## Controls

- Gadget law: 300 random bases and gadgets with arbitrary offsets, no failure; closed forms of (B) for
  `delta <= 5`; sink and source rules (verify_s22.py V1-V4).
- Three counts agree on 77 random towers (subset count, transfer-matrix count, recursion); families D1-D3 on
  702 tournaments by the transfer-matrix count (V6, V7).
- F1-F4 on all 1,099 labelled tournaments of order at most 5 (2,166 unit sets) and random twins (V5).
- No pair `(B, O)` of order at most 7 (62,422 pairs, all offsets) has a ray `r + 2^t M` through 7 or 21, as
  THM-1370 demands; 13, 19, 27 and 1039 have 15, 17, 110 and 156 such rays (rays_through.py, runs/rays_through.out).
- The subset programme pairstate.c was checked on transitive tournaments of orders 12 to 27: it is correct to
  order 25 and wrong at 27 (allocation size); certs_verify.py refuses orders above 24. The certificates have order
  at most 14.

## What is not proved

- Whether every odd `c` has a seed (HYP-9658 (A)) and an anti-seed (HYP-9662 (A)).
- Whether towers reach every odd `c`. At a fixed depth the reachable `c` form a thin set (sums of a bounded number
  of terms `2^x 3^y` with coefficients from 3, 5, 13).
- Rays of offset `+1` or `-1` do not exhaust the spectrum. Example (PROVED): `15 = 1 + 7 * 2 = 2^4 - 1` is a path
  count on no such ray. A seed has slope `2^j` (`|O| <= 1`), `6 * 2^j` (core `C_3`) or at least 10 (by F1 and F2 a
  core with `H >= 5` gives `M >= 10`), never 7. An anti-seed has `M >= 2 H(B)` by F1, so `M <= 8` forces
  `H(B) = 1` or `3`; by F2 its slope is then `3 * 2^j`, `5 * 2^j` (transitive `B`) or `9 * 2^j`, `15 * 2^j`
  (`B = TT => C_3 => TT`, one vertex of `O` in the 3-cycle and the other next to it), never 8, 4, 2 or 1.
  At least 3,722 of the 4,095 odd numbers in `[3, 8191]` lie on a certified ray (ray_cover.py).
- The Chan-Pak conjecture.
