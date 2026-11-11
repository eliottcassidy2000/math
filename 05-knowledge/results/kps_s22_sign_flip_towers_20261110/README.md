# Sign-flip towers: constructing seeds for the 2-adic lines of the path spectrum

kind-pasteur-2026-11-10-S22. Owner directive: *construct a seed for every odd `c`; see how prior work on the numbers
2, 3, 5, 7, 13, 19, 27 and 1039 relates to the symmetry structure of tournaments, primes and Collatz-like
constructions; think of path counts in k-vertex extensions and how Chan-Pak could be proved by turning them into
Hamiltonian path counts of tournaments themselves.*

Canon: [THM-5130](../../../01-canon/theorems/THM-5130-sign-flip-towers-gadget-transfer-law-and-seeds-of-both-signs.md),
[HYP-9662](../../hypotheses/HYP-9662-two-sided-seeds-and-the-reach-of-sign-flip-towers.md). Continues
[THM-5126](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md)
and [HYP-9658](../../hypotheses/HYP-9658-unit-offset-seeds-of-the-two-adic-lines.md).

## 0. Result in one page

| statement | status |
|---|---|
| A pair `(B, O)` has a state `(r, M, h, E)`; putting any gadget `(X, P)` at its tail and adding a closer acts on the state by an explicit 4x4 integer matrix | PROVED (row `h` KNOWN: Session216; concurrent: ONE_EXIT.md and the marked tensor of Session261, which give row `h` and the existence of a transport, not the matrix) |
| For a transitive gadget and two vertices at distance `delta`, the offset is reversed for every base iff `2^delta - 3^(delta-1) = 1`, i.e. `delta` = 1 or 2 | PROVED |
| Towers: slopes obey `M_(k+1) = alpha 2^(m+a_k) M_(k-1) + 2 gamma 2^(i-1) 3^m M_k - [i=1] E_k`, `(alpha, gamma)` = (3, 5) or (5, 13) | PROVED |
| `1 + q 3^j 2^s` (`s >= 3`) and `1 + (2 q 3^j - 1) 2^s` (`s >= 1`) are path counts, `q` in {15, 25, 39, 65}; so are `2^s (3 2^m + 5 2^K 3^m) - 1` and `2^s (5 2^m + 13 2^K 3^m) - 1` | PROVED |
| Every odd `c <= 4097` has an anti-seed (ray `c 2^s - 1`) | FINITE-EXACT (2,049 certificates recounted) |
| Every odd `c <= 4097` has a seed (ray `1 + c 2^s`) | KNOWN (THM-5126 + `odd_c_seed_symmetry_20261110`, which found order-11 roots for the 33 values of HYP-9658 (B) first); here a second set of 2,049 certificates, the 33 by towers |
| A unit set lies in one strong component; all other components are single vertices | PROVED |
| One pair of order 8 has offset `-1` and slope `2^10`: tournaments on `9 + t` vertices with `2^(10+t) - 1` paths | FINITE-EXACT + KNOWN closing law (the numbers `2^k - 1` themselves are KNOWN: THM-4748) |
| Every odd `c` has a seed; Chan-Pak | OPEN |

Continuation of 2026-11-11 (source interfaces, the domination-defect inequality of Sessions 260-261, the
Euler-function criterion, two papers): [SOURCES_AND_DOMINATION.md](SOURCES_AND_DOMINATION.md).

Credit, after the sweep of 2026-11-11: while this session ran, other sessions pushed the strong-core structure of
unit sets (Session250, Session258), the twin identity, the pair bracket in one carrier and the slopes 15, 25, 39, 65
times `3^b` (Session260), the lines `b 3^a - 1` (endpoint_coverage_compiler) and order-eleven roots for the 33 lines
(odd_c_seed_symmetry). Rows of the table above that are also theirs: "a unit set lies in one strong component"
(KNOWN, theirs first) and the families with `q` (CONCURRENT). Not found elsewhere: the explicit transport matrix
(rows for offset, slope and end count), the sign flip for an arbitrary base (pre-registered 2026-11-10), the tower
recurrence, the certified list of anti-seeds.

What changed in understanding: a seed is not an isolated object found by census. Seeds and anti-seeds beget each
other through one operation that works over every base, and the only tail gadgets of order at most 7 that keep
`|r| = 1` for every base are governed by the two solutions of `2^delta - 3^(delta-1) = 1`.

## 1. Setting

As in THM-5126. For a pair `(B, O)`: `h = H(B)`, offset `r = h - W_B(O)`, slope `M = U_B(O)`, `E` = number of
Hamiltonian paths of `B` ending in `O`. The closed tournament `T_t(B, O) = (B => TT_t => v, v -> O)` has
`H = r + 2^t M` (KNOWN). `r = 1`: seed of the line `1 + c 2^s` (`c` = odd part of `M`); `r = -1`: anti-seed of the
line `c 2^s - 1`. Every odd `N > 1` lies on one line of each kind, and on one of them with exponent 1.

The second kind of line is forced on us by the first result below, and it fits the two forbidden values:
`7 = 1 + 3 * 2 = 2^3 - 1` and `21 = 1 + 5 * 4 = 11 * 2 - 1`. So no seed has slope 3, 5 or 10 and no anti-seed has
slope 1, 2, 4 or 11 (THM-1370 read on the rays).

## 2. The gadget law (THM-5130 (A))

Take a pair `(B, O)`, put a tournament `X` at the tail, add the closer `v` (`B => X => v`, `v -> O`), and choose a new
out-set `P` inside `X`. The new state is `G(X, P) (r, M, h, E)^T`:

    r+ = r_X r + kappa_X M,                    r_X = H(X) - W_X(P),  kappa_X = Pi_2(X) - 2 W_X(P) - 2 W3_X(P)
    M+ = u (r + h) + (4u - 2 sigma + 2 u3) M - sigma E
    h+ = H(X) r + Pi_2(X) M
    E+ = eps (r - h) + nu M.

Proof idea: contract `X` to a vertex `z`. A visit of `z` is entered from `B` and left towards `v`, and `v` occurs
once, so a path visits `z` at most twice and a pair of paths at most three times. Counting the walk systems by their
number of visits gives the coefficients `r + 2M`, `M` (one walk), `r + 2M`, `2M` (two walks glued at `z`),
`2M + h + r - E`, `4M + h + r`, `2M` (two walks, the second starting at `z`); the constants of `X` count the blocks.

Row `h` is Session216's star theorem and a case of the substitution law THM-3121 (4). The point of the full
matrix is that the *offset and the slope of the next closer* are transported too, so closers can be stacked.

## 3. Which gadgets keep a unit offset: 2^delta - 3^(delta-1) = 1

For `X = TT_a` and `P = {i, i + delta}`:

    r+ = rho_delta r + 2^(a-delta) (2^delta - 1 - 3^(delta-1)) M,     rho_1 = rho_2 = -1,  rho_delta = 1 - 2^(delta-1).

The bracket counts, with signs, the ways of distributing the `delta - 1` vertices between the two chosen ones:
`2^delta` from all covers by two paths, `3^(delta-1)` from the three-block covers. It vanishes exactly for
`delta = 1` (`2 - 1 = 1`) and `delta = 2` (`4 - 3 = 1`). At `delta = 3` it is `8 - 1 - 9 = -2`: Catalan's pair
`3^2 - 2^3 = 1` is the first distance at which the pair stops being neutral.

So a transitive pair at distance 1 or 2 turns every seed into an anti-seed and every anti-seed into a seed. The
exhaustive scan of all gadgets of order at most 7 (`neutral.py`) finds `kappa = 0` only for `|P| = 2`, always with
`r_X = -H(X)`: no gadget of order at most 7 keeps the sign for every base, and the transitive pairs at distance
1, 2 are the only ones that keep `|r| = 1` for every base.

For these two distances, with `alpha = 2^delta + 1` and `gamma = 3^delta + 2^delta` (the pairs (3, 5) and (5, 13)):

    r+ = -r,  h+ = r + 2^a M,  M+ = alpha 2^m (h + r) + 2 gamma 2^(i-1) 3^m M - [i=1] E,
    E+ = alpha 2^(i-1) M + [m=0] (r - h).

## 4. Towers and k-vertex extensions

Iterating gives a **tower**: a transitive chain `Z_0 => Z_1 => ... => Z_(k+1)` and closers `v_1, ..., v_(k+1)`, each
closer beating exactly two vertices of the previous tail (at distance 1 or 2), everything after its own tail, and the
later closers. It is a `(k+1)`-vertex extension of a transitive tournament.

Two facts answer the question about k-vertex extensions.

1. *Cover counts of one level are path counts of the next.* The slope `M_j` counts covers of the base `B_j` by two
   paths, the second starting in the out-set. The base of the next level is the closed tournament itself, so
   `H(B_(j+1)) = r_j + 2^a M_j`: the count that "qualifies" in the one-vertex extension of `B_j` becomes, exactly,
   the Hamiltonian path count of a tournament, and then feeds the next slope.
2. *The slopes satisfy a three-term recurrence.* Since `h_(j+1) + r_(j+1) = 2^a M_j`,

       M_(j+1) = alpha 2^(m + a_j) M_(j-1) + 2 gamma 2^(i-1) 3^m M_j - [i = 1] E_j.

   A tower slope is a continuant-like sum of terms `2^x 3^y` whose coefficients are products of 3, 5 and 13.

What this does and does not give for Chan-Pak. Every odd `N` is `1 + c 2^s` and `c' 2^s' - 1`; a seed for `c` with
small 2-adic valuation or an anti-seed for `c'` does the job. On `[3, 8191]` at least 3,722 of the 4,095 odd numbers
lie on a certified ray (`ray_cover.py`). But the rays of offset `+1` and `-1` are not everything: 15 is a path count
on no such ray (proof in THM-5130). So a proof along these lines needs all offsets, i.e. the full portrait
`{(M, r)}`; the gadget law applies to every offset (`r+ = -r` holds for all `r`), which is why it is stated for
states and not for seeds. The density form is HYP-9662 (B): if the odd slopes of seeds have positive lower density,
so does the spectrum.

## 5. Explicit families (THM-5130 (D))

- One closer (KNOWN, THM-4748): `3 * 2^n - 1` and `5 * 2^n - 1`.
- Two closers over a fan `(TT_n, {p})`:
  `H = 2^t [(2^delta + 1) 2^(m+1) + 2 (3^delta + 2^delta) 2^(i-1+n-p) 3^m - [i=1][p=n]] - 1`.
- Two closers over a transitive pair: `H = 1 + 2^t [2 q 2^(m_1+i_2-1) 3^(m_2) - [i_2=1][m_1=0]]`,
  `q = (2^delta_1 + 1)(3^delta_2 + 2^delta_2)` in {15, 39, 25, 65}. So the lines `c = 15 * 3^j, 25 * 3^j, 39 * 3^j,
  65 * 3^j` are in the spectrum from `s = 3` on, and the lines `c = 30 * 3^j - 1, 50 * 3^j - 1, 78 * 3^j - 1,
  130 * 3^j - 1` entirely.
- Depth three: the slope can be a power of two. The identity behind the smallest case is `5^3 + 3 = 2^7`
  (`M_3 = 2^9 (25 M_1 + 234)` with `M_1 = 646`): four closers, `N` vertices, `2^N - 1` paths for every `N >= 24`.

## 6. Seeds and anti-seeds for every odd c <= 4097 (THM-5130 (E))

`certs_plus.txt`, `certs_minus.txt`: one pair per odd `c <= 4097` and sign, least order found, each rebuilt and
recounted (`certs_verify.py` + `pairstate.c`, 0 mismatches of 4,098).

| | pairs of order <= 9 | one level | two levels | odd slope | largest order |
|---|---|---|---|---|---|
| seeds (`r = +1`) | 1,969 | 63 | 17 | 1,225 | 14 |
| anti-seeds (`r = -1`) | 1,959 | 88 | 2 | 1,187 | 14 |

In all, 2,016 seed and 2,013 anti-seed certificates have order at most 9 (47 and 54 of the towers do as well); the
other 33 and 36 have order 10 to 14.

The 33 exceptional values of HYP-9658 (B) (no seed among unit-offset roots of order at most 10) all have tower
seeds; `B+` has order 10 (2023, 3001), 11 (245, 2037, 3033, 3985), 12 (16 values), 13 (9), 14 (2). Credit: the
session of `../odd_c_seed_symmetry_20261110/` found unit-offset roots of order eleven for all 33 by a direct search
and pushed first; HYP-9658 (B) is settled there. The towers give the same lines by a different mechanism, at rooted
orders 11 to 15, and the anti-seed row of the table, as a certified list, is what is new here (single anti-lines
were known: `c` = 1, 3, 5 in THM-4748, 9 and 15 in the packets of 2026-11-11).

## 7. Structure of unit sets (THM-5130 (F))

- `U_B(o) = H(B) + sum over out-neighbours x of o of Q_B(x, o)`, so `M >= |O| H(B)`: for `|O| >= 2` the first
  tournament of a ray has at least `4 H(B) - 1` paths.
- Two vertices in different strong components have `W_ef >= H(B)`. So a unit set (`|O| >= 2`) lies in one strong
  component `Q`, every other component is a single vertex, and `M_B = 2^(vertices after Q) M_Q`. Seeds are
  *strong unit cores* with a transitive prefix and suffix. Anti-seeds are not: the pair of a transitive tournament
  has its two vertices in different components.
- Twin lemma: twins `e -> f` replacing a vertex `y` of `K` have `r = H(K[TT_2 at y]) - 2 H(K)`.
- `W_ef` does not see the orientation of the arc `ef`.
- Typical size: in a random-like tournament `W_ef` is about `4H/n`, so a unit set has about `sqrt(n/2)` vertices
  (heuristic). Up to order 9 only sizes 2 and 3 occur.

## 8. The eight numbers

Typed: **[T]** a role inside a theorem of this note, **[D]** a datum (certificate), **[P]** prior note, **[N]**
coincidence without mechanism.

| number | here | prior notes |
|---|---|---|
| 2 | [T] base of every ray `r + 2^t M`; two-path covers | |
| 3 | [T] `3 = 2^1 + 1` (first `alpha`); the third block of a cover gives the powers `3^m`; `2^delta - 3^(delta-1) = 1` | [P] `3^k` orders of the Sierpinski-Steiner-Collatz note (THM-4754..4757) |
| 5 | [T] `5 = 2^2 + 1 = 3 + 2` (second `alpha`, first `gamma`); anti-seed line `5 * 2^s - 1`; [D] its seed line starts at 41 because slopes 5 and 10 are forbidden by 21 | [P] genus-zero levels {2, 3, 5, 7, 13} (GENUS_BRIDGE) |
| 7 | [T] `7 = 2^3 - 1`, and `2^3 - 1 != 3^2`: the first distance where a transitive pair is not neutral; the missing first point of the line `1 + 3 * 2^s` and the missing point of the line `2^s - 1`; [D] anti-seed of slope `7 * 32` = one tower level over the seed of the line of 21 | [P] THM-1370, THM-4748 |
| 13 | [T] `13 = 3^2 + 2^2` (second `gamma`); `13 = 1 + 2 * 6`, first point of the ray of `C_3`; anti-seed slope `3 + 5 * 2` | [P] first one-apex miss (Session191); Ellison exceptions {13, 14, 16, 19, 27} (THM-5109) |
| 19 | [T] `19 = 5 * 2^2 - 1`: on the ray of the distance-2 pair of `TT_3` (9, 19, 39, ...); [D] seed of slope `19 * 2^7` (order 9); pure towers of depth 2 give slope `19 * 2^6` (run E1) | [P] orders `2 * 3^k + 1`; Ellison; THM-135 |
| 27 | [D] anti-seed of odd slope 27 (order 6); `27 = 3 + 2^3 * 3`; on no certified ray of offset `+1`, `-1` (it is `3 * 9`) | [P] Collatz root-27 / Antihydra (THM-4771) |
| 1039 | [D] `1039 = 130 * 2^3 - 1` (anti-seed of slope 130, order 7) `= -17 + 2^4 * 66` (Session195's ray) `= 15 + 2^3 * 128`; seed of odd slope 1039 (order 8), so `1 + 1039 * 2^s` is always a path count | [P] Session195; prime clocks `p_*(208) = 1039`; [N] 208 is the slope of the order-8 seed of `c = 13`; [N] `1040 = 2^4 * 5 * 13` |

Honest reading: 2, 3, 5 and 13 are structural here (they are the constants of the recursion), 7 is structural as
the failure of `delta = 3` and as the common root of the two lines `c = 3` and `c' = 1`, 19 is the second point of
an explicit ray, 27 and 1039 enter only through certificates.

## 9. Pre-registered runs

See `PREREGISTRATION.md` (pushed before each run). E1 (pure towers, depth <= 4): `N_2 = 805`, `N_4 = 1612` of 2,049.
E2 (towers over pairs of order <= 7, two levels): seeds for 1,872 of 2,049, among them 32 of the 33 exceptional values
(the 33rd, 3773, needs a pair of order 8; 3427 gets a small seed from a pair of order 9).

| | prediction | outcome |
|---|---|---|
| R1 | `300 <= N_4 <= 1400` | FAILED (1,612) |
| R2 | `N_4 >= 1.5 N_2` | CONFIRMED (1,612 against 1,207.5) |
| R3 | no power of two as slope at odd depth <= 3 | FAILED (15 states of depth 3; verified as tournaments) |
| R4 | at least 3 of the 33 exceptional `c` gain a seed | CONFIRMED (32) |
| R5 | least uncovered `c >= 3` is at most 13 | CONFIRMED (5) |
| R6a | at least 4 of 5, 11, 23, 73, 85, 91, 95, 127 at depth <= 6 (small box) | FAILED (3: 5, 73, 85) |
| R6b | 5 does not occur | FAILED (slope `5 * 2^32` at depth 6) |
| R6c | `N_6' >= 1500` | CONFIRMED (1,592; `N_2' = 614`, `N_4' = 1357`) |

Decision rule: `300 <= N_4 < 1900`, so the towers are reported as one mechanism among several. E4 (G1-G4, sources
and the domination ray): SOURCES_AND_DOMINATION.md section 8; G1 and G3 were not blind (audit).

## 9a. Audit

Independent audit: PASS WITH CORRECTIONS ([audit/AUDIT.md](audit/AUDIT.md), programmes and outputs in `audit/`).
All proofs of THM-5130 were checked; the gadget law (5,876 cases), all 4,098 certificates, the runs E1-E3, the neutral
gadgets of order at most 7, the grid sources and the 43-grid argument were recomputed with independent code. Majors:
(M1) F2 (a), (b) were typed "not found elsewhere" although Sessions 250, 258 and 216 have them; (M2) the credit for
the gadget law omitted ONE_EXIT.md and the marked tensor of Session261; (M3) the E4 pre-registration said that
nothing had been computed while two of its four predictions were decidable from data files of the same commit.
All majors and the 26 minors are applied; see MISTAKE-759 in `01-canon/MISTAKES.md`.

## 10. Connections (typed)

- **Catalan / Pillai.** Map: distance of the pair -> the unit equation `2^delta - 3^(delta-1) = 1`. Preserved:
  exactness (it is an identity of cover counts). Test: `verify_s22.py` V3 for `delta <= 5`.
- **Collatz-type maps.** On odd parts a tower level acts by `(a, b) -> (b, alpha a 2^y + gamma 3^m b)`,
  `(b, alpha a + gamma 3^m b 2^y)` or `(b, (alpha a + gamma 3^m b) / 2^v)`: two-term analogues of `x -> (3x+1)/2^v`.
  Analogy only; the 2-adic bookkeeping and the correction `E` are part of the state, and no dynamics is imported.
- **Sums of terms `2^x 3^y`** (Erdos-Lewin): which `c` are tower slopes is a representation problem of this type
  with structured exponents. Loss: coefficients 3, 5, 13 and the nesting; open (HYP-9662 (C)).
- **Session216 / Session197 / THM-3121**: the state space and linearity of `H` are theirs; **THM-4748**: one closer.
- **Session249** (Mahler source `R_*` of the conjectured spectrum): `R_*` asks for a ray from `s = 1` for every
  `c` other than 3 and 5; among the least-order certificates 1,225 of the 2,049 odd `c <= 4097` have a seed of odd slope.
- **Session195 A/K/D grammar**: its rays `66 * 2^k - 17` and `88 * 2^k - 13` are rays of pairs with offsets `-17`,
  `-13`; the forest grammar has no unit-offset state beyond fans (`exploratory/ozj.py`, slopes up to 65,536).

## 11. Search log (including dead ends)

1. Duplicating a vertex of a seed (`twinscan.c`): 596 twinnable pairs of order <= 8; the nontrivial ones duplicate a
   vertex of in-degree 1; families `6 + 7 * 2^a`, `16 + 24 * 2^a`, `4 + 42 * 2^a`. Explained by the substitution law,
   not by the towers.
2. Two-module families with offset `-1` (`mfan.py`): led to the second kind of line.
3. Twins of a transitive module: unit exactly when the module ray has offset `-1` -> twin lemma -> gadget law.
4. One-vertex operations that preserve a unit pair for every base (`unitops.py`, order <= 5): only sink and source.
5. Gadgets with a non-zero `kappa` and a Bezout condition `r_X r + kappa_X M = +-1`: none found with transitive `X`.
6. Tool failure caught by a control: `pairstate.c` is wrong at order 27 and valid only to order 24 (array size);
   a known transitive case exposed it. Certificates are limited to order 24 and have order at most 14.

## 12. Reproduction

    cd scripts
    python verify_s22.py                       # 18 blocks, under a minute, no input files
    python sources_check.py                    # S1-S9, no input files
    python mersenne_towers.py                  # depth-3 towers with slope a power of two
    gcc -O2 -o pairstate pairstate.c && python certs_verify.py ../certs_plus.txt && python certs_verify.py ../certs_minus.txt
    python ray_cover.py
    # class lists (classgen.c of 05-knowledge/results/first_creators_genealogy_20261107/): ./classgen 8 cls ; ./classgen 9 cls9
    python neutral.py 7                        # needs cls_1..7.txt in the working directory
    gcc -O2 -o portrait portrait.c ; ./portrait m 1 < cls_m.txt > pt_m.txt    # pairs with offset +-1
    gcc -O2 -o tower_reach tower_reach.c ; ./tower_reach 1 10 4 4 4097        # E1 ; mode 2 reads states "r g M E"
    python make_certs.py 1 4097 ; python make_certs.py -1 4097               # needs pt_1..9.txt
