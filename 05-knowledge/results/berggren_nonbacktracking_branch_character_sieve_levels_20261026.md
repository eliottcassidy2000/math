# Legal-path recovery, the non-backtracking Berggren tree, sieve levels and the branch character

2026-10-26, mac-mini-2026-10-25-pyth (worktree `codex/session-pyth-cuboids`).

**Status.**
- **PROVED, elementary:** THM-4807 and THM-4808.
- **CITED:** classical sieve theory and the concurrent proofs, identified in each section.
- **FINITE-EXACT and NUMERICAL:** HYP-9418. Nine pre-registered predictions (4 sieve ratios, 2 height sums,
  3 spectral statistics), all CONFIRMED.
- THM-4807 and THM-4808 are both independently audited, with all corrections applied. The THM-4808 corrections
  came after a pre-audit push and are logged as MISTAKE-678.

No statement here proves twin-prime infinitude, abc, or new Hamilton-spectrum gaps.

## 0. The request, the merge and the outcome

The owner asked for six things:
- continue from the twin-prime tree isomorphism;
- merge recent agent work;
- grow prime sieves past the fixed-modulus obstruction;
- push the F8 branch parity toward exact connections with abc, F1 and the 7/21 gaps;
- use the observation that an arc absent from every Hamilton path sits across a detectable vertex cut to recover the
  missing orientations;
- write a unifying account.

The same prompt reached at least four sessions. All of them landed on `origin/main` today, and this note credits each:
- guard_reconstruction;
- Session70 (alignment, guards, carriers);
- S61 (THM-4806, HYP-9417);
- growing_sieve_path and sieve_path_merge.

The unifying account is
[07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md](../../07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md).

## 1. Hamilton paths determine the tournament: the elementary route

[THM-4807](../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md)
proves all-order injectivity of `T ↦ HP(T)` from Rédei and Camion alone.

**The cut vertex is the tail.** Let `u → v` be an arc of a strong `T` that lies on no Hamilton path. Then:
- `T − u` is not strong, and `v` lies in its terminal component. A middle-component detour excludes every other
  position.
- `u → v` is the chord of a Hamilton path `uP'` from `u` to `v` that uses only traced arcs.
- No Hamilton path runs from `v` to `u`.
- Reversing `u → v`, together with any further untraced arcs, closes `uP'u` into a Hamiltonian cycle.

This yields three results:
- **An explicit separating word** for any two distinct tournaments, in three cases (α, β, γ). The γ case needs the
  least tail index and then the largest head index.
- **An endpoint decoder.** Inside a strong block, orient an untraced pair from the start to the end of a legal word.
  Exactly one direction qualifies.
- **An exact tail-deleted characterisation** of untraced arcs.

The four proofs compare as follows:
- The concurrent proofs (guard_reconstruction G1, Session70 G4/G5, growing_sieve_path P1–P4) all use the
  Busch–Jacobson–Reid separator `z` (Busch's thesis, Theorem 4.12). Their decoders use two-step counts or
  deleted-vertex reachability on the traced core.
- S61's single-flip route closes once its "strong form" (HYP-9417 (1b)) holds. That strong form is THM-4807 (b).
- All four decoders agree on every class with `n ≤ 9`, and on every strong class at `n = 10` (audit).
- This closes vitali-response C5's open question for the named guards. It does not close the fresh, unnamed
  re-marking questions of HYP-9416.

**Census, FINITE-EXACT** (the `n = 10` values come from the audit):

| n | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|
| untraced arcs in strong classes | 1 | 2 | 11 | 39 | 200 | 1207 |
| strong classes carrying one | 1 | 2 | 9 | 33 | 167 | 1072 |

The endpoint relation alone first collides at `n = 4`, so the decoder needs the traced core too. The traced core alone
suffices; its untraced orientations are read by any of the three non-local rules.

## 2. The Berggren tree is a non-backtracking walk

[THM-4808](../../01-canon/theorems/THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md)
(G) gives the structure.

**The operator.** On Euclid slopes `t = m/n`, write `U(t) = 2 − 1/t`, `A(t) = 2 + 1/t`, `D(t) = t + 2`. In theta-group
form these are `D∘S`, `D∘S∘ι` and `D`. Let `V = {t, −t, 1/t, −1/t}` and `τ(t) = 2 − t`. Then:
- the multiset `[Ut, At, Dt, τt]` is `D(V·t)`;
- `B = P_U + P_A + P_D` is exactly the Hashimoto non-backtracking operator of the graph on `V`-orbits.

**Over Q.** The graph has two components, exchanged by `C(t) = (t+1)/(t−1)` (THM-3357's intertwiner). The Euclid
component (slopes with `m + n` odd) is an edge-indexed tree whose covering tree is 4-regular:
- its vertices are the unsigned primitive triangles, plus the degenerate triangle `(1,0,1)`;
- its edges are the Berggren parent–child pairs, the edge `(1,0,1)–(3,4,5)`, and a half-edge at `(1,0,1)`;
- `(1,0,1)` is the root's fourth neighbour, and it is the pole `t = ∞` of the F8 marking.

Each primitive triple has three children and one parent. That ternary branching is exactly non-backtracking
continuation. Over Q this is, in other coordinates, Alperin's enumeration by reduced words in `Γ(2)` (Amer. Math.
Monthly 112 (2005)). The Berggren word is the even continued fraction of `m/n` (Kraaikamp–Lopes; Romik). What this
lane adds is the mod-`q` Hashimoto form and its spectra (audit: not found in the literature).

**Mod q.** The quotient graph `G_q` has about `q/4` vertices.
- **Ihara–Bass:** every eigenvalue `λ ≠ ±1` of `B` satisfies `λ² − μλ + 3 = 0`, with `μ` an adjacency eigenvalue.
- **Exact, all primes `q < 331`** (Sturm counts in the audit): every eigenvalue other than `3` and `±1` has
  modulus `√3`. The Ihara–Bass identity itself is PROVED: `BL = LA − T` and `BT = 3L`.
- **Ramanujan up to 1000:** `G_q` is Ramanujan except at `q = 331, 569, 719, 809, 887`.
- **Pre-registered range `1000 < q < 4000`:** 42 of 382 primes are exceptions (0.110), the worst ratio is 1.0105, and
  the largest nontrivial `|μ|` is 3.50. All three predictions were confirmed.
- **Uniform expansion** is CITED: Selberg's 3/16 theorem with the Brooks–Burger transfer (Lubotzky 1994, Ch. 4),
  passed to `G_q` by a Dirichlet-form comparison. The Ramanujan property is only nearly true. The exceptional shares
  (3% and 11%) lie far below the random-regular-graph limit of about 31% (Huang–McKenzie–Yau), so the statistics of
  `G_q` are open (HYP-9418 (4)).
- **The signed operator `B_χ`** has exact exceptions 251, 479, 863 and 947 below 1000 (audit).

**What this explains.** A sieve on the tree sees the tree through these quotients.
- **Fixed modulus.** Strong approximation and CRT make every finite quotient surjective. A fixed modulus therefore
  sees only congruence-forced predicates:
  - mod 3 forces the single edge of `J`;
  - a fixed modulus loses twin-ness (the CRT residue obstruction);
  - by THM-4808 (c), it also loses the branch character.
- **Growing modulus.** Expansion equidistributes depth-`k` residues in every class. That is the input an affine sieve
  needs (Bourgain–Gamburd–Sarnak; Kontorovich on Pythagorean orbits, CITED as a framework). On the D ray it reduces to
  Eratosthenes.

## 3. Growing sieves: the exact window, its density and what lies beyond

**The exact window.** growing_sieve_path W1–W3 (CITED) proves it:

    π₂(H) = R_⌊√(H+2)⌋(H) + π₂(⌊√(H+2)⌋).

- Prime-equality exceptions are restored.
- The square threshold is sharp.
- The compatible wheel path `x ≡ −1` represents the degenerate signed triple `(−1, 0, 1)`, the pole again.

**The density at the exact level.** The naive independent-residue density overcounts by

    ρ(X) = (X/2) ∏_{3≤q≤√X} (1 − 2/q) / π₂(X)  →  4e^(−2γ) = 1.2609…

This is Mertens' paradox in dimension two. All four pre-registered predictions were confirmed:

| X | ρ(X) | predicted range |
|---|---|---|
| 10⁸ | 1.11158 | [1.100, 1.130] |
| 10⁹ | 1.13123 | [1.120, 1.145] |
| 10¹⁰ | 1.14480 | [1.130, 1.160] |

The values also increase and stay below the limit, and the Hardy–Littlewood–Mertens model matches to within 0.001.
Our sieve to `10¹⁰` reproduces OEIS A007508.

So the sieve at its exact level decides twin-ness but cannot count twins from local independence: the
inclusion–exclusion error is as large as the main term. This is where the parity barrier lives.

**Beyond: classical transports, CITED.**
- **Brun.** Over the vertices of `J`, `Σ (1/r + 1/ρ) = B₂ ≈ 1.90216`. This converges, while the same sum over the whole
  D ray diverges.
- **Chen (1973).** There are infinitely many D-ray triangles with prime inradius and exradius a prime or a product of
  two primes ("Chen triangles"). A sieve can certify this much and, by the parity barrier, cannot go further.
- **Zhang, Maynard–Tao and Polymath 8b.** There are infinitely many pairs of prime-leg U-ray vertices `P(p)`, `P(p')`
  within ray distance `(p' − p)/2 ≤ 123`.

None of these decides whether `J` is infinite.

## 4. The F8 branch parity: exact connections and exact losses

THM-4808 parts (a)–(f), all PROVED:

- **(a) Pole frame.** `χ = αn − γm`, where `(α, γ)` is the image of the pole `(1, 0)`. The bit is the orientation of the
  frame formed by the transported degenerate triangle and the transported root.
- **(b) Automaton and closed form.** `#A` counts the even heads of a digit automaton on the regular continued
  fraction of `m/n`. In closed form, `χ = (−1)^(m + c₀ + ⋯ + c_L)` (audit, via THM-4057's Stern–Brocot depth
  cocycle). No residue character works, by (c).
- **(c) Every modulus.** `(16M+2, 4M+1) ≡ (2, 1) (mod M)` has word `A D^(M−1) U D` and `χ = −1`. So no residue datum,
  and no Dirichlet or Hecke character of finite conductor, determines `χ`. This extends the packet's mod-7 hostile to
  every modulus and fixes the Euclid vector itself, not just the conic point.
- **(d) Sign law, CITED: THM-3357 (19)–(20) and THM-3341 (24).** `sign(a − b)` sees only the trailing `A`-run. On the
  middle ray, `χ = −(a − b)`, the Pell norm sign in `Z[√2]`.
- **(e) Level identity, CITED: the scalar form of THM-3357 (37)/(40), via `U − A + D = I`.**
  `Σ_{|w|=k} χ x^#D y^#U = (x + y − 1)^k`. At every depth there is exactly one more even triple than odd. `χ` is
  THM-3357's `(−1)^#M` Walsh character.
- **(f) Cusp pairing.** Every non-root parabolic ray pairs with another of opposite character, with interleaved
  heights. Hence `Σ_{c≤X} χ = 1 + N_U + N_D + E(X)`, where `N_U + N_D = (1/2 + 1/√2)√X + O(1)` counts the two root
  cusps. `E(X)` oscillates:
  - `S/√X = 1.2421` at `2·10⁹` and `1.2194` at `5·10⁹`, both pre-registered and confirmed;
  - `1.2941` at `10¹⁰`;
  - the conjectured logarithmic mean is `1.2071` (HYP-9418 (2));
  - at 700 log-spaced heights in `[10⁴, 10¹¹]`, the remainder `E/√X` has logarithmic mean `−0.0043`. It is
    positive at 353 heights and negative at 345, changes sign 112 times, and has `max |E|/√X = 0.346`
    (`parity_height_logmean.c`, not pre-registered).

The main term `(1/2 + 1/√2)√X` of the parity bias comes from the two root rays. These are the D ray, carrying the
twin images and the centre identity, and the U ray, carrying the prime legs and the classical abc family
`P(3^(2^j))`. The rays the isomorphism `L(H) ≅ J` was built on are exactly the uncancelled cusps.

That the paired remainder `E(X)` contributes nothing on average is HYP-9418 (2), not a theorem. At full resolution
(audit, every `X ≤ 10⁹`): `max |E|/√X = 0.368`, and `E` changes sign 13,538 times.

**Verdicts on the requested connections.**

| Target | Exact content | Typed verdict |
|---|---|---|
| 7/21 gaps | Mod 7 with the pole fixed, `χ` is the converse bit of Paley 7 (order-21 group; finite_field §5, growing_sieve_path field). `H(Tᵒᵖ) = H(T)`, so no H-count statement, THM-1370 included, sees `χ`. | Exact non-connection; anything further is NUMEROLOGY |
| abc | Two abc families (the twin centre identity and the `9^(2^j)` family) lie on the `χ = +1` root rays, but they are not the only carriers. For the squared packet, `ADU = (16, 9)` gives `(175, 288, 337)`, with `c² = 113569 > rad = 70770` and `χ = −1` (audit; growing_sieve_path field §4 realises all four `(χ, hit)` combinations). `χ` does not determine abc-hit status, and neither does the radical determine `χ`. | Exact loss in both directions (this lane's pre-audit "χ = +1 on the abc carriers" was false; MISTAKE-678) |
| F1 | `χ` is the monoid evaluation at `(1, −1, 1)` (growing_sieve_path field §5). (c) shows that no fixed-level residue datum, so no finite-conductor multiplicative evaluation, sees it. | Exact loss: two free monoids, incompatible gradings |
| Twin primes | The twin images live on a `χ = +1` cusp, and the twin predicate needs the growing level `√x`. Both are invisible to every fixed modulus. `χ` is invisible at every level below the integer itself. | Shared obstruction type, different growth needs |

## 5. Merge record (read at `origin/main`, 2026-10-26)

**Incorporated or credited:**
- guard_reconstruction (G1–G3);
- Session70 alignment, guards and carriers;
- session71 status reconciliation;
- S61 THM-4806 and HYP-9417, including the strong form at `c71ab459f0`;
- growing_sieve_path and sieve_path_merge;
- unnamed-carrier (fresh re-marking; HYP-9416 stays OPEN beyond its universe).

**This lane's HYP-9416 audit** was applied: MISTAKE-676, 57 edits. That HYP was not C5's quotient, and its item (1) is
a trivial complete invariant.

## 6. Reproduction and audits

**Directories:**
- `04-computation/experiments/hamilton_path_recovery_20261026/`: THM-4807 scripts and the audit in `audit_thm4807/`;
- `04-computation/experiments/berggren_sieve_parity_20261026/`: THM-4808 and HYP-9418 scripts, plus
  `PREREGISTRATION.md`, committed at `2fe36c8088` and in a later commit, both before the runs.

**Runtimes on the mac-mini:**
- `verify_reconstruction.py` takes 13 s for `n ≤ 8` and 23 min for `n = 9`;
- `twin_sieve_ratio.c` takes 10 s to `10¹⁰`;
- `parity_height_big.c` takes 92 s to `10¹⁰`;
- `berggren_ramanujan_scan.py` takes 43 s for `1000 < q < 4000`;
- `berggren_ihara.gp` takes 6 s.

**Audits:**
- **THM-4807:** independent audit, confirmed at every order. One major and twelve minor corrections were applied
  before the first push; see `audit_thm4807/REPORT.md`.
- **THM-4808:** independent audit at `berggren_sieve_parity_20261026/audit_thm4808/REPORT.md`. Every proved part was
  confirmed, with 34 edits applied:
  - the graph over Q has two components, and the Euclid one is an edge-indexed tree;
  - prior art is now cited: THM-3357, THM-3341, Alperin, the even continued fraction;
  - the abc row is corrected;
  - the bias sentence is scoped;
  - Ihara–Bass and the box are now PROVED;
  - the closed form `χ = (−1)^(m + Σcᵢ)` is added.
  These are logged as MISTAKE-678.
