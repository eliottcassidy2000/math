# The seven-point sphere repaired, the Fano plane as a degenerate signature simplex, and equicentric tournaments

2026-11-07, mac-mini-2026-11-07-sevensphere. Owner directive (re-sent from 2026-10-21): "consider how tournament
structure relates to the open mathematical questions we investigate and see how the fano plane is analogous to the
below, and keep merging in a diverse set of ideas for inspiration on continued hypotheses", followed by the pasted
seven-point sphere of a tetrahedron.

Canon:
- [THM-4910](../../01-canon/theorems/THM-4910-the-seven-point-sphere-holds-with-the-anticomplementary-vertex-eight-points-on-a-skew-quadrilateral.md)
  (PROVED): the seven-point sphere holds when (C) is the anticomplementary vertex; eight points on a skew
  quadrilateral of centres.
- [THM-4911](../../01-canon/theorems/THM-4911-the-signature-simplex-of-a-tournament-equicentric-equilateral-skew-hadamard-and-the-fano-line-code.md)
  (PROVED): the signature simplex of a tournament; equicentric and equilateral tournaments; skew-Hadamard; the Fano
  line code.
- [HYP-9500](../hypotheses/HYP-9500-strong-equicentric-tournaments-exist-at-every-even-order-a-weak-skew-hadamard-shadow.md)
  (OPEN, large parts proved): strong equicentric tournaments at every even order.

Claims are typed PROVED, KNOWN, CITED, FINITE-EXACT, STRUCTURAL, ANALOGY, NUMEROLOGY, HEURISTIC or NEGATIVE.

## 0. What was already done (credited, not re-derived)

The same directive reached three lanes on 2026-10-21, and all three were read before any claim here.
- **session23b** (`tetrahedral_inversion_spheres_20261021/`):
  - read (C) literally as the vertex `V0` and proved the pasted claim false;
  - found the sphere `S1` through `B, C, D, E, F` and the exceptional locus;
  - gave a counterexample with an equilateral face.
- **opus-2026-10-21-S49**, THM-4742 (audited):
  - the six-label sphere `S2 ∋ A, B, D, E, F, H'`, the pencil with `S1`, and the inversion `ι` of radius `√(n+1)R`;
  - the defect `Σ = Σ|V0Vi|²` that keeps `V0` off `S2`, and the dimension range `n ≥ 3`;
  - the Fano comparison of the "seventh incidence".
- **The October Fano and tournament side.**
  - THM-4743 and Session24: converse Gluck.
  - session23b `fano_tournament_move_quotients_20261021`: the six edges of a tetrahedron plus a formal point form the
    Fano plane; Fano gauges.
  - mac-mini-2026-10-21-fanopsp: THM-4746/4747 Camion pseudoprimes, Paley/Kantor.
  - S49's §4 ledger of bridges.
- **Older.**
  - THM-1950 (klein; the inverse response `s(T)`, its bounds and the velocity law).
  - THM-211 (Paley design); THM-447 (skew-Sylvester double); THM-4557; THM-3729 (rooted Pfaffians).
  - Reid–Brown (doubly regular tournaments and skew-Hadamard matrices; classical).

This session works only on what the October lanes left open:
- the reading of (C);
- S49's open questions 2–3;
- a geometric home for tournament regularity in which the seven-point sphere and the Fano plane both appear.

No same-directive packets were pushed on 2026-11-07 before this session's claim (2c707d746b).

## 1. The seven-point sphere is true with the anticomplementary vertex (THM-4910)

**The reading.** The pasted (C) says "the vertex of the tetrahedron complementary to V1V2V3". October read it as
`V0`, and the claim fails because `Σ > 0`. Read instead as the vertex of the **anticomplementary** tetrahedron,

    C* = V0 + 4(G − V0) = 4G − 3V0 = V1 + V2 + V3 − 2V0,

the claim holds for every tetrahedron in the regime (`O, V0, G, O'` not coplanar), and in every dimension `n ≥ 3`
with `4 → n + 1`.
- This is the same "four times" construction as the pasted (D) and (F): `D` from `O` through `G`, `F` from `O`
  through `O'`, and now `C*` from `V0` through `G`.

**Why it is true.** Take the skew quadrilateral of centres `O → G → G' → O'` (circumcentre, centroid, facet centroid,
facet circumcentre).
- THM-4742's sphere `S2` is the sphere whose power at each corner is `c·(R² − |OP|²)`, with `c = 4, 3, 2, 3` for a
  tetrahedron (`n+1, n, n−1, n` in general).
- Vieta on each side gives the two points of `S2` on it:

  | side | points | meaning |
  |---|---|---|
  | `OG` | `D, E` | the simplex's ×(n+1) point and the inverse of `G` |
  | `OO'` | `F, B` | ×(n+1) of `O'` and its inverse (the pole of the face) |
  | `GG'` (the median) | `C*, A` | ×(n+1) from `V0`, and the circumsphere chord |
  | `O'G'` (face Euler line) | `H', X` | the face's own ×n point (orthocentre) and the inverse of `G'` in the face circle |

- The pasted seven are these eight minus `X`.
- **The defect read again.** `V0` lies on the median, outside `S2`. Its power is `V0A · V0C* = Σ_i |V0Vi|²`, which
  is exactly THM-4742's defect. The literal (C) is "the wrong end of the chord".

**Checks** (exact; `eight_point_sphere_exact.out`):
- symbolic for all `n` in THM-4742's Gram basis;
- 150 rational tetrahedra, and 40/20/10 rational simplices for `n = 4/5/6`.

In every exact case the eight labels are distinct and lie on `S2`, the corner powers are as stated, and the defect
matches.

**What does not extend** (NEGATIVE, `flag_negative.out`). The corner-power rule does not continue down a flag of
faces for `n ≥ 4` (residuals 0.12–0.19). The eight-point sphere is a one-facet phenomenon.
- S49's question 2 (the `2(n+1)` spheres through `D, E`): a numerical scan of ten candidate coincidences found none
  beyond the shared points `D, E` and the coplanar centres (all in the bisector plane of `DE`). It stays open.

## 2. The Fano plane: complementary against anticomplementary, and the Fano tournament

- **Two homotheties about the centroid** (ANALOGY with an exact core).
  - Fano's real picture is a triangle, its centroid, and its *complementary* (medial) triangle. That triangle is the
    homothety of ratio `−1/2` about `G`. The seventh line, through the midpoints, is drawn as the circumcircle of the
    complementary triangle, i.e. the nine-point circle. It is straight only in characteristic 2.
  - The seven-point sphere's seventh point is a vertex of the *anticomplementary* simplex, the homothety of ratio `−n`
    about `G`.
- **Mod 2 the two readings agree** (STRUCTURAL, THM-4910(5)).
  - In barycentric coordinates `C* = −2V0 + V1 + V2 + V3 ≡ V1 + V2 + V3 (mod 2)`. That is the PG(3,2) point `{1,2,3}`,
    the complement of `{0}` in the subset model of PG(3,2).
  - So the pasted phrase "complementary to V1V2V3" names the right point over GF(2), the field where Fano's seventh
    line is straight. Over `R` it must be read as anticomplementary.
  - That this explains the owner's wording is an interpretation, not a claim.
- **The Fano plane is a signature simplex** (PROVED, THM-4911(5)).
  - For the Paley tournament QR7, the signature vector `p_i = (I+M)e_i` of each vertex is `−1` exactly on the Fano
    line `i + {1,2,4}`.
  - The seven signature vectors are the incidence vectors of the seven Fano lines, pairwise at Hamming distance 4,
    forming a regular 6-simplex in the 7-cube.
- **The irony** (STRUCTURAL). The seven-point sphere needs separated centres (`O ≠ G`, non-coplanar). QR7's simplex
  is equilateral, so all its centres coincide. The Fano tournament is therefore exactly where the seven-point sphere
  does not exist. More generally the sphere sees precisely a tournament's deviation from regularity.

## 3. Tournament structure and an open problem (THM-4911, HYP-9500)

**The signature simplex.** Send vertex `i` to its signature vector `p_i ∈ {±1}^n`: `+1` on itself and its
in-neighbours, `−1` on its out-neighbours. All `p_i` lie on the sphere of radius `√n`.
- **The inverse response as a distance.** THM-1950's `s(T) = uᵀ(I+MᵀM)⁻¹u` (also `F_T(1)` of THM-4900) is
  `1/dist(0, aff{p_i})²`, and the circumradius is `r_T² = n − 1/s(T)`.
  - With THM-4900(3), `s(DT) = s(T)`: the skew-Sylvester double keeps the signature simplex at the same distance from
    the origin.
  - Doubling also preserves skew-Hadamard (KNOWN: THM-447(2)), since `I + D(M) = [[H, H], [−Hᵀ, Hᵀ]]` with `H = I + M`.
- **Edges are separation counts.** `|p_i − p_j|² = 4(1 + δ_ij)`, where `δ_ij` counts the vertices beating exactly one
  of `i, j`.
- **The regularity hierarchy is a symmetry hierarchy.**
  - Equicentric (centroid = circumcentre) means `M²u ∥ u`.
    - For odd `n` this is exactly regular, by a Pfaffian kernel argument: every signed sub-Pfaffian of a `±1` skew
      matrix is odd.
    - Non-strong equicentric tournaments are exactly `R_a ⇒ R_b`.
    - Almost-regular ones exist exactly when `n ≡ 2 (mod 4)`.
    - `λ ≡ n − 1 (mod 4)`.
  - Equilateral is exactly doubly regular (`c = −1`) or skew-Hadamard (`c = 0`). The proof is a two-line eigenvalue
    argument.
- **The open problem.** By Reid–Brown, the skew-Hadamard conjecture (OPEN) says that every `n ≡ 0 (mod 4)` carries a
  tournament with pairwise orthogonal signature vectors.
- **HYP-9500**, a weak shadow. Every even order `n ≥ 6` should carry a *strong equicentric* tournament (orthogonality
  replaced by the single eigenvector condition), and the `λ`-values should fill a whole progression. Proved so far:
  - `n ≡ 2 (mod 4)`;
  - every `n ≡ 0 (mod 4)` with a skew-Hadamard matrix (switch a bordered doubly regular tournament at one vertex).

  Open at the remaining orders.
- **The seven-point sphere of a tournament** (THM-4911(6)). For every non-equicentric tournament and vertex `v` in the
  regime, the eight-point sphere of the signature simplex has `Pow(O_T) = n(n − 1/s(T))` and `Pow(p_v) = 4(n − 1 + Δ_v)`.
  Here `Δ_v` counts the directed 2-paths ending or starting at `v`.

## 4. Merging ideas: a ledger for continued hypotheses

| idea | status | what transfers |
|---|---|---|
| Seventh incidences. Fano's line needs char 2. The sphere's literal `V0` needs `Σ = 0` (impossible over `R`); its anticomplementary `C*` holds identically. Positive Collatz cycles stay open. | ANALOGY (S49's ledger, updated) | Find the right point before calling an incidence obstructed. In two of the three cases the obstruction was a misreading or a characteristic. |
| Corner weights `c = n+1−α−β` on the Boolean square (centroid?, facet?) | STRUCTURAL | A third axis (anticomplement: `h_(G,−n)` maps `O, G, G', O'` to `D, G, V0, O'_W`) gives a cube of centres. The two squares do not share one sphere (`D ∈ S2` while `Pow_(S2^W)(D) = (n+1)n²R² ≠ 0`). |
| Signature simplex and THM-1950's `s` | PROVED | THM-1950's velocity addition under `⇒` becomes a rule for composing distances: with `d² = 1/s`, `d²(A ⇒ B) = (d_A²d_B² + 1)/(d_A² + d_B²)`, the coth-addition law. |
| Equicentric ⊂ equilateral ⊂ skew-Hadamard (as conditions on `I + MᵀM`) | PROVED + OPEN | HYP-9500 is a tournament-only step towards the skew-Hadamard conjecture. |
| QR7 = Fano line code; the regular simplex | PROVED | Doubly regular tournaments are the designs whose signature simplex is regular (classical). |
| `Pow(p_v) = 4(n − 1 + Δ_v)` with 2-path counts | PROVED | A new tournament invariant attached to each vertex by the sphere. Its relation to `H(T)` is untested (HEURISTIC). |
| The 7/8 counts (seven labels plus `X`; Fano's 7 points and PG(1,7)'s 8, both acted on by the group of order 168) | NUMEROLOGY | No mechanism; recorded so it is not re-derived as structure. |

## 5. Ledger

| claim | type |
|---|---|
| `C* ∈ S2`: the pasted seven-point sphere holds with the anticomplementary vertex, all `n ≥ 3` | PROVED (new) |
| Eighth point `X`; corner powers `c·(R² − OP²)`; side-pair structure; defect `= V0A·V0C*` | PROVED (new) |
| `S2` through `A, B, D, E, F, H'`; the `ι`-invariance; `Σ`; all `n`; the literal reading false | KNOWN (THM-4742, session23b) |
| Mod-2 reduction of `C*` to `{1,2,3}` | STRUCTURAL |
| No flag extension (`n = 4, 5, 6`) | NEGATIVE (numerical) |
| Signature simplex: cube, sphere, `δ`-edges, `s = 1/dist²`, `r² = n − 1/s` | PROVED (new in this form) |
| Equicentric: odd `n` ⟺ regular; non-strong ⟺ `R_a ⇒ R_b`; almost regular ⟺ `n ≡ 2 (mod 4)`; `λ ≡ n−1 (mod 4)` | PROVED (new) |
| Equilateral ⟺ `c ∈ {−1, 0}` (doubly regular or skew-Hadamard) | PROVED; KNOWN in substance (Reid–Brown) |
| QR7's signature vectors are the Fano lines | PROVED; KNOWN in substance (THM-211) |
| Eight-point sphere of the signature simplex; `Pow(p_v) = 4(n−1+Δ_v)`; `Pow(O_T) = n(n − 1/s)` | PROVED |
| HYP-9500 | OPEN (proved for `n ≡ 2 mod 4` and for skew-Hadamard orders) |
| Equicentric counts 1, 2, 1, 5, 3, 14, 15, 165 (`n = 3..10`) | FINITE-EXACT |

## 6. Next obligations

- **HYP-9500.**
  - A direct construction of strong equicentric tournaments at `n ≡ 0 (mod 4)` (blocks with `x ∈ {±1, ±3}`).
  - Local search at `n = 12, 16, 20` for the full `λ`-progression.
- **The sphere as an invariant.** Does `Pow(p_v)`, or the radius of `S2(T, v)`, separate tournaments that `s(T)`
  does not? Does it correlate with `H(T)`?
- **S49's question 2 for signature simplices.** The `n` vertex choices of one tournament give `n` eight-point spheres
  through the common `D_T, E_T`. Is their configuration a tournament invariant?
- **The two homotheties.** Is there a single statement covering both the complementary circle (Fano's seventh line,
  the nine-point circle) and the anticomplementary sphere?

## 7. Files

`04-computation/experiments/sevensphere_anticomplementary_20261107/`:
- `eight_point_sphere_exact.py`: THM-4910 (P1)–(P3), exact.
- `flag_negative.py`: THM-4910(6).
- `signature_simplex.py`: THM-4911(1)–(6); census through `n = 10`, about 2 minutes; needs nauty `gentourng`.
- `switching_paley.py`: HYP-9500's switching construction.
