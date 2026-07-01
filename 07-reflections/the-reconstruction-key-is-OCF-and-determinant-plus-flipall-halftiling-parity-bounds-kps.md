# The reconstruction key is (OCF, determinant): the first OCF-cospectral twin at n=6 is det(I+S)-visible; and flip-all on the half-tiling gives the odd-parity lower bounds

*kind-pasteur-2026-07-01-S15. Two threads meet. (1) The twin separator: the n=6 reconstruction twins are resolved not by any single invariant but by the **pair** (OCF independence polynomial `I(Ω,x)`, determinant `d=det(I+S)`) — the two orthogonal coordinates of `the-determinant-lens`. The hard twin `{4,6}` is genuinely OCF-cospectral (same `I(Ω,x)`, score, `|Aut|`, cycle spectrum) yet `d`-visible. (2) The owner's identity "flipping all half-tiles = flipping all m tiles" makes the complement `φ` a fixed-point-free flip-all on the `2^D` half-tilings, and combined with the odd grid-symmetric count per SC class yields clean parity lower bounds.*

## The twin separator (which invariant splits the n=6 twins)

Six twin-pairs share `(category, blue-deg, black-deg, H)`. Testing invariants per twin node:

| invariant | pairs separated / 6 |
|---|---|
| `|Aut|` | **0** (all twins share `|Aut|`) |
| score sequence | 2 |
| `c₃` | 1 |
| `c₅` | 4 |
| `I(Ω,x)` (OCF poly) | 4 |
| **`d = det(I+S)`** | **5** |
| `cpA` (adjacency char-poly) | 5 |
| `cpS` (skew spectrum) | 5 |
| **any of the above** | **6/6 ✓** |

- `|Aut|` is useless (twins are automorphism-matched — mostly `|Aut|=1`).
- The OCF poly `I(Ω,x)` and `c₅` each get 4/6.
- **`d=det(I+S)` gets 5/6 — including the one pair `I(Ω,x)` misses.**

## The OCF-cospectral pair {4,6}: H-blind, d-visible

The pair `{4,6}` (both NS/pure-black, `blue-deg 0`, `black-deg 18`, `H=9`) is **indistinguishable by the entire OCF/cycle machinery**: identical score `(0,2,2,3,4,4)`, `|Aut|=1`, `c₃=3`, `c₅=1`, and the same independence polynomial `I(Ω,x)=1+4x`. This is the **first OCF-cospectral pair** — two non-isomorphic tournament classes the project's crown invariant cannot tell apart, arriving exactly at n=6 (the reconstruction break).

They are split by the **determinant coordinate**: `d(node 4)=64`, `d(node 6)=96` (and correspondingly the skew spectra `cpS` differ, `…39,0,9` vs `…55,0,25`). This is precisely the `the-determinant-lens` thesis in action: **`H` and `d` are two nearly-orthogonal smooth coordinates on the metagraph** (`R²≤0.004`), and where the path-count axis `H` is blind, the switching/determinant axis `d` sees. The reconstruction obstruction is not that the metagraph lacks a fingerprint — it is that **no single classical invariant suffices; the OCF and the determinant are jointly required.**

## The reconstruction key

`(I(Ω,x), d)` is a **complete node-fingerprint at n=6** (separates all 6 twin-pairs; the rest were already distinct). Since it is injective, knowing each node's `(I(Ω,x), d)` plus the connection pattern between fingerprint-values reconstructs the metagraph at n=6. So the answer to "what closes reconstruction beyond `(category, degree, H)`" is: **add the determinant `d`.** `H` alone stalls at the OCF-cospectral wall; `(H, d)` clears it — at least at n=6. (Whether `(I(Ω,x), d)` stays injective for `n≥7` is the natural next question; cospectral families may reappear, needing a third coordinate.)

## Flip-all on the half-tiling = the complement (the owner's identity)

The tiles split under the grid reflection `σ` into `f` fixed (anti-diagonal) + `p=(m−f)/2` swapped pairs; a **half-tiling** is a choice per `σ`-orbit, `D=f+p=⌊(n-1)²/4⌋=A002620` cells (THM-549/550), and grid-symmetric tilings = half-tilings. The owner's observation:

> **Flipping all `D` half-tiles = flipping all `m` tiles = the complement `φ`.**

Because flipping every `σ`-orbit (the `f` fixed cells and both of each of the `p` pairs) flips all `m` tiles. So on the grid-symmetric (blue) tilings, `φ` acts as **flip-all on the `2^D`-cube of half-tilings** — and flip-all is **fixed-point-free** for `D≥1`. Hence the `2^D` grid-symmetric tilings pair perfectly under `φ` into `2^{D-1}` **blue lines** (matching the count). The blue subgraph is thus the **complement structure of the half-tiling cube `Q_D`**, quotiented by iso — the exact recursion "blue = the pairing process one fold down."

## Odd grid-sym parity ⟹ lower bounds

Each SC class holds an **odd** number of grid-symmetric tilings (THM-281). Combined with flip-all being fixed-point-free:

- **`#SC` is even (clean re-proof via the half-tiling):** `2^D = Σ_{SC classes} g_C` with every `g_C` odd, so `#SC ≡ 2^D ≡ 0 (mod 2)` for `n≥3`. (Sharper than the `2^m` version: the sum is over the `2^D` half-tilings, and `#SC ≤ 2^D`.)
- **Blue-other-degree `≥ 1` (odd) for every SC node — no blue-isolated SC class.** A blue self-loop consumes 2 of a class's grid-sym tilings (a `φ`-pair inside the class); since `g_C` is odd, `g_C − 2·(#blue self-loops)` is odd `≥ 1`, so at least one grid-sym tiling has its `φ`-image in *another* class. **The odd count forces the SC spine to be blue-connected outward:** flip-all cannot close an odd set within one class. This is a structural minimum-degree bound on the blue spine, straight from parity.
- **`#pure-blue ≥ 1`** (the transitive class, `g=1`), and every pure-blue node's blue-degree equals its odd tiling count.

## Two twin mechanisms, cleanly separated

The two threads explain the two *kinds* of twin:
- **Black (NS) twins** — like `{4,6}` — have no grid-sym tilings (blue-deg 0); parity says nothing about them, and they need the **determinant `d`** (the switching axis) to separate. `d` is their reconstruction coordinate.
- **The SC/blue spine** is governed by the **odd grid-sym parity** (half-tiling flip-all): blue-degrees odd, spine blue-connected, `#SC` even. Mixed twins (e.g. `{12,22}`, `{43,44}`) sit here and are separated by `c₅`/`I(Ω,x)` and also `d`.

So reconstruction decomposes: the **blue spine** by parity/half-tiling recursion, the **black sea twins** by the determinant. `(I(Ω,x), d)` is the joint key.

## Honest status & next

- **Computed (n=6):** the twin-separation table; `{4,6}` OCF-cospectral, `d`-visible; `(I(Ω,x),d)` separates 6/6; the flip-all identity and the parity lower bounds (`#SC` even, blue-other-deg ≥1 odd) are proved.
- **Open:** does `(I(Ω,x), d)` stay injective for `n≥7` (or do larger cospectral families force a third coordinate)? Is the blue-spine minimum-degree bound improvable to *connectivity*? Does the half-tiling recursion make the blue counts satisfy THM-550's recurrences exactly?
- **Reframe:** "does `H` close reconstruction" → *no, but `(H, d)` does at n=6*; the metagraph's two orthogonal coordinates are jointly a fingerprint, and the OCF-cospectral wall is where you first need both.

— Related: `does-H-close-reconstruction-a-realization-degeneracy-metric-suite-kps.md` (the twins, the metric suite), `the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md` (`d ⊥ H`, the switching axis), `buckets-and-pairs-…` (R1 `#SC` even, R3 blue=half-tiling), THM-549/550 (half-tiling), THM-281 (SC sizes odd), THM-282 (blue=SC spine), HYP-3809 (mac-mini, SC twin-pairing). Script: `04-computation/merged_metagraph_twin_separator_kps.py` (+ .out). Not a HYP reservation.
