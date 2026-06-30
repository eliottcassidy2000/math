# APEX CORE ATLAS — a field guide to the Z_p cores, their gaps, and their graphs

*klein-2026-06-29-S21. A finely-detailed reference map of the apex problemscape: the cores O subset Z_p, the cyclotomic gap g(O) (THM-590), and every clean structural feature found. For myself and future agents. Scripts: apex_core_atlas / apex_cross_prime_family / apex_core_graph_map_and_profiles / apex_arc_fejer_kernel _klein.py.*

The central object: for a core `O subset Z_p`, the **apex gap**
`g(O) = min_{k != 0} |sum_{x in O} zeta^{kx}|^2` (`zeta = e^{2pi i/p}`).
By THM-590 / HYP-3604, `g(O)` is exactly the **smallest eigenvalue of the autocorrelation circulant**
`C(O)_{ij} = a((i-j) mod p)`, `a(d) = #{x in O : x+d in O}` — a Bochner-PSD (SOS, sigma-even) matrix.

---

## I. The five gap classes in Z_7 (exhaustive, all 127 nonempty cores)

| gap value | name | # cores | sizes | flat? | circulant graph |
|---|---|---|---|---|---|
| `0` | ZERO / **cusp** | 1 | {7} | yes | `7J` (rank-1, all-resonance) |
| `4cos^2(3pi/7)=0.198062` | **DOUBLET** (binding) | 42 | {2,5} | no | `2I+A(C_7)=Q(C_7)` signless Laplacian of the **7-cycle** |
| `0.307979` | mid (arc) | 42 | {3,4} | no | interval autocorr `[3,2,1,0,0,1,2]` |
| `1` | unit | 14 | {1,6} | yes | `I` (empty graph) |
| `2` | FLAT / difference-set | 28 | {3,4} | yes | `2I+(J-I)=I+J` ≈ **K_7** (Paley/Fano) |

Multiplicities `{2,42,42,14,28}` sum to `128 = 2^7` (incl. `empty`).

---

## II. The classification rules (fine patterns)

- **[R-flat] Flat spectrum <=> INTEGER gap.** A core has a flat nonzero spectrum (`|Ohat(k)|^2` constant
  over `k != 0`) iff it is a perfect difference set, a singleton, or the full set — and exactly these have
  gap in `{0, 1, 2}` (integers). The NON-flat cores (doublets, generic triples) have the IRRATIONAL gaps
  `{0.198, 0.308}` in `Q(cos 2pi/7)`. **Integer gap = flat = difference-set-like; irrational gap = the
  binding non-flat cores.**
- **[R-diffset] gap = (p+1)/4 (=2 at p=7) <=> perfect (p,(p-1)/2,*) difference set** (QR/QNR & translates;
  the (7,3,1) Fano planes). These are the 28 cores of gap 2: 14 triples + 14 quadruples. Equivalently the
  circulant graph is `K_7` (`I+J`): maximally spread.
- **[R-spread] The gap ORDERS cores by SPREAD (it is a concentration index).**
  LOW gap = CONCENTRATED core; HIGH gap = SPREAD core:
  `doublet (0.198) < arc-3/4 (0.308) < singleton (1) < difference-set (2)`, with full-`Z_7` = `0` (degenerate).
  **The binding apex obstruction is the LEAST-spread core (the doublet) — the OPPOSITE of the random/Paley
  core.** The hardest case is maximal concentration, not maximal equidistribution.
- **[R-doublet] every doublet `{a, a+d}` (any difference `d != 0`) gives the SAME gap `4cos^2(3pi/7)`;**
  the worst Fourier mode is `k*` with `d k* ≈ ±(p±1)/2` (for `d=1`, `k*=3,4`); `|1+zeta^3|^2 = 2+2cos(6pi/7)`.

---

## III. Cores are circulant graphs; ARCS are the Fejér–Bochner cores

The autocorrelation `a(d)` makes `C(O)` a (weighted) **circulant / Cayley graph** on `Z_p`, and `g(O)` is its
signless-Laplacian-type least eigenvalue. Named cases:
- **doublet -> 7-cycle `C_7`**: `g = lambda_min(Q(C_7)) = 2 - 2cos(pi/7) = 4 sin^2(pi/14)`. Positive **iff
  C_p non-bipartite iff p ODD** (HYP-3604/3606). This is why the apex obstruction is positive: 7 is odd.
- **difference set -> `K_7` (`I+J`)**: `g = (p+1)/4`. Maximally spread.
- **singleton -> empty graph `I`**: `g = 1`. **full `Z_7` -> `7J`** (rank-1): `g = 0` (the cusp).
- **m-ARC `{0,..,m-1}` -> Dirichlet/Fejér core**: `|Ohat(k)|^2 = (sin(pi m k/p)/sin(pi k/p))^2` = the
  **Fejér kernel** (verified p=7,11,13). So the arc cores ARE the floor owners' Fejér–Bochner minorant
  cores (S75e). The **doublet is the 2-arc**, and its gap `4 sin^2(pi/2p)` is the Fejér kernel's minimum.
  The binding apex atom = the smallest Fejér-kernel value = the 2-arc.

---

## IV. The cross-prime family (apex prime p)

- **Doublet (binding) gap `= 2 - 2cos(pi/p) = 4 sin^2(pi/2p) ~ pi^2/p^2`.** Values:
  `p=3:1, 5:0.382, 7:0.198, 11:0.081, 13:0.058, 17:0.034, 19:0.027`. **Smaller odd apex prime = larger
  floor**; the obstruction decays like `pi^2/p^2`. `p=7` (LRC(14)) is the 3rd-largest.
- **p mod 4 split (Paley / the nu_2 theme).** For `p = 3 mod 4` (3,7,11,19) the QR is a **Paley difference
  set** => FLAT spectrum, gap `(p+1)/4`. For `p = 1 mod 4` (5,13,17) `-1` is a QR, the QR is symmetric and
  is NOT a difference set (2 spectral values). **`7 = 3 mod 4` is the Paley/flat-QR/difference-set case**
  (matches CLAUDE.md `nu_2 = 0 <=> apex = 3 mod 4 <=> Paley`).
- **# distinct gap values grows**: `1, 2, 4, 14, 35, ...` for `p = 3,5,7,11,13` (apex spectral complexity).
- **Bracket**: the gap spectrum lives in `[4 sin^2(pi/2p) (doublet, min), (p+1)/4 (Paley diff-set, max)]`
  for `p = 3 mod 4`.

---

## V. The descent GAP-PROFILE (a new trackable)

For a covering `S`, the 2-adic descent (THM-580) yields cores `O_0, O_1, ...`; the **gap-profile** is the
sequence `[g(O_0), g(O_1), ...]`. Computed:
```
consec{1..13}        [0.0,  0.308, 0.198, 1.0]    min-pos 0.198, depth 4
tightest{1..12,182}  [1.0,  0.308, 0.198, 1.0]    min-pos 0.198, depth 4
skip12{1..11,13,84}  [0.0,  0.308, 0.198, 1.0]    min-pos 0.198, depth 4
binding {1..13}\7    [1.0,  0.308, 0.198, 1.0]    min-pos 0.198, depth 4
consec{1..7}         [      0.308, 0.198, 1.0]    min-pos 0.198, depth 3
```
- **Near-universal TAIL `[..., 0.308 (arc), 0.198 (doublet), 1.0 (singleton)]`**: the descent always ends
  in the small odd cores `{1,3,5}->{1,3}->{1}`, so the tail is fixed; the binding is the doublet at the
  penultimate level.
- **The HEAD (level 0) distinguishes coverings**: `0.0` = the dense-consecutive coverings whose odd part
  fills `Z_7` (pass through the **cusp**, mod-7 covering at the top); `1.0` = the co-singleton coverings.
- The TOTAL floor is the product over the profile (`-> 0` for deep chains, klein-S16); the per-level binding
  atom is universally the doublet `0.198`.

---

## VI. New trackables (defined here)

1. **gap-profile** `[g(O_j)]_j` — per-covering descent signature (universal tail, varying head). §V.
2. **concentration/spread index** — `g(O)` itself orders cores; low = concentrated (binding), high = spread.
   Normalized: `g(O) / ((p+1)/4)` in `[0,1]` for `p=3 mod 4`. §II[R-spread].
3. **core shape type** — doublet / arc / V / difference-set / singleton / full (the structural class). §III.
4. **flat-vs-nonflat** — integer gap (flat, difference-set-like) vs irrational gap (binding non-flat). §II.
5. **level-0 apex type** — cusp (`g(O_0)=0`, dense fills `Z_7`) vs proper (`g(O_0)>0`). §V.

---

## VII. The grader space — multi-axis, not one line (klein-S22, HYP-3610)

The gap `g` is only ONE grader. A family of spread/concentration metrics splits the apex into `>= 3`
INDEPENDENT axes (Spearman over all 127 cores):
- **Axis 1 = `g` (gap)** — worst-mode concentration / the certificate. Concentrated pole = **doublet** (odd
  cycle `C_p`); spread pole = difference set.
- **Axis 2 = `F` (flatness `GM/AM`) ~ -`D` (difference-set DEFECT `Var(a(d))`)** — global difference-set-ness.
  Only **0.50** correlated with `g`; concentrated pole = the **ARC/interval** (a DIFFERENT core than `g`'s).
- **Axis 3 = `IPR` ~ entropy ~ odd-girth ~ -`|O|`** — effective support / size.

Disagreement (the proof it is multi-axis): `doublet` minimizes `g` (0.198) but `arc-4` minimizes `F` (0.5,
defect 0.667). The difference set is the common spread/random pole. CONJECTURE: `g` = the existence/
certificate axis (sigma-odd), `F` = the equidistribution/measure axis (sigma-even).

**Universal graders** (need only a spectrum or an associated graph -> grade every object):
- `F` (flatness, any spectrum) — structured<->random.
- `g` (gap / least eigenvalue, any spectrum) — certificate / worst mode.
- **odd-girth / cyclicity** (any associated graph) — intransitivity / cycle-scale. VERIFIED cross-object:
  on tournaments the n=4 Klein-four classes grade `{inf,3,3,3}` (transitive = orderable `inf` cusp = the
  `g=0` analog; rest = Condorcet `3`); cyclicity `{0,1,1,2}`. The apex cusp (full `Z_7`, `g=0`) and the
  transitive tournament (orderable, odd-girth `inf`) are the SAME degenerate pole. New invented metrics:
  difference-set defect, Cayley odd-girth, binding-mode `k*`.

## The one-paragraph picture

The apex problemscape is the lattice of `Z_p` cores graded by a single concentration index `g(O)` = the
least eigenvalue of the core's autocorrelation circulant graph. The spectrum is bracketed by the
**doublet** (2-arc, the smallest Fejér-kernel value `4 sin^2(pi/2p)`, the binding obstruction, the
signless-Laplacian of the odd cycle `C_p`, positive iff `p` odd) at the concentrated end, and the **Paley
difference set** (`K_p`, gap `(p+1)/4`, flat, exists iff `p = 3 mod 4`) at the spread end, with singletons
(`I`, gap 1) and the full-set cusp (`7J`, gap 0) as the flat-integer poles. The 2-adic descent threads
every covering through this lattice with a near-universal gap-profile tail ending at the binding doublet.

See THM-590, HYP-3604/3606 (the least-eigenvalue certificate), HYP-3598/3600 (the finite families),
HYP-3612 (the gap localizes to level 0), HYP-3611 (this atlas's new patterns/trackables).
