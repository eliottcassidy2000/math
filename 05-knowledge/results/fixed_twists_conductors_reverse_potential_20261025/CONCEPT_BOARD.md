# Concept board (second session half)

Portfolio: **Anchor** — the normseven packet's two stated targets (fixed-twist residual data,
conductor classes on non-UFD bases). **Niche** — implementing the reverse valuation potential on the
production counting graph. **Wildcard** — Antihydra as a valuation barrier, with the LCA reading.

| Concept | Object / operation | Preserved | Lost / sidecar | Decisive test |
|---|---|---|---|---|
| Twisted root algebra | `R[u]/(u_i^2 - d_i L_i)`, normalize | characters, valuations | normality; keep `c_z` per character | `(uv/3)^2 = 35` over `Z[1/2]` |
| Port classes | `c_z^{-1} u_z` over Dedekind `R` | Steinitz class per character | principal generators | `Z[√-21]`: `u√-21/3` integral, `u/3` not |
| Compatibility bank | single-root tests, mixed-root tests | `B`-membership after every continuation | the carry (3) is additive, not monoidal | scalar `a` needs `15 | a` for `(15,21,55)` |
| Fixed-twist lifts | squareclass vector at a base point | the three pair ports, the conic | interior points | Hilbert symbols vs brute force, 300 twists |
| Reverse potential | suffix minima on the ordinary graph | every continuation, `B = floor(N/2)` | `v_2(a(N))` (carries) | naive initial subtraction has a negative edge |
| Fixed precision | edges of reduced cost `< K` | `Z_terminal mod 2^K` | validity for larger `K` | `K=1`: 550 states at `N=100` |
| Antihydra barrier | walk `W_k = sum (2v_i-3)` | halting ⟺ `W = -6` | actual orbit vs model | `phi^{-6}` is model-only |
| Valuation LCA | trie depth `1 + sum v_i` | cylinders are balls | basin membership | Kempe half-length `v_2(y+1)` |

Connection contracts filled: conductor valuation matrix → BOUNDARY (4)–(6) (map: set all `c_z = R`);
central conic → THM-4765 Walsh characters (map: evaluate the three pair characters on the twist);
reverse potential → REVERSE_POTENTIAL (map: identity on the graph; verified `B`);
Antihydra → THM-2228 (map: negation); cylinders → Session56 trie (map: depth `= 1 + sum v`).
Stopping reasons: the six low-order ledger cells (framework not re-implemented); Petersen beyond
THM-4529/261/071 (no map); the mixed-order collisions of HYP-9374 (untouched this half).
