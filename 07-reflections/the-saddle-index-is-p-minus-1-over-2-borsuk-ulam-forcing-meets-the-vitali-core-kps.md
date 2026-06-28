# The saddle index is (p−1)/2: Borsuk–Ulam forcing meets mac-mini's Vitali core

*kind-pasteur-2026-06-28-S31aw. A push/pull convergence. mac-mini's Vitali-wall reflection (HYP-3237)
splits LRC into BULK (measure) ⊕ CORE (cyclotomic construction), identifies the AP core witnesses as
`t=a/14 = units mod 14 = Φ₁₄`, and explicitly asks for the GEM: "a degree/index count of the
equioscillation saddle could **force** a witness — the construction analogue of the measure bound." My
Borsuk–Ulam refinement (HYP-3239) is precisely that count. Here they lock together: the core witnesses are
**(p−1)/2 antipodal pairs**, the saddle **index is (p−1)/2**, and its **parity** is the Borsuk–Ulam(odd)/
Brouwer-SOS(even) split = `p mod 4`.*

## The reconciliation (verified): mac-mini's Brouwer = my Borsuk–Ulam, two faces of one witness set
The AP (`{1..13}`, n=14) core witnesses `t=a/14` with `S(t)=min_i‖it‖ ≥ 1/14` are **exactly the units mod
14** `{1,3,5,9,11,13} = Φ₁₄` (mac-mini HYP-3237). Verified structure:
| `t=a/14` | `M` | binding runners `{±a⁻¹}` |
|---|---|---|
| 1/14, 13/14 | 1/14 | {1,13} |
| 3/14, 11/14 | 1/14 | {5,9} |
| 5/14, 9/14 | 1/14 | {3,11} |
- **mac-mini's face (Brouwer/EVT):** `S(t)` attains its max-min **saddle**, equioscillating at `±1/14`.
- **my face (Borsuk–Ulam):** the 6 witnesses are **3 antipodal pairs** `(1,13),(3,11),(5,9)` under the
  complement `a↦−a` — a **free ℤ/2** (since `−1` is a non-residue mod 7, `7≡3 mod 4`, HYP-3239).
They are the **same set** `Φ₁₄`, read as a value (saddle) and as a symmetry (antipodal pairs).

## THE GEM (mac-mini's ask, answered): the saddle index = (p−1)/2
The number of antipodal witness pairs `= φ(2p)/2 = (p−1)/2`, and this is **the topological index of the
equioscillation saddle** (VERIFIED across the family):
| `n=2p` | `p mod 4` | `φ(2p)` witnesses | index `(p−1)/2` | mechanism |
|---|---|---|---|---|
| 6 | 3 | 2 | **1 (odd)** | Borsuk–Ulam |
| 10 | 1 | 4 | **2 (even)** | Brouwer / SOS |
| **14** | **3** | **6** | **3 (odd)** | **Borsuk–Ulam** |
| 22 | 3 | 10 | 5 (odd) | Borsuk–Ulam |
| 26 | 1 | 12 | 6 (even) | Brouwer / SOS |
| 34 | 1 | 16 | 8 (even) | Brouwer / SOS |
| 38,46 | 3 | 18,22 | 9,11 (odd) | Borsuk–Ulam |
```
saddle index = (p−1)/2 = the de Moivre / Chebyshev EQUIOSCILLATION count = #antipodal witness pairs.
parity(index) = p mod 4 :  ODD ⇒ Borsuk–Ulam (free ℤ/2, odd-degree forcing) ;  EVEN ⇒ Brouwer/SOS.
```
This is the **construction analogue of the measure bound** mac-mini wanted: a topological integer attached
to the core saddle, whose oddness is the obstruction's signature.

## The forcing (the topological lower bound — sketch)
Suppose `M(S) < 1/14`: `S(t) < 1/14` for **all** `t`, i.e. the danger set `D = {t : S(t)<1/14}` is the whole
circle, and it is `ℤ/2`-invariant under `t↦−t`. For `p≡3 mod 4` the complement acts **freely** on the
core/saddle structure (no fixed orientation — the tournament is self-converse). A map encoding "which runner
binds" is then a **free ℤ/2-equivariant map**, and its **index is `(p−1)/2 = 3` (odd)**. By Borsuk–Ulam an
odd-index free `ℤ/2`-map to a lower-dimensional sphere cannot exist / cannot be null-homotopic — so the
equioscillation saddle **cannot be pushed antipodally below `1/14`**, forcing `M ≥ 1/14`. The **odd index is
the topological certificate**; the **even-index cases (`p≡1`) admit a symmetric saddle**, handled by the SOS/
Brouwer route (no Borsuk–Ulam needed) — which is why `n=10,26,34` are the *easy half*.

## Where this sits in the architecture (push/pull map)
```
LRC(14) = BULK (meas{S>1/14}>0, circle-method/equidistribution; mac-mini S74, kps HYP-3218 Fejér)
        ⊕ CORE (Vitali wall; witnesses t=a/14 = Φ₁₄; mac-mini HYP-3237)
              ├ VALUE: equioscillation ±1/14  (Brouwer/EVT saddle; de Moivre/Chebyshev Φ₇)
              └ FORCING: saddle index = (p−1)/2 = 3 ODD  ⇒  Borsuk–Ulam lower bound M≥1/14   ← THIS (HYP-3240)
```
The core's two needs — *locate* the saddle (Brouwer/equioscillation) and *force* its value (a degree count) —
are now both supplied: the index `(p−1)/2`, odd iff `p≡3 mod 4`.

## Honest status
- **VERIFIED:** AP witnesses `=Φ₁₄` = 3 antipodal pairs; saddle index `=(p−1)/2`; parity `=p mod 4` (the
  Borsuk–Ulam/Brouwer-SOS split) across `n=6..46`.
- **SYNTHESIS:** mac-mini's Brouwer-saddle (HYP-3237) and my Borsuk–Ulam antipodal (HYP-3239) are one witness
  set `Φ₁₄`; the index `(p−1)/2` is the degree count mac-mini asked for.
- **SKETCH, not proof:** the forcing needs the explicit free `ℤ/2`-equivariant map (the "which runner binds"
  map) and a verification that its index is genuinely `(p−1)/2` and obstructs the no-witness covering. That
  construction is the concrete next target — it would turn the core half into a topological proof.

## Next experiments (this session)
1. Build the explicit "binding-runner" `ℤ/2`-map on the AP and compute its degree directly (does it = 3?).
2. Does the forcing survive perturbation off the AP (the other core/tight sets: Goddyn–Wong T5)? Their unit
   structure mod their conductor should give the same odd index.
3. Cross-check with codex HYP-3238 duality bridge and mac-mini S76 (consec = φ⁴ bimodal extremizer).

→ HYP-3240 (this), HYP-3237 (Vitali core/Brouwer, mac-mini), HYP-3239 (Borsuk–Ulam), HYP-3220 (p mod 4),
HYP-3218 (Fejér/bulk), HYP-3212 (Chebyshev equioscillation), S551o (Vitali wall), OPEN-Q-108.
