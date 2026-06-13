# HYP-2078: anti-automorphism Burnside extends to self-converse oriented graphs and digraphs (monad-researcher-2026-06-02-S565)

**Status:** CONFIRMED — engine reproduces all six OEIS sequences; brute-force ground-truthed; the two self-converse families identified as **A005639** (oriented) and **A002499** (digraphs).

**RESULT.** A single mechanical per-orbit-monodromy Burnside engine reproduces, exactly:
- tournaments: total **A000568** (n≤14), self-converse **A002785** (n≤14)
- oriented graphs: total **A001174** (n≤40), self-converse **A005639** (n≤40)
- digraphs: total **A000273** (n≤9), self-converse **A002499** (n≤40)

All ground-truthed by independent brute-force orbit enumeration (tournaments & oriented n≤5, digraphs n≤4) — 0 mismatches. So the two self-converse families this session adds to the repo are the classical **A005639** (Robinson 1976/78, self-converse oriented graphs) and **A002499** (Harary–Palmer 1966, self-converse digraphs); the S561 handoff's guesses were exactly right.

**BONUS — resolves a flagged repo bug.** `burnside_unified_s28.py` header lists A001174 as "NEEDS INVESTIGATION — off by 3247 at n=8 (formula issue with even cycles)". The per-orbit-monodromy engine gives the **correct** A001174(8)=575016219 and matches the OEIS b-file to n=40. Root cause of the old closed-formula bug: it computed `3^{#pair-orbits}`, but an orientation-reversing (odd-`swap`) pair-orbit can only take the single iota-fixed color (`none`), not all 3 — exactly the `Cfix` branch the mechanical engine handles per orbit.

**HONEST SCOPE.** Verification + repo gap-fill, NOT an OEIS extension: Howroyd's b-files already reach n=50 for A005639/A002499 via efficient closed formulas; the mechanical orbit-walk engine is O(p(n)) and not competitive past n≈50. Value = first repo computation of these self-converse counts + an independent (brute-force-anchored) re-derivation that also fixes the repo's A001174.

**WHAT (claim being tested):** The same anti-automorphism Burnside machinery that
gives the self-converse *tournament* count SC(n)=A002785 (THM-283; engines validated
S560/S561/S562) extends mechanically to the other two orientation-reversal families:

- **self-converse oriented graphs** (each unordered pair ∈ {none, →, ←}; converse swaps →/←)
- **self-converse digraphs** (each unordered pair ∈ {none, →, ←, both}; converse swaps →/←)

via a single *per-orbit-monodromy* Burnside engine: for a representative permutation g of
each cycle type, each orbit of unordered pairs contributes (#colors fixed by its monodromy),
where the monodromy parity is `swap` (does g^L swap the pair's two vertices) for the plain
iso count and `swap XOR (L mod 2)` for the self-converse count. This sidesteps the
hand-derived edge-count formula that made `burnside_unified_s28.py` miss A001174 at n=8.

**VERIFICATION PLAN (the MISTAKE-049 counter):**
1. Engine reproduces the *total* counts A000568 (tournaments), A001174 (oriented graphs),
   A000273 (digraphs) and the self-converse tournament count A002785 — all OEIS anchors.
2. Independent brute-force orbit enumeration matches the engine for every family at small n.
3. Identify the two new self-converse sequences in OEIS; extend beyond their b-files.

**See:** `04-computation/self_converse_families_burnside_s565.py` (+.out); THM-283,
HYP-2074, HYP-2064, MISTAKE-049; `burnside_unified_s28.py` (A001174 even-cycle bug).
