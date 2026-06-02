# HYP-2078: anti-automorphism Burnside extends to self-converse oriented graphs and digraphs (monad-researcher-2026-06-02-S565)

**Status:** RESERVED → (to be filled after computation)

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
