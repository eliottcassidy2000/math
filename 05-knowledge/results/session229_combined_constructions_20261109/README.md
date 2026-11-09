# Session229: combined constructions toward Chan–Pak

**COMPLETE research session; independently audited. Chan–Pak completeness remains OPEN.**

Target: actual tournament realization of every positive odd Hamilton-path count except 7 and 21. This session combines older binary core constructors with newer sidecar fibers, graft inverses, and recursive path-cover compilation. It proves scoped constructive results and precise barriers, not the all-odd theorem or a new unrestricted coverage cutoff.

## Strongest constructive results

1. **Universal adjacent-count repairs.** Preserving the complete endpoint interface and changing only its T coordinate yields consecutive odd outputs under every common legal binary prefix/suffix. A five-vertex core supplies a pair; two six-vertex core topologies supply a triple. In particular, for every b≥0 there are literal strong tournaments with counts 258·2^b+21,258·2^b+23,258·2^b+25. The maximal fixed-interface run lengths through core orders 4, 5, 6 are exactly 1, 2, 3 in the complete labeled census. [Proof](dynamic/PROOF.md), [audit](dynamic_audit/AUDIT.md).

2. **Unsaturated cores with unit switches.** The H95/Q1276 and H389/Q5600 cores have gcd one and full-interface gap-two switches. The first supplies strong adjacent-count pairs at every order≥9. Unsaturation belongs to the core and need not survive star closure: the H681 output has gcd(H,Q)=3. A rigid H29/Q320 block separately gives four equal-alpha transfers and 4^k distinct full transfers at a common one-sided coordinate. Transfer distinctness is not automatically scalar distinctness. [Proof](rich/PROOF.md), [audit](rich_audit/REVIEW.md).

3. **Larger exact recursive language.** Arbitrary three-block C3 synthesis retains the sharp node bound H≥2^(n−2)+1. Its complete grammar through 4095 has 620 counts, compared with 577 previously; adding the inherited H19 seed gives 660. The explicit operation C3[TT2,TT2,A] escapes the old grammar at 137. A structural argument leaves 29 absent at every recursive depth in the enlarged language. A general finite cycle-certified template compiler is also proved. [Proof](recursive_complete/PROOF.md), [independent complete compiler](ternary_audit/REPORT.md).

## Exact calculus and remaining barriers

- **Nearest-neighbor defects and bounded carries.** Every nonempty intervening block kills long-range equal-alpha interactions; only adjacent quadratic terms remain. Equal-U words have an exact finite carry alphabet, with alpha-word rigidity when gcd(H,Q)=1. Literal strong H9/Q84 carry witnesses have counts 13149 and 16719 at the same sidecar. All masks of the fixed H19 strong-flank template still give only O(X^(3/4)) values≤X; mixing the named H19/H9 blocks remains O(X^(15/16)). [Proof](overlap/PROOF.md), [audit](overlap_audit/REPORT.md).
- **Finite all-gadget physical descent.** Every proper graft has explicit target bounds and strictly smaller source and gadget Hamilton counts. The inherited physical-state oracle yields a terminating scalar inverse over all gadgets, including singular ones. Prime targets normalize to reachable-root sources. This proves a finite decreasing search, not that every admissible target has a successful branch. [Proof](inverse/PROOF.md), [audit](inverse_audit/REPORT.md).
- **Why simple amplification stops.** Repeating central unit switches scales their spacing by outside Hamilton counts. Two exact one-vertex extension censuses retain their gcd-one filters and show that the named seeds do not yield an unsaturated universal triple that way. The negative statements apply only to those explicit universes. [Extension proof](rich/EXTENSION_PROOF.md), [independent audit](extension_audit/REPORT.md).

[Final synthesis](SYNTHESIS.md) explains the connections and remaining existence obligation. The next needed theorem is an overlap/coverage argument for the repair windows, a growing T interval that avoids product spacing, or an all-target physically admissible inverse selection. None is inferred from gcd one, large fibers, finite reachability, or scalar descent.

## Verification and final status

`python3 -I -B verify_all.py` replays 36 canonical checks: isolated normal/optimized Python and both C++ census implementations at -O0/-O3. All receipts match. `sha256sum -c SHA256SUMS` binds the final packet. [Verification scope](VERIFICATION.md), [replay receipt](REPLAY.json), [final scope review](FINAL_REVIEW.md), and [publication record](PUBLICATION.md) are included.

The historical candidate/pending headings in inverse/PROOF.md and recursive_complete/PROOF.md are superseded by their linked independent PASS audits. The first-stage recursive packet remains frozen separately. MISTAKE-736 records the saturation and convention repairs. The concurrent first-realizer reservation was preserved, and a blank-line-only navigation adjustment restored the repository documentation budget without changing its mathematics.

The session has ended. Unrestricted Chan–Pak completeness remains OPEN.
