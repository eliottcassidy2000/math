# Independent audit: singleton–triangle return termination

2026-11-04. **PASS within the stated grammar and clock-uniform notion.** No mathematical correction requested.

Audited `../spectrum/PROOF.md`, `verify.py`, `candidates.json`, and `receipt.json`. Replayed producer normal and optimized executions: outputs agree and reproduce all reported counts.

## Mathematical review

1. The six component symbols exhaust marked singleton/C3 strong components up to marked isomorphism. A nonempty transitive prefix and the stated root interface are strongly connected exactly when the terminal component has a port. The source attachment architecture and signature definitions are preserved.
2. The F/h recurrence yields nonoverlapping one- and two-bit chunks. Positive F, exact triangle count, and a nonzero final chunk bound total width by binary length plus one. The only padding possibility is the high zero of a final C1. Thus the inverse truly covers every length, not merely a bounded graph census. Exact chunk boundaries recover the word uniquely.
3. For nonzero-modulo-three target constants, an initial zero-port triangle is impossible since it multiplies c by three. Removing initial S0 blocks preserves c and divides F by powers of two. A first marked symbol gives v2(F) in {0,1}. The two normalized inverse lists therefore cover all fixed integer shifts, including negative shifts at sufficiently large legal clocks. No single-clock conclusion is smuggled into the uniform statement.
4. An independent enumerator chooses disjoint two-bit intervals rather than using the producer's recursive digit walk. It reproduces normalized candidate counts 5,10,18,29,50 for odd F=21,63,189,567,1701. A closed signature formula, summing same-block, adjacent-block, and opposite-order joint contributions, independently reproduces the full candidate signatures. Every one of the 112 candidates and all 280 digit comparisons agree. The resulting language is exactly {empty,2,4,22,42,422,4222}. Therefore neither original sibling has a +4 child, and no unsearched deeper state can occur.
5. The 3-adic proof correctly distinguishes same-block, adjacent-block, and opposite-order joint terms. Only adjacent triangle blocks lose two factors of three; their contribution is a_i*a_(i+1). Subtracting J yields the negative sign. The k=1 correction includes the same-triangle a(a−1) term and both adjacent singleton bits. The k=2 coprimality iff and k>=3 obstruction follow. The old coarse k−2 bound is explicitly credited and is not repackaged as new.
6. The full seven-state table, displayed H values at t=1, terminal strongness, and accepted zero-shift representatives agree with the independent reconstruction and producer literal graph DP.

## Independent executable controls

Run `python3 audit.py` and `python3 -O audit.py` in this directory. Their JSON outputs byte-match.

- 112 distinct normalized candidate words and all signature coordinates.
- 280 successor comparisons and exact seven-state language.
- 1,792 zero-prefix normalization/triangle-obstruction checks.
- 624 leading-coefficient checks with mixed singleton/triangle words, varying marks, two through fourteen triangles, and clocks 1,2,9.

The producer's independent literal subset DP additionally checks 882 small component signatures, 492 inverse roundtrips, 72 full TT1 attachments, all 112 actual candidates, and literal attachments at two clocks for each accepted state. These finite controls corroborate, but do not replace, the inverse and normalization completeness arguments.

## Scope boundary

This proves termination of one explicit singleton/C3 ordered-component return grammar under fixed clock-uniform ternary digit equations. It is not a global Hamilton-path count spectrum obstruction, does not cover larger strong components, does not forbid individual-clock coincidences, and does not establish old-name-preserving/native-operation legality. The author's final scope paragraph respects these boundaries. No global conjecture closure is endorsed by this audit.
