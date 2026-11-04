# Independent audit of the exact three-carry classification

2026-11-04. PASS. Audited `../other/PROOF.md`, `verify.py`, and `verification.json` as hashed in this packet's receipt.

## Symbolic review

All six carry recurrences follow correctly from the physical low-bit recursion. In particular the two less immediate odd cases are

- D1(2n+1) = 3[j(n)+n+p] − j(3n+2) = D2(n)+3n+3p;
- D2(2n+1) = 3[j(n)+n+p] − [j(3n+2)+3n+2+p] = D2(n)+2p−2.

The three zero base cases are correct. Every stated interval inclusion in the simultaneous strong induction is valid, including n=0 and the use of p<=n in the nonzero even D1 case. Every recursive argument is smaller at positive input.

Substitution at 4l+1 and 4l+3 gives the displayed formulas exactly. The former is strictly negative. The latter exceeds6 for l>=5, and the five exact residual values are2,6,8,12,8. Therefore D0=4 is impossible and D0=6 occurs exactly at odd part7.

The clock coefficient comparison correctly yields 2^k g=3f and k=a−b for source and target two-adic orders a,b. Arbitrary positive and negative fixed integer shifts are covered; the stated lower bound on t keeps both prefixes nonempty. The converse strongness is justified by the actual terminal suffix ports, not an abstract signature.

The claim that no two-step word over {2,4} stays in the chart is valid: the sole allowed first step has odd source7 and odd target21, and the target cannot take another letter. The nontransitive +2 entry, prescribed-single-clock returns, and native legality are explicitly excluded. No global spectrum or two-branch closure conclusion is asserted.

## Independent controls

The audit uses literal ordered subset covers (not either author's binary recurrence) to compute J on1,539 masks, then checks all six carries for256 inputs, simultaneous bounds on512 inputs, and digit classification on511 positive inputs. The author checker was rerun normally and under optimization: outputs match exactly and agree semantically with the published receipt. These controls audit the proof; the infinite conclusion is supplied by the symbolic induction above.

No defect or phase is reset. The physical joint pair formula is the correct named cover count; this obstruction does not construct a new child and consequently claims no phase transport or native guard witness.
