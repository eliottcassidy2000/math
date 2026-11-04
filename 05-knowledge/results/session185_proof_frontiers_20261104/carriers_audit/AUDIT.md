# Independent audit: growing five-vertex carriers

**PASS, no mathematical corrections required.** Audited `../carriers/PROOF.md`, `verify.py`, `inherited.py`, `bank.json`, and the six physical escape witnesses. Source hashes are recorded in `normal.json`; normal and `python -O` receipts match, as do the author's independently rerun normal/optimized receipts.

## Independent checks

`check.py` does not import producer helpers. It uses forward-extension endpoint DP, reverse adjacency for start counts, and explicit subset convolutions for all seven signature coordinates. It reconstructs every labelled tournament through order five and canonicalizes every strong graph. The labelled strong counts are 1,0,2,24,544 for orders one through five. The resulting nine types exactly match the supplied bank. All 218 marked/unmarked bank rows have independently recomputed signatures and no duplicate type/mask keys.

All 867 candidates in the constant-bounded fixed-front inverse are recomputed from literal full suffix rows and port masks, rather than merely applying the producer's composition recurrence. Their full root attachments are checked strong. No candidate has constant −35 or −29 and odd coefficient 63. This verifies the actual finite certificate on which the all-length statement rests.

The positive normalized example has the stated strong five-vertex front and signature (11,100,12,27,0,1,4). Its first physical suffix, with tail mask645, has F64512 and c−21849. Its seventeen-vertex attachment at clock1 is strong and has H107175. Integrality and valuation for every family index follow from 2^20≡1 modulo25 and the exact algebra in the proof, rather than extrapolation from twelve examples.

All six SCC7/8 witness signatures, strongness assertions, and literal full attachment counts at clocks1,2,3 pass independent recomputation. This audit checks those physical numbers, not the companion packet's complete address tree.

## Proof review

1. The tail scan uses low-to-high port bits in the ordered transitive tail, with F=u and initial-start port bit u mod2. The joint term e(u mod2)+Gu is correctly retained in the ordered sum.
2. A nonempty marked front has G>0: a singleton path starting at any port, together with a Hamilton path on its complement, contributes to the defining two-slot count. Every remaining subtracted term is nonnegative, so u≤floor((h−J−C)/G). The empty tail is inspected separately. Canonical positive u represents every transitive tail that meets the terminal-port strongness guard, including initial unmarked tail vertices.
3. Equality of positive odd parts is necessary and sufficient for the coefficient equation at some fixed integer dyadic clock shift. Clock positivity restricts usable source clocks, as stated. Thus the inverse covers every tail length and shift in its stated fixed-front class.
4. The SCC5 family starts with a marked strong component; its unbounded v2(F) is not leading unmarked-singleton padding. Its constant tends to negative infinity, so the example does not realize the missing fixed constants.
5. For the auxiliary fixed-F bound, the last marked SCC has F≥1 and the preceding prefix satisfies p≥2^|P|, giving the stated total-order bound. This is correctly distinguished from an all-clock bound.

## Scope retained

The negative conclusion covers one **marked strong front of order at most five**, followed by an arbitrary transitive tail. It does not cover unmarked fronts, arbitrary SCC≤5 component words, all high-valuation carriers, native-operation guards, infinite branch closure, or global spectrum completeness. The packet states these limitations accurately. No forbidden central ledger was read and no software was installed.
