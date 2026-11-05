# Independent audit of the two six-vertex additions

**PASS.** Session 195, 2026-11-05. The frozen source proof has SHA-256

    3d5ee515a2465790dd01bd8402ef2d6952b84feb16de77c5a9b24a3abc25ec3c

`receipt.json` pins all twelve reviewed source artifacts, including the proof, two implementations, all inventories, both target witnesses, and the single-addition and bounded-spectrum certificates.

## Proof audit

The one-seam decomposition correctly gives all four coordinates `(h,e,R,d)` for every physical rooted input, with empty path slots and the singleton correction accounted for. The two displayed matrices are actual six-vertex pieces, not formal numerical profiles. The H67 construction starts at a singleton. The displayed H173 graph uses D and then the second piece; using the isomorphic rooted singleton-gadget predecessor also gives the same state.

The bounded search is complete at arbitrary composition depth. Every non-D operation has output H at least the previous R, while D leaves H unchanged and increases R. A target witness can have terminal D operations removed. Therefore an R>X state cannot be an essential preterminal state in a successful count-at-most-X derivation, and terminal values must still be recorded when their R exceeds X. Merging repeated four-coordinate states is valid because every continuation depends exactly on those coordinates. The finite retained domain and visited-state rule handle neutral loops. No graph-order cutoff is needed.

The minimality claim is correctly restricted to the fixed original 65-piece bank plus K_old,D. All gadgets of orders at most five are already represented. The exhaustive six-vertex inventory includes the deletion coordinate in its signature. Testing each additional signature through target173 is sufficient to decide its ability to reach both 67 and 173 at arbitrary depth. No single six-vertex matrix succeeds, while the displayed pair does. This does not rule out a larger single gadget or another operation architecture.

One wording clarification was requested and applied: maximum gadget order six is necessary and sufficient **to cover 67 and 173, and the stated certified interval**, not a proof of global spectrum coverage. No mathematical correction was needed.

## Independently reproduced computations

The auditor wrote a forward-extension path DP, separately from the source's terminal-removal recurrence, and regenerated the complete inventories from labelled graphs:

- all 1,099 labelled graphs through order 5 give exactly 65 distinct full eight-entry signatures;
- all 32,768 labelled six-vertex graphs and 196,608 marks give exactly 280 distinct full signatures;
- every generated signature set equals the frozen source inventory.

A separately written four-tuple state search, without the source's16-bit packing, gives:

    original bank, bound 173: 1,267 states; missing 7, 21, 67, 173;
    enlarged bank, bound 1000: 21,353 states; missing 7, 21;
    enlarged bank, bound 20000: 2,584,754 states; missing 7, 21.

All280 one-addition searches were rerun and compared candidate-by-candidate with the source certificate. Exactly6 additions reach 67, exactly 5 reach 173, and none reaches both. Every retained-state count agrees as well.

Independent Python graph construction and forward DP check both new gadget signatures, the literal H67/H173 output rows, all four output coordinates and strongness. Literal permutation/cut enumeration independently agrees on both gadgets and both outputs. The transfer is then checked on every labelled input tournament through order 4 and every root, giving 570 physical composition checks. The two stated feedback-vertex deletion controls also pass.

C++ O2 and O0 replays agree byte-for-byte on both inventories, every single-addition result, and both enlarged-bank census totals. Python normal and optimized-mode outputs and receipts also agree after the wording clarification; all checks use explicit exceptions.

## Reproduction

From this audit directory:

    g++ -O2 -std=c++17 -Wall -Wextra -pedantic audit.cpp -o audit
    ./audit 20000 > independent.stdout
    python3 audit.py --source PATH_TO_THE_TWO_GADGET_PACKET
    python3 -O audit.py --source PATH_TO_THE_TWO_GADGET_PACKET

The C++ executable generates `independent-base.txt`, `independent-six.txt`, and `independent-single.txt`. The Python script compares these outputs to the pinned source artifacts and rejects source changes during its run.

The exact enlarged grammar has every required odd value through 20,000. Whether it attains every larger admissible odd integer remains **OPEN**. This audit makes no universal ternary-return or unrestricted-spectrum claim.
