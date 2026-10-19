# A minimum-endpoint interval cycle through 482,3362,6242

**2026-10-19. Exact finite certificate.**

## The result

There is a cyclic ordering of every integer in `[1,6242]`, each exactly once, beginning

    482,3362,6242,4783,3498,6106,4919,4685,3596,6008,
    5017,3264,2665,4391,5213,5812,2469,4587,1342,5714,

such that every adjacent sum, including the closing sum, is a square.

The prescribed AP seed is retained at its original labels and as a consecutive path: 482+3362=62² and 3362+6242=98². The extra graph edge 482–6242 has root 82 but is not traversed by the larger cycle. The endpoint 6242 is the smallest possible endpoint of an interval `[1,N]` containing this exact seed, simply because one seed vertex is 6242. No claim of a difficult lower-bound theorem is attached to that observation.

The canonical one-integer-per-line cycle has 30103 bytes and SHA-256

    a892ae08c5304c35fae319dfcd9429faa5ea37b2f21e89d3a715e17912ceca6d.

Gerbicz's 2018 construction already proves ordinary square-sum interval Hamiltonicity for all n>=32, as recorded in his [OEIS comment](https://oeis.org/A090461) and [original discussion](https://www.mersenneforum.org/node/17331/page2). The extra property certified here is the designated consecutive seed. No claim that this finite prescribed-seed instance was previously unknown is made, and no all-large-N prescribed-seed theorem follows from this one witness.

## 1. A compact input and a controlled split lift

`single_cut_input_127.json` gives an explicitly checked square-sum path through `[1,127]`. It starts 69,127, ends 10, and contains 7,2 at one-based positions 102,103.

For each integer c in [−24,24], lift its j-th entry a_j to

    49a_j+(−1)^j c.

Internal edge sums are multiplied by 49 and remain squares. For a fixed a_j, the 49 choices of c give exactly `[49a_j−24,49a_j+24]`. Thus all lifted blocks together partition `[25,6247]`, irrespective of the position parity of a_j.

Cut the lifted 7–2 edge in each block. Its two segments have endpoint pairs

    (3381−c,343+c),       (98−c,490−c).

The final minus sign uses the odd input length 127. Cutting changes no vertex coverage. The compact word `final_split_p1_d24_cutodd_reverse.json` uses each of the 98 segments and each singleton 1,...,24 once, in specified orientations. Every connector sum is square. The word starts at 1 and ends at 8, so its closure is square as well.

This yields a cycle through `[1,6247]`. It contains the designated seed because the second segment for c=8 ends at 482 and joins the first segment for c=19, which begins 3362,6242:

    490−8=482,
    3381−19=3362,
    49·127+19=6242.

This is a one-input lift. No unused companion path or unproved nice-pair existence is needed for this finite construction.

The controlled split has a structural purpose. In an intact-block lift with input endpoints 10 and 69, the 27 ports in `[3379,3405]` have no square-sum connector among the old ports and singletons. The new quotient-2 ports repair this obstruction via z+(59²−z)=59², with 59²−z in `[76,102]`. The listed glue word verifies global assembly, rather than assuming that the local repairs automatically join into one cycle.

## 2. Removing the five excess labels

Rotate and orient the 6247-cycle to start 482,3362,6242, then delete 6243,...,6247. This leaves exactly five nonsquare cyclic adjacencies. The file `seeded_cycle_6242_repair_trace.json` supplies 24 exact two-edge reversals which repair them.

A move cuts cyclic edges a–b and c–d and reverses the intervening segment, replacing those edges by a–c and b–d. The vertex set does not change, and all other undirected adjacencies remain unchanged. The trace lists the two cut indices, both removed and added edge pairs, and the number of nonsquare adjacencies before and after each move. That number never increases and ends at zero. Neither protected edge 482–3362 nor 3362–6242 is ever cut.

The local search used to discover these 24 moves is not part of the proof. The independent checker replays the fixed trace, checks every move and the exact final cycle, and requires no random choices or optimizer.

## 3. Reproduction

Keep the three JSON inputs and the verifier in the same directory. Run:

    python verify_seeded_interval_cycle_6242_20261019.py
    python -O verify_seeded_interval_cycle_6242_20261019.py

Both produce byte-identical receipts. The verifier uses only the Python standard library, with explicit exceptions instead of removable assertions. Its 12998 checks include the input path, every lifted and final cyclic edge, exact interval coverage, seed preservation, all 24 repair moves, and both canonical cycle hashes.

To regenerate the optional compressed final witness, run:

    python verify_seeded_interval_cycle_6242_20261019.py --write-cycle

`seeded_interval_cycle_6242.txt.gz` decompresses to the canonical one-integer-per-line ordering. The plaintext hash above identifies the mathematical witness independently of gzip header differences.

This packet establishes the finite minimum-endpoint example. General prescribed-seed extension, the missing residue/parity templates of a uniform induction, and common-center AP-petal compatibility remain distinct questions.
