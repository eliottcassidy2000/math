# Width eight is necessary for the staggered criterion

**FINITE-EXACT, independently audited. F(7,137) is false; F(8,137) remains OPEN.**

At base X=153 and residue r=1, every shift j=0,...,7 is impossible at every original source order. Shift j=8 is realized by the explicit rooted tournament in SHARP_NORMAL.json, with (h,s,L)=(77,25,85) and (B,s)=(137,25). Thus any true F(J,137) needs J>=8.

This is a counterexample to a narrower sufficient profile condition, not a gap in the Hamilton-path spectrum. The repaired physical construction has exactly 461 Hamilton paths. No eventual threshold is proved.

Read PROOF.md for the complete source-order reduction and enumeration. The primary bank contains one representative of each of the 1,113 isomorphism classes of strong tournaments with H<=153 through order 10, and finds 535 strong rooted profiles. A separate rooted generation agrees. Independent closure checks cover all 281,628 mixed extensions, with direct Hamilton-word counts and exact isomorphism backtracking.

The audit also independently certifies the complete 480-pair box B<=137, s<=26; its proof, all witnesses, and replay files are preserved under audit/reviewed-complete/. This is a separate finite-box conclusion, not an infinite tail.

Reproduce the focused result:

    g++ -O3 -std=c++17 enumerate.cpp -o enumerate
    cd run153
    ../enumerate 153 26 > replay.json
    cmp replay.json NORMAL.json
    cd ..
    python3 check_sharp.py > replay.json
    cmp replay.json SHARP_NORMAL.json

The independent audit is in audit/REVIEW.md, with its own exact manifest and replay instructions. Every reviewed entry is retained. Published dependency provenance is in dependencies/PROVENANCE.json.

The proof promotion changes only its status line. The separately logged CODE_GUARD.diff prevents unsupported 64-bit canonical encoding above order 11 and rejects nonpositive or excessive endpoint caps. The declared 153/26 and 137/26 runs reproduce unchanged bytes after this guard. This release certifies only the declared finite runs; the unguarded frozen originals are retained for audit, not advertised as arbitrary-bound tools.
