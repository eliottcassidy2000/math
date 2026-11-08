# Set-root graft constructions and the finite-menu boundary

**PROVED, independently audited PASS. Chan–Pak remains OPEN.**

- A positive cover formula extends single-root grafts to a set R of roots, with sharp maximum degree min(|R|,|A-R|)+1.
- For interfaces with min(|R|,|A-R|)<=r, exactly r+1 rational linear gadget-summary coordinates are necessary and sufficient; physical strong-compatible contexts witness minimality.
- The output is strong exactly when A's first SCC meets R and its last SCC meets A-R. Every such graft satisfies the sharp bound H>=2^|B|+1.
- A fixed finite source bank and a k-type ordered gadget menu give only O((log X)^k) counts below X. Arbitrary outer ordered sums still give at most exp(O((log X)^(k/(k+1)))) distinct counts. This is a scoped language boundary, not an obstruction to general tournaments or recursively adaptive grafts.
- Explicit strong families include -5*3^m+4*12^m+2*33^m, and binom(k+l,k) pairwise nonisomorphic strong tournaments sharing the complete rooted graft state and count 3*3^k5^l+12*12^k32^l.

The general cover calculus and ordered-sum transform are inherited from THM-3121 and THM-3166; the single-root state transfer and context equivalence are inherited from Sessions225–226. The proof distinguishes these dependencies from its set-root specialization and consequences.

Files:

- `PROOF.md`: all-order proofs and scope boundaries.
- `verify.py`: standalone standard-library exact checker; all guards remain active under optimization.
- `CHECK.json`, `CHECK.optimized.json`: matching generated receipts with literal witnesses.

Reproduce:

    python3 -I -B verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.json

The source/gadget universe is all labelled tournaments of orders 1 through 4, all source root subsets, excluding only 4-by-4 pairs: 16,814 literal pairs, 7,934 strong. There are no inherited connectivity filters. Further controls cover 121 ordered products, 16 sharp-interface cases, three explicit three-mode graphs, physical minimal-summary rank matrices through r=4, and eight full-state collision-family witnesses through order15.

An initial draft incorrectly asserted equality for an unequal-shore leading coefficient. The literal verifier rejected it before publication. The correct statement is a positive lower bound for unequal shores and exact equality for equal shores, because unequal shores can support maximal-run paths with only one occupied boundary slot. No other theorem uses that rejected equality.
