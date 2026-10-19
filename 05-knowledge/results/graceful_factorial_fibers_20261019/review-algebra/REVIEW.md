# Independent algebra review: factorial frozen-tail fibers

19 October 2026. **PASS** for the three stated theorems and their scope. Reviewer read the proof and wrote a separate state-graph checker without reading the author's verifier.

The arithmetic forest construction is valid for every k>=1 and m>=k^2+2k: each replacement endpoint i(k+2) lies outside S by its residue i modulo k+1, and its low endpoint i attaches once to the retained star. These are all non-singleton vertices. Every free-color edge touching a singleton is one of exactly its two stated contacts; singleton-singleton contacts cannot occur because their spacing exceeds k. Contraction therefore forces a perfect matching of colors to singletons. This establishes the full k!2^k count and k! disjoint cubes, rather than only a lower bound.

The width-two classification is exhaustive because two bijections differing at at most two color positions are identical or differ by a transposition. The first-touch lower bound is sound: a color fixed initially and not yet touched remains a fixed point in the CURRENT relative permutation, since no move has changed its value. Its first transposition necessarily merges its singleton cycle with another cycle. Each such decreasing-cycle transposition can first-touch at most two of the tracked colors. Thus 2b>=t really gives q=r+2b>=r+t. Moves that never transpose the other f-t colors can repair at most two of them at a time. No sign-correction shortcut through permutation detours beats the claimed lower bound.

The upper bound sets signs on final transposition participation and handles initially fixed errors separately. The diameter argument covers s=0, s>=2, k=1, and the special k=2 complete graph correctly. No issue found in the relative-permutation convention tau^{-1} pi.

Independent exact BFS checks every target from the identity through k=5 and every target from an additional nonidentity signed-permutation origin at k=5: 8,122 distances total, all equal to the formula. Receipt: `independent_metric_receipt.json`; replay: `python3 independent_metric_check.py` and `python3 -O independent_metric_check.py`.

These checks concern the abstract state graph; the arithmetic realization was verified directly from its formulas and proof. The cited literature attribution was not independently researched in this audit. General graceful-tree existence, released-upper-tail one-color connectivity, and general upper-tail diameter are appropriately left open.
