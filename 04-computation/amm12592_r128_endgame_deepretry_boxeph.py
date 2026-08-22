#!/usr/bin/env python3
"""Deep retry on 'unknown' endgame residuals (relaxation passed, construction
failed) dumped by amm12592_r128_endgame_algebra_boxeph.py hunts.  Retries the
e>=1 sandwich with a much larger DFS budget and a wider scan grid.  Any
success is a genuine S_2 membership witness (verified exactly); prints where
in a hunt it would have closed.  Usage: deepretry.py <unknowns.json> ..."""
import json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_endgame_algebra_boxeph as ea

def main():
    tot = won = 0
    for path in sys.argv[1:]:
        samples = json.load(open(path))
        print(f"{os.path.basename(path)}: {len(samples)} unknown residuals")
        for s in samples:
            tot += 1
            sig, da, db = s["sigma"], s["da"], s["db"]
            w = ea.absorb2_solve(sig, da, db, scan_grid=(4, 10), budget=300000)
            tag = f"deg={len(sig)-1} L1={sum(abs(x) for x in sig)}"
            if w:
                won += 1
                ok = ea.verify_pair(sig, w[0], da, w[1], db)
                print(f"  RESOLVED via {w[2]} ({tag}) verify={ok}")
            else:
                print(f"  still unknown ({tag})")
    print(f"deep retry: {won}/{tot} resolved")

if __name__ == "__main__":
    main()
