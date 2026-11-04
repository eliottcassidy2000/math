#!/usr/bin/env python3
"""fo_rule_modelcheck3.py -- HYP-9463 (4): uniform quantifier-free rules with at most THREE atoms S(c_i x + d_i),
|c_i|, |d_i| <= BOX (default 4), any Boolean function of the three atom bits (256 truth tables), tested as in
fo_rule_modelcheck.py (symmetric output, well-defined on dilation classes, bijection onto multiplier-even classes)."""
import sys, itertools, time
import numpy as np
from fo_rule_modelcheck import setup

def main():
    args = sys.argv[1:]; BOX = 4
    if args and args[0] == "--box": BOX = int(args[1]); args = args[2:]
    primes = [int(a) for a in args] or [13, 17]
    cs = [c for c in range(-BOX, BOX + 1) if c != 0]; ds = list(range(-BOX, BOX + 1))
    atoms = [(c, d) for c in cs for d in ds]
    tts = [np.array([(i >> j) & 1 for j in range(8)], dtype=np.uint8) for i in range(256)]
    survivors = None
    for p in primes:
        t0 = time.time(); D = setup(p); A = D["skew_arr"]; pow2 = D["pow2"]; neg = [(-x) % p for x in range(p)]
        cls_sym = D["cls_sym"]; skew_cls = D["skew_cls"]; nskew = D["nskew"]; even_cls = D["even_cls"]
        reps = np.unique(skew_cls, return_index=True)[1]
        perms = {(c, d): [((c * x + d) % p) for x in range(p)] for (c, d) in atoms}
        passing = []; tested = 0
        for t in itertools.combinations_with_replacement(atoms, 3):
            a1 = A[:, perms[t[0]]].astype(np.int64); a2 = A[:, perms[t[1]]].astype(np.int64); a3 = A[:, perms[t[2]]].astype(np.int64)
            idx = 4 * a1 + 2 * a2 + a3
            for i in range(256):
                tested += 1
                out = tts[i][idx]
                if out[:, 0].any(): continue
                if not np.array_equal(out, out[:, neg]): continue
                masks = (out.astype(np.int64) * pow2).sum(axis=1)
                cls = np.array([cls_sym[int(m)] for m in masks])
                ok = True
                for k in range(nskew):
                    sel = cls[skew_cls == k]
                    if sel.min() != sel.max(): ok = False; break
                if not ok: continue
                img = set(cls[reps].tolist())
                if len(img) == nskew and img == even_cls: passing.append(t + (i,))
        print("p=%d: t_p=%d e_p=%d rules tested=%d passing=%d (%.0fs)" % (p, nskew, len(even_cls), tested, len(passing), time.time() - t0), flush=True)
        if passing[:6]: print("   examples:", passing[:6], flush=True)
        S = set(passing); survivors = S if survivors is None else survivors & S
        print("   passing at every prime so far:", len(survivors), flush=True)

if __name__ == "__main__":
    main()
