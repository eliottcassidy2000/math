#!/usr/bin/env python3
# Auditor "hypotheses": independent re-derivation of the remainder of the product law at 12 from the saved
# histograms (different code path from even_moduli_census.py), sheets 0 and 1, X = 2^22 .. 2^28.
import numpy as np, math, json
divs = [d for d in range(1, 316) if 315 % d == 0]
for E in (22, 24, 26, 28):
    hist = np.array(json.load(open("even_moduli_hist_%d.json" % E))["hist"], dtype=np.int64)
    out = []
    for sh, cs in ((0, -1), (1, 0)):
        T = hist[sh].sum()
        D3 = sum(hist[sh, i].sum() for i, d in enumerate(divs) if d % 3 == 0)   # smaller state of type 3 (w.r.t. m = 3)
        D1 = T - D3
        N = np.arange(2520)
        cnt12 = np.array([hist[sh][:, N % 12 == c].sum() for c in range(12)], dtype=float)
        worst = 0.0
        for c in range(1, 12, 2):
            pred = 0.5 * (D1 / 3 + (D3 if (c - cs) % 3 == 0 else 0))
            worst = max(worst, abs(cnt12[c] - pred))
        # the remainder at 3 alone
        cnt3 = np.array([hist[sh][:, N % 3 == c].sum() for c in range(3)], dtype=float)
        w3 = max(abs(cnt3[c] - (D1 / 3 + (D3 if (c - cs) % 3 == 0 else 0))) for c in range(3))
        out.append((12 * worst / T, 3 * w3 / T, 9 * D3 / T))
    print("X = 2^%d: remainder at 12 (x12): %.4f / %.4f ; remainder at 3 (x3): %.4f / %.4f ; F_3 * 9: %.3f / %.3f" % (
        E, out[0][0], out[1][0], out[0][1], out[1][1], out[0][2], out[1][2]))
