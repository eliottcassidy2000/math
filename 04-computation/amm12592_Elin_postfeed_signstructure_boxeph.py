"""LANE D2 (E-lin) -- exact instrumentation of the POST-FEED autonomous
phase of the T6 clamped-Pascal junk flow at linear slack.

Motivation: the proved skeleton (T4b + T6a/T6b + exact i_feed) shows that at
D0 >= ceil(eps* R), eps* = 2(1-g-g^2)/(3+2g) ~ 0.02118, no death can occur
during the feed phase.  The ONLY remaining obligation for E-lin is the
autonomous (feed-free) flow from the feed-end state: no death + capture by
row R-2.  This script measures, exactly, the structure any envelope proof
must exploit:

  per sampled row: junk support [tmin, tmax], L1 bits, front gap d - tmax,
  ALTERNATION fraction (fraction of adjacent nonzero-cell pairs with
  opposite signs, distance-1 pairs only), number of sign blocks,
  cap-margin at the front (bits of cap(front) - bits of junk(front)),
  and the L1 decay rate over the post-feed phase.

Everything exact; floats only in *reported ratios* (never decisions).
"""
import os, sys, json, importlib.util, contextlib, io
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
spec = importlib.util.spec_from_file_location(
    "fast", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fast = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fast)
floor_gamma_star = fast.floor_gamma_star
two_G_coeffs = fast.two_G_coeffs
initial_junk = fast.initial_junk


def run_instrumented(R, D0, sample_every):
    g = two_G_coeffs(R)
    d = floor_gamma_star(R) + D0
    j, junkL1, c0 = initial_junk(R, d)
    rows = []
    i_feed_obs = None

    def snap(i, d, j, junkL1, fed):
        ts = sorted(j)
        if not ts:
            return {"i": i, "d": d, "empty": True, "fed": fed}
        signs = [(t, 1 if j[t] > 0 else -1) for t in ts]
        adj = [(a, b) for (a, sa), (b, sb) in zip(signs, signs[1:])
               if b == a + 1]
        alt = sum(1 for (a, b) in adj
                  if j[a] * j[b] < 0)
        blocks = 1 + sum(1 for (ta, sa), (tb, sb) in zip(signs, signs[1:])
                         if sa == sb and tb == ta + 1)
        # sign blocks: count maximal runs of same sign among adjacent cells
        nsb = 1
        for (ta, sa), (tb, sb) in zip(signs, signs[1:]):
            if sa != sb:
                nsb += 1
        front = ts[-1]
        capf = 2 * comb(d - 1, front)
        return {"i": i, "d": d, "tmin": ts[0], "tmax": front,
                "gap": d - front, "nnz": len(ts),
                "junkL1_bits": junkL1.bit_length(),
                "adj_pairs": len(adj), "alt_adj_pairs": alt,
                "alt_frac": (round(alt / len(adj), 4) if adj else None),
                "n_sign_runs": nsb,
                "front_junk_bits": abs(j[front]).bit_length(),
                "front_cap_bits": capf.bit_length(),
                "fed": fed}

    rows.append(snap(0, d, j, junkL1, True))
    outcome = None
    for i in range(1, R - 1):
        dn = floor_gamma_star(R + i) + D0
        delta = dn - d
        w = {}
        if delta == 0:
            for t, v in j.items():
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + v
        else:
            for t, v in j.items():
                w[t] = w.get(t, 0) + v
                w[t + 1] = w.get(t + 1, 0) + 2 * v
                w[t + 2] = w.get(t + 2, 0) + v
        fed = False
        if dn + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn + i]; fed = True
        if delta == 1 and dn - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[dn - 1 + i]
            w[1] = w.get(1, 0) + g[dn - 1 + i]; fed = True
        if not fed and i_feed_obs is None:
            i_feed_obs = i - 1
        jn = {}
        junkL1 = 0
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            P = comb(dn - 1, ta)
            Pprev = comb(dn - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    if u != v:
                        jn[t] = v - u
                        junkL1 += abs(v - u)
                Pprev = P
                P = P * (dn - 1 - t) // (t + 1) if t < dn - 1 else 0
        if dn in jn:
            outcome = {"status": "DIE", "row": i}
            rows.append(snap(i, dn, jn, junkL1, fed))
            break
        j, d = jn, dn
        if (i % sample_every == 0) or (i_feed_obs is not None
                                       and abs(i - i_feed_obs) <= 2):
            rows.append(snap(i, d, j, junkL1, fed))
        if not j and not fed and d + i > R - 1:
            outcome = {"status": "CLOSED", "capture_row": i}
            break
    return {"R": R, "D0": D0, "i_feed_obs": i_feed_obs, "outcome": outcome,
            "rows": rows}


if __name__ == "__main__":
    jobs = [(2048, 128, 32), (2048, 64, 32), (2048, 38, 32), (2048, 37, 32),
            (4096, 256, 64)]
    out = {}
    path = os.path.join(RESULTS,
                        "amm12592_Elin_postfeed_signstructure_boxeph.json")
    for R, D0, se in jobs:
        r = run_instrumented(R, D0, se)
        key = f"R{R}_D0{D0}"
        out[key] = r
        json.dump(out, open(path, "w"), indent=1)
        # console digest: post-feed L1 trajectory
        pf = [row for row in r["rows"] if not row.get("empty")
              and not row["fed"]]
        digest = [(row["i"], row["junkL1_bits"], row["gap"],
                   row["alt_frac"]) for row in pf[:40]]
        print(key, r["outcome"], "i_feed=", r["i_feed_obs"], flush=True)
        print("  (i, L1bits, gap, altfrac):", digest, flush=True)
    print("DONE ->", path, flush=True)
