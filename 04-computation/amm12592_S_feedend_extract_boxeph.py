"""Lane E2 (Hypothesis S via super-block impossibility) -- feed-end state
extractor.

Runs the certified T6 clamped-Pascal junk flow (rule A, plain clamp) from row 0
through the FIRST AUTONOMOUS ROW i0 (the first row whose load receives no feed,
= i_feed_obs + 1 in the D2 feed-end-state convention), and dumps the exact junk
state there.  This is the input state for the post-feed super-block certificate
`amm12592_S_superblock_certificate_boxeph.py`, which replaces the remaining
~79% of the rows of a full Hypothesis-S verification by a proof.

Everything exact int; no floats in any decision (the engine's float seed in
floor_gamma_star is corrected by exact Lucas/Fibonacci comparisons).
Transport/feed/clamp logic is a line-for-line port of the certified fast engine
`amm12592_transient_fast_junkflow_boxeph.py` (T2+T6 conjugacy); regression
against the 16 published feed-end records (m, j0, d, first_autonomous_row,
junkL1 bits, N/C flags, floor(log2 A_t) profile) is done by the certificate
script.

Usage:  python3 amm12592_S_feedend_extract_boxeph.py R D0 [outpath]
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

spec = importlib.util.spec_from_file_location(
    "fastflow", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fastflow = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fastflow)
two_G_coeffs = fastflow.two_G_coeffs
initial_junk = fastflow.initial_junk
floor_gamma_star = fastflow.floor_gamma_star


def extract_feedend(R, D0, outpath=None, flush_every=256):
    """Run rows 0 .. i0 (first autonomous row) exactly; return state dict."""
    t0 = time.time()
    g = two_G_coeffs(R)
    d_prev = floor_gamma_star(R) + D0
    j, junkL1, c0 = initial_junk(R, d_prev)
    minus2 = 1 if c0 == -2 else 0
    front0 = max(j) if j else -1
    last_fed = 0 if (d_prev + 0 <= R - 1) else None   # row 0 is always fed-ish
    i0 = None
    prog = []
    for i in range(1, R - 1):
        d = floor_gamma_star(R + i) + D0
        delta = d - d_prev
        assert delta in (0, 1), (i, d, d_prev)
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
        if d + i <= R - 1:
            w[0] = w.get(0, 0) + g[d + i]; fed = True
        if delta == 1 and d - 1 + i <= R - 1:
            w[0] = w.get(0, 0) + g[d - 1 + i]
            w[1] = w.get(1, 0) + g[d - 1 + i]; fed = True
        jn = {}
        junkL1 = 0
        c0 = 0
        if w:
            ts = sorted(w)
            ta, tb = ts[0], ts[-1]
            assert tb <= d, ("support beyond bottom cell", i, tb, d)
            P = comb(d - 1, ta)
            Pprev = comb(d - 1, ta - 1) if ta >= 1 else 0
            for t in range(ta, tb + 1):
                v = w.get(t, 0)
                if v:
                    lo, hi = -2 * P, 2 * Pprev
                    u = min(hi, max(lo, v))
                    assert (v - u) % 2 == 0 and v % 2 == 0, ("parity", i, t)
                    if u != v:
                        jn[t] = v - u
                        junkL1 += abs(v - u)
                    if t == 0:
                        c0 = u
                Pprev = P
                P = P * (d - 1 - t) // (t + 1) if t < d - 1 else 0
        assert d not in jn, ("DEATH during feed phase -- Theorem B violated?",
                             R, D0, i)
        if c0 == -2:
            minus2 += 1
        j = jn
        d_prev = d
        if fed:
            last_fed = i
        else:
            # first autonomous row reached; feed can never fire again because
            # d_i + i is strictly increasing (certified below).
            i0 = i
            # d_i + i is strictly increasing, so d_{i'} + i' >= R + 1 for all
            # i' > i: neither feed branch (d+i <= R-1; d+i <= R with delta=1)
            # can ever fire again.
            assert d + i >= R, ("autonomy check failed", i, d)
            break
        if i % flush_every == 0:
            prog.append({"i": i, "d": d, "front": (max(j) if j else None),
                         "junkL1_bits": junkL1.bit_length(),
                         "elapsed_s": round(time.time() - t0, 1)})
            if outpath:
                json.dump({"R": R, "D0": D0, "partial_row": i,
                           "progress": prog}, open(outpath + ".partial", "w"))
    assert i0 is not None, "never reached autonomy (impossible)"
    m = max(j) if j else -1
    state = {
        "R": R, "D0": D0,
        "i_feed_obs": last_fed,            # last fed row (D2 convention)
        "first_autonomous_row": i0,        # = i0; state below is j_{i0}
        "d": d_prev,                       # degree at row i0
        "m": m,
        "j0": str(j.get(0, 0)),
        "minus2_rows_so_far": minus2,
        "junkL1_bits": junkL1.bit_length(),
        "N_all_negative": all(v < 0 for v in j.values()),
        "C_contiguous_0_to_m": (sorted(j) == list(range(0, m + 1))),
        "elapsed_s": round(time.time() - t0, 1),
        "junk": {str(t): str(v) for t, v in sorted(j.items())},
    }
    if outpath:
        json.dump(state, open(outpath, "w"))
        try:
            os.remove(outpath + ".partial")
        except OSError:
            pass
    return state


if __name__ == "__main__":
    R, D0 = int(sys.argv[1]), int(sys.argv[2])
    out = (sys.argv[3] if len(sys.argv) > 3 else
           os.path.join(RESULTS,
                        "amm12592_S_feedend_R%d_D0%d_boxeph.json" % (R, D0)))
    st = extract_feedend(R, D0, out)
    print(json.dumps({k: v for k, v in st.items() if k != "junk"}, indent=1))
