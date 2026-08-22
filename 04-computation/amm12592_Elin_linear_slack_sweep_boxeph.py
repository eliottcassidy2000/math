"""LANE D2 (E-lin) -- exact linear-slack sweep of plain rule A.

For each dyadic R, run the certified fast T6 junk-flow engine at
D0 = ceil(eps*R) for eps in {1/32, 1/16} plus hazard probes (large D0 at
R = 512 to test whether rule A itself is monotone in D0 -- NOT implied by
the LIFT lemma, which lifts witnesses, not rule behavior).

Everything exact (engine is int-only, certified bit-identical to the slow
path at 3 implementations).  Partial results flushed after every run.

Validity guard: the engine's initial_junk (T4 closed form) requires
d_0 <= R-2 (the hockey-stick step needs deg <= R-2); asserted per run.
"""
import sys, os, json, importlib.util, contextlib, io, time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")
OUT = os.path.join(RESULTS, "amm12592_Elin_sweep_boxeph.json")

spec = importlib.util.spec_from_file_location(
    "fast", os.path.join(HERE, "amm12592_transient_fast_junkflow_boxeph.py"))
fast = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(fast)
floor_gamma_star = fast.floor_gamma_star


def ceil_div(a, b):
    return -((-a) // b)


def one(R, D0, results):
    d0 = floor_gamma_star(R) + D0
    assert d0 <= R - 2, ("engine validity: need d_0 <= R-2", R, D0, d0)
    t0 = time.time()
    res = fast.run_fast(R, D0, keep_rows=False)
    rec = {"R": R, "D0": D0, "d0": d0, "eps_num": D0, "eps_den": R,
           "outcome": res["outcome"], "capture_row": res["capture_row"],
           "minus2_rows": res["minus2_rows"], "front0": res["front0"],
           "T6b_bound": res["T6b_bound"],
           "front0_eq_R-2-d0": (res["front0"] == R - 2 - d0),
           "elapsed_s": round(time.time() - t0, 1)}
    results.append(rec)
    json.dump(results, open(OUT, "w"), indent=1)
    print(f"R={R:6d} D0={D0:5d} (eps={D0/R:.5f} d0={d0}): "
          f"{res['outcome']}  capture={res['capture_row']} "
          f"minus2={res['minus2_rows']} front0={res['front0']} "
          f"({rec['elapsed_s']}s)", flush=True)
    return rec


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "main"
    results = []
    if os.path.exists(OUT):
        results = json.load(open(OUT))
    done = {(r["R"], r["D0"]) for r in results}
    if which == "main":
        jobs = []
        for R in (128, 256, 512, 1024, 2048, 4096, 8192):
            for eps_den in (32, 16):
                jobs.append((R, ceil_div(R, eps_den)))
        # hazard probes: is rule A monotone in D0?  large-slack behavior
        jobs += [(512, 64), (512, 128), (512, 204),
                 (1024, 256), (1024, 409),
                 (2048, 256), (2048, 512)]
        for R, D0 in jobs:
            if (R, D0) in done:
                continue
            one(R, D0, results)
    elif which == "big":
        # R = 16384 at eps = 1/32 and 1/16 (background-friendly)
        for R, D0 in ((16384, 512), (16384, 1024)):
            if (R, D0) in done:
                continue
            one(R, D0, results)
    print("DONE", which, flush=True)
