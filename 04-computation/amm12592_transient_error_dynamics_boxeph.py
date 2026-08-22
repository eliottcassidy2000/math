"""ANGLE B2 (transient theorem) -- exact error-coordinate dynamics of attractor rule A.

Everything exact (int). p = x, q = 1-x, E_m = -1+x+...+x^m (E_{-1}=0).

RULE A (plain, THM-3329):
  sigma_{-1} = q^{R-1};  row i:  ideal_i = sigma_{i-1} - x E_{R-2-i}
  Delta_i = AdmClamp_{d_i}(trunc_{d_i}(ideal_i));  need x | sigma_{i-1}-Delta_i;
  sigma_i = (sigma_{i-1}-Delta_i)/x.  Closure iff sigma_{R-1}=0.

ERROR COORDINATES (Lemma T2, this session):
  e_i := sigma_i - E_{R-2-i}   (E_{-1}:=0), so e_{-1} = q^{R-1}-E_{R-1} = 2 G_R.
  For 0 <= i <= R-2 the rule is EXACTLY conjugate to
     w_t   = decode_{d_i}(trunc_{d_i}(e_{i-1}))_t          (t = 0..d_i)
     c_t   = clamp of w_t into the asymmetric even box
                 [-2 binom(d_i-1,t), +2 binom(d_i-1,t-1)]
             (parity fix relative to the shifted ballot corner; at dyadic R
              the parity step provably NEVER fires -- Lemma T3)
     survival iff [x^0](e_{i-1} - c_poly) = 0, equivalently e_{i-1}[0] in {0,2}
     e_i   = (e_{i-1} - c_poly)/x
  with Delta_i = (2x-1) + c_poly, i.e. cellwise delta_t = b_t + c_t where
  b_t = binom(d-1,t)-binom(d-1,t-1) is the ballot vector (decode of 2x-1).
  Row R-1: Delta_{R-1} = -1 + c, decode of -1 is the corner -binom(d,t),
  c_t in [0, 2 binom(d,t)] even; closure iff e_{R-1} = 0.

This module: (1) runs the error dynamics with full instrumentation,
(2) proves the conjugacy by exact block-for-block comparison against rule A
    (adm_clamp implementation) and against stored witness JSONs,
(3) dumps per-row traces (junk front, bulge depth, clamp counts, ballot word)
    for the transient-bound analysis.  Partial results saved every 16 rows.
"""
import sys, os, json, time, io, contextlib, importlib.util
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.path.join(os.path.dirname(HERE), "05-knowledge", "results")

from amm12592_allR_family_toolbox_boxeph import (padd, psub, pneg, pmul, pshift,
    pscale, qpow, bern_to_poly, poly_to_bern, ballot, admissible, epoch_sum)

spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)

def Em(m): return [-1] + [1]*m if m >= 0 else []

def exact_profile(R, D0=0):
    return [iref.floor_gamma_star(R + i) + D0 for i in range(R)]

def two_G(R):
    """e_{-1} = q^{R-1} - E_{R-1}, exact."""
    return psub(qpow(R - 1), Em(R - 1))

# ---------------------------------------------------------------- rule A ref
def adm_clamp_ref(P, d):
    """Reference AdmClamp (tozero parity fix), identical to rule-A script."""
    cells = poly_to_bern(P, d)
    out = []
    for t, v in enumerate(cells):
        cap = comb(d, t)
        w = min(cap, max(-cap, v))
        if (w - cap) % 2:
            w = w - 1 if w > 0 else w + 1
        out.append(w)
    return out

def run_ruleA(R, D0=0):
    """Plain rule A with the EXACT floor profile; returns (blocks|None, msg)."""
    pr = exact_profile(R, D0)
    sig = qpow(R - 1)
    blocks = []
    for i in range(R):
        d = pr[i]
        ideal = psub(sig, pshift(Em(R - 2 - i), 1))
        delta = adm_clamp_ref(ideal[:d + 1], d)
        t = psub(sig, bern_to_poly(delta, d))
        if t and t[0] != 0:
            return None, f"DIE row {i} const bits {abs(t[0]).bit_length()}"
        sig = t[1:] if t else []
        blocks.append(delta)
    return (blocks, "CLOSED") if sig == [] else (None, "residual nonzero")

# ------------------------------------------------- error-coordinate dynamics
def ebox(d, t):
    """Asymmetric even box for c_t at degree d: (lo, hi)."""
    return (-2 * comb(d - 1, t), 2 * (comb(d - 1, t - 1) if t >= 1 else 0))

def run_error_dynamics(R, D0=0, trace_path=None, flush_every=16):
    """Error-coordinate run.  Returns dict with blocks (delta cells), trace,
    outcome.  Asserts (dyadic R): parity fix never needed on any cell."""
    pr = exact_profile(R, D0)
    e = two_G(R)
    dyadic = (R & (R - 1)) == 0
    blocks, trace = [], []
    parity_fires = 0
    t0 = time.time()
    outcome = None
    for i in range(R):
        d = pr[i]
        last = (i == R - 1)
        # corner: ballot b (rows <= R-2)  /  bottom corner -binom (last row)
        if not last:
            corner = [comb(d - 1, t) - (comb(d - 1, t - 1) if t else 0)
                      for t in range(d + 1)]
        else:
            corner = [-comb(d, t) for t in range(d + 1)]
        w = poly_to_bern(e[:d + 1], d)  # handles len <= d+1
        c = []
        nclamped = 0; tmin = None; tmax = None; junkL1 = 0
        for t in range(d + 1):
            cap = comb(d, t)
            if not last:
                lo, hi = ebox(d, t)
            else:
                lo, hi = 0, 2 * cap
            v = w[t]
            u = min(hi, max(lo, v))
            # parity fix (translate of tozero on delta = corner + c):
            if (u - lo) % 2:
                parity_fires += 1
                delta_v = corner[t] + u
                u = (u - 1) if delta_v > 0 else (u + 1)
            if u != v:
                nclamped += 1
                junkL1 += abs(v - u)
                if tmin is None: tmin = t
                tmax = t
            c.append(u)
        cpoly = bern_to_poly(c, d)
        rem = psub(e, cpoly)
        e0 = rem[0] if rem else 0
        # instrumentation on e BEFORE shift
        row = {"i": i, "d": d, "nclamped": nclamped, "tmin": tmin, "tmax": tmax,
               "junkL1_bits": junkL1.bit_length(),
               "c0": c[0] if c else 0,
               "e_in0": (e[0] if e else 0)}
        if e0 != 0:
            row["death_const_bits"] = abs(e0).bit_length()
            trace.append(row)
            outcome = {"status": "DIE", "row": i,
                       "const_bits": abs(e0).bit_length()}
            break
        e = rem[1:] if rem else []
        # post-shift stats
        supp = [j for j, v in enumerate(e) if v]
        row["degE"] = (supp[-1] if supp else -1)
        row["minsupp"] = (supp[0] if supp else -1)
        row["L1_bits"] = sum(abs(v) for v in e).bit_length()
        # junk front: least j with |e[j]| > 2
        front = next((j for j, v in enumerate(e) if abs(v) > 2), -1)
        row["front"] = front
        blocks.append([corner[t] + c[t] for t in range(d + 1)])
        trace.append(row)
        if trace_path and (i % flush_every == 0 or i == R - 1):
            json.dump({"R": R, "D0": D0, "partial_row": i, "trace": trace,
                       "parity_fires": parity_fires,
                       "elapsed_s": round(time.time() - t0, 1)},
                      open(trace_path, "w"))
    if outcome is None:
        outcome = {"status": "CLOSED"} if e == [] else \
                  {"status": "OPEN_RESIDUAL",
                   "L1_bits": sum(abs(v) for v in e).bit_length()}
    res = {"R": R, "D0": D0, "outcome": outcome, "trace": trace,
           "parity_fires": parity_fires, "dyadic": dyadic,
           "elapsed_s": round(time.time() - t0, 1),
           "blocks": blocks if outcome["status"] == "CLOSED" else None}
    if trace_path:
        json.dump({k: v for k, v in res.items() if k != "blocks"},
                  open(trace_path, "w"))
    return res

# ------------------------------------------------------------- verification
def verify_conjugacy(R, D0=0):
    """Exact block-for-block equality: error dynamics vs plain rule A."""
    ra, msg = run_ruleA(R, D0)
    ed = run_error_dynamics(R, D0)
    if ra is None or ed["outcome"]["status"] != "CLOSED":
        return {"R": R, "D0": D0, "ruleA": msg, "err_dyn": ed["outcome"],
                "match": (ra is None and ed["outcome"]["status"] != "CLOSED")}
    same = (len(ra) == len(ed["blocks"]) and
            all(ra[i] == ed["blocks"][i] for i in range(R)))
    # full witness verification of the error-dynamics blocks
    pr = exact_profile(R, D0)
    adm = all(admissible(ed["blocks"][i], pr[i]) for i in range(R))
    idt = epoch_sum(R, ed["blocks"], pr) == qpow(R - 1)
    return {"R": R, "D0": D0, "ruleA": msg, "err_dyn": "CLOSED",
            "blocks_equal": same, "admissible": adm, "identity": idt,
            "parity_fires": ed["parity_fires"],
            "match": same and adm and idt}

def verify_vs_witness_json(R, D0, path):
    """Error dynamics blocks must equal a stored rule-A witness JSON."""
    ed = run_error_dynamics(R, D0,
        trace_path=os.path.join(RESULTS,
            f"amm12592_transient_trace_R{R}_D0{D0}_boxeph.json"))
    W = json.load(open(os.path.join(HERE, path)))
    if ed["outcome"]["status"] != "CLOSED":
        return {"R": R, "D0": D0, "err_dyn": ed["outcome"], "match": False}
    same = all(ed["blocks"][i] == W["blocks"][i] for i in range(R))
    return {"R": R, "D0": D0, "err_dyn": "CLOSED", "witness_equal": same,
            "parity_fires": ed["parity_fires"], "match": same,
            "trace_rows": len(ed["trace"])}

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "small"
    if mode == "small":
        out = []
        for R in (8, 16, 32, 64, 128):
            v = verify_conjugacy(R, 0)
            print(f"R={R:4d} D0=0: {v}", flush=True)
            out.append(v)
        json.dump(out, open(os.path.join(RESULTS,
            "amm12592_transient_conjugacy_small_boxeph.json"), "w"), indent=1)
    elif mode == "R256":
        v = verify_vs_witness_json(256, 1,
            "amm12592_witness_R256_ruleA_D0_1_boxeph.json")
        print(v, flush=True)
        json.dump(v, open(os.path.join(RESULTS,
            "amm12592_transient_conjugacy_R256_boxeph.json"), "w"), indent=1)
    elif mode == "R512":
        v = verify_vs_witness_json(512, 8,
            "amm12592_witness_R512_ruleA_D0_8_boxeph.json")
        print(v, flush=True)
        json.dump(v, open(os.path.join(RESULTS,
            "amm12592_transient_conjugacy_R512_boxeph.json"), "w"), indent=1)
    elif mode == "die512":
        # death forensics: D0 = 5, 6, 7 at R = 512 (unknown region of the scan)
        for D0 in (5, 6, 7):
            r = run_error_dynamics(512, D0, trace_path=os.path.join(RESULTS,
                f"amm12592_transient_trace_R512_D0{D0}_boxeph.json"))
            print(f"R=512 D0={D0}: {r['outcome']} ({r['elapsed_s']}s)",
                  flush=True)
    elif mode == "die256":
        r = run_error_dynamics(256, 0, trace_path=os.path.join(RESULTS,
            "amm12592_transient_trace_R256_D00_boxeph.json"))
        print(f"R=256 D0=0: {r['outcome']} ({r['elapsed_s']}s)", flush=True)
    elif mode == "R128trace":
        r = run_error_dynamics(128, 0, trace_path=os.path.join(RESULTS,
            "amm12592_transient_trace_R128_D00_boxeph.json"))
        print(f"R=128 D0=0: {r['outcome']} ({r['elapsed_s']}s)", flush=True)
