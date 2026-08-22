"""AMM 12592 / post-THM-3302: close epoch R = 128 at the gamma* FLOOR profile
d_i = floor(gamma*(128+i)), D0 = 0  (gamma* = log(phi)/log(sqrt 5)).

This is the R = 64 winning recipe (amm12592_r64_floor_solver_boxeph.py),
extended:
  * l1deg beam:  rank states by (residual degree, FULL residual L1) -- the key
    that closed R = 64 (winning trajectories shed L1 monotonically).
  * EXACT capacity pruning (new): the residual after row i must satisfy
        sigma_i = sum_{j>i} x^{j-i-1} Delta_j,
    and for any admissible degree-d block  |[x^m] Delta| <= C(d,m) 2^m
    (proof: [x^m] sum_k delta_k x^{d-k}(1-x)^k bounded by
     sum_{k>=d-m} C(d,k) C(k,d-... ) = C(d,d-m) 2^m = C(d,m) 2^m).
    Hence |sigma_i[m]| <= cap_i[m] with the recursion
        cap_{R-2} = B_{R-1},   cap_i[m] = B_{i+1}[m] + cap_{i+1}[m-1],
        B_j[m] = C(d_j, m) 2^m  (m <= d_j).
    Also |sigma_i(1)| <= R-1-i since Delta_j(1) = delta_{j,0} = +-1 forced.
    Any state violating these can NEVER close: pruning is exact, so it only
    frees beam slots for viable states (still: search negatives are never
    infeasibility evidence, THM-3029).
  * DEEPENED ENDGAME: exhaustive completion over the last L = 2..5 rows
    (rows R-L..R-2 steered over small target grids with capacity pruning,
    row R-1 by exact Bernstein decode), from beam snapshots at rows
    R-6..R-3.
  * checkpoints every 30 rows to amm12592_r128_beam_checkpoint.json.

Exact arithmetic only (int); floats never touch anything.
The gamma* proxy GS floors EQUAL the true gamma* floors for all m <= 512
(proven upstream); R = 128 needs m <= 255.

Usage:
  python3 amm12592_r128_floor_solver_boxeph.py --hunt            # full escalation
  python3 amm12592_r128_floor_solver_boxeph.py --hunt --phases 400
  python3 amm12592_r128_floor_solver_boxeph.py --ref64           # R=64 engine validation
  python3 amm12592_r128_floor_solver_boxeph.py                   # referee saved witness
"""
import argparse, json, os, sys, time
from fractions import Fraction
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
GS = (5979874356654402, 10**16)

def prof(R, g1=GS[0], g2=GS[1], D0=0):
    return [(g1 * (R + i)) // g2 + D0 for i in range(R)]

# ---------------------------------------------------------------------------
# core cell algebra (verbatim semantics from amm12592_r64_floor_solver_boxeph)
# ---------------------------------------------------------------------------
_QK = [[1]]
def qk(k):
    while len(_QK) <= k:
        r = _QK[-1]
        nr = [0] * (len(r) + 1)
        for i, c in enumerate(r):
            nr[i] += c
            nr[i + 1] -= c
        _QK.append(nr)
    return _QK[k]

def qpow(n):
    return list(qk(n))

def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a

_CB = {}
def cb(d):
    r = _CB.get(d)
    if r is None:
        r = [comb(d, k) for k in range(d + 1)]
        _CB[d] = r
    return r

def step(sigma, d, targets):
    """One residual step at degree d.  targets[j] steers sigma_next[j]
    (None = greedy clamp).  Returns (deltas, sigma_next) or None.
    Semantics identical to the R=64 solver; slice ops for speed."""
    if not sigma or abs(sigma[0]) != 1:
        return None
    res = list(sigma) + [0] * (d + 2 - len(sigma) if len(sigma) < d + 2 else 1)
    if len(res) < d + 1:
        res += [0] * (d + 1 - len(res))
    caps = cb(d)
    deltas = [0] * (d + 1)
    v = sigma[0]
    deltas[d] = v                      # forced: B_{d,d} = (1-x)^d, offset 0
    t = qk(d)
    n = min(len(t), len(res))
    res[:n] = [a - v * b for a, b in zip(res[:n], t[:n])]
    for j, tgt in enumerate(targets):
        k = d - 1 - j
        if k < 0:
            return None
        cap = caps[k]
        if tgt is None:
            want = res[j + 1]
            want = max(-cap, min(cap, want))
            if (want - cap) % 2:
                want = want - 1 if want - 1 >= -cap else want + 1
        else:
            want = res[j + 1] - tgt
        if abs(want) > cap or (want - cap) % 2:
            return None
        deltas[k] = want
        if want:
            t = qk(k)
            off = d - k
            res[off:off + k + 1] = [a - want * b
                                    for a, b in zip(res[off:off + k + 1], t)]
    for k in range(d - 1 - len(targets), -1, -1):
        cap = caps[k]
        want = res[d - k]
        vv = max(-cap, min(cap, want))
        if (vv - cap) % 2:
            vv = vv - 1 if vv - 1 >= -cap else vv + 1
        deltas[k] = vv
        if vv:
            t = qk(k)
            off = d - k
            res[off:off + k + 1] = [a - vv * b
                                    for a, b in zip(res[off:off + k + 1], t)]
    if res[0] != 0:
        return None
    return deltas, trim(res[1:])

def final_decode(sigma, d):
    """Exact Bernstein-d decode of sigma; deltas iff admissible + exact."""
    if len(sigma) - 1 > d:
        return None
    res = list(sigma) + [0] * (d + 1 - len(sigma))
    caps = cb(d)
    de = [0] * (d + 1)
    for k in range(d, -1, -1):
        cap = caps[k]
        want = res[d - k]
        if abs(want) > cap or (want - cap) % 2:
            return None
        de[k] = want
        if want:
            t = qk(k)
            off = d - k
            res[off:off + k + 1] = [a - want * b
                                    for a, b in zip(res[off:off + k + 1], t)]
    return de if not trim(res) else None

# ---------------------------------------------------------------------------
# exact capacity pruning
# ---------------------------------------------------------------------------
def coeff_caps(R, d):
    """caps[i][m] = exact upper bound on |[x^m] sigma_i| over all admissible
    completions by rows i+1..R-1.  caps[i] has finite support; any state
    exceeding it (or the degree support) can never close."""
    B = [[comb(d[j], m) << m for m in range(d[j] + 1)] for j in range(R)]
    caps = [None] * (R - 1)
    caps[R - 2] = list(B[R - 1])
    for i in range(R - 3, -1, -1):
        prev = caps[i + 1]
        bj = B[i + 1]
        n = max(len(bj), len(prev) + 1)
        cur = [0] * n
        for m, v in enumerate(bj):
            cur[m] += v
        for m, v in enumerate(prev):
            cur[m + 1] += v
        caps[i] = cur
    return caps

def sig_ok(sig, cap, rem):
    """Exact viability: coefficientwise |sigma| <= cap, and |sigma(1)| <= rem
    (rem = number of remaining rows; each contributes delta_{.,0} = +-1)."""
    if len(sig) > len(cap):
        return False
    for m, c in enumerate(sig):
        if c and (c if c > 0 else -c) > cap[m]:
            return False
    s = sum(sig)
    return -rem <= s <= rem

# ---------------------------------------------------------------------------
# exact verification (search-independent)
# ---------------------------------------------------------------------------
def admissible(delta, dd):
    if len(delta) - 1 > dd:
        return False
    return all(abs(v) <= comb(dd, k) and (v - comb(dd, k)) % 2 == 0
               for k, v in enumerate(delta))

def epoch_identity(R, sol, d):
    acc = [0] * (R + max(d) + 2)
    for i, de in enumerate(sol):
        for k, v in enumerate(de):
            if v:
                t = qk(k)
                off = i + d[i] - k
                for s in range(k + 1):
                    acc[off + s] += v * t[s]
    return trim(acc) == trim(qpow(R - 1))

def verify_witness(R, blocks, profile):
    ok_adm = all(admissible(blocks[i], profile[i]) for i in range(R))
    ok_id = epoch_identity(R, blocks, profile)
    return ok_adm, ok_id

def eff(R, d):
    return max(Fraction(d[i], R + i) for i in range(R))

# ---------------------------------------------------------------------------
# beam engine
# ---------------------------------------------------------------------------
CKPT = os.path.join(HERE, "amm12592_r128_beam_checkpoint.json")
WITNESS = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")

def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f)
    os.replace(tmp, path)

def l1(sig):
    return sum(c if c > 0 else -c for c in sig)

def run_beam(R, d, caps, beam, ctrl, span, dedup, snap_rows, phase_tag,
             ckpt_path, log, prune=True):
    """Rows 0..R-3 by l1deg beam.  Returns (states@row R-3, snaps, hist)."""
    opts = [None] + list(range(-span, span + 1))
    grids = [tg for tg in product(opts, repeat=ctrl) if tg[0] in (1, -1)]
    states = [([], qpow(R - 1))]
    snaps, hist = {}, []
    t0 = time.time()
    for i in range(R - 2):
        di = d[i]
        cap_i = caps[i] if prune else None
        rem_i = R - 1 - i
        nxt = []
        for acc, sig in states:
            for tg in grids:
                r = step(sig, di, tg)
                if r is None:
                    continue
                de, ns = r
                if not ns or (ns[0] != 1 and ns[0] != -1):
                    continue
                if cap_i is not None and not sig_ok(ns, cap_i, rem_i):
                    continue
                nxt.append((acc + [de], ns))
        if not nxt:
            log(f"  row {i:3d}: BEAM DIED (all children pruned/failed)")
            hist.append({"row": i, "n": 0})
            return None, snaps, hist
        nxt.sort(key=lambda st: (len(st[1]), l1(st[1])))
        seen, uniq = set(), []
        for a, sg in nxt:
            key = tuple(sg[:dedup])
            if key in seen:
                continue
            seen.add(key)
            uniq.append((a, sg))
        states = uniq[:beam]
        b = states[0][1]
        l1s = [l1(sg) for _, sg in states[:5]]
        hist.append({"row": i, "n": len(states), "deg": len(b) - 1,
                     "l1_top5": [str(v) for v in l1s]})
        log(f"  row {i:3d}: kept={len(states):5d} bestdeg={len(b)-1:3d} "
            f"top5L1={l1s} t={time.time()-t0:7.1f}s")
        if i in snap_rows:
            snaps[i] = list(states)
        if ckpt_path and (i % 30 == 29 or i == R - 3):
            write_json(ckpt_path, {
                "R": R, "phase": phase_tag, "row": i,
                "elapsed_s": round(time.time() - t0, 1), "hist": hist,
                "states": [{"acc": a, "sigma": sg}
                           for a, sg in states[:40]]})
            log(f"    [checkpoint written at row {i}]")
    return states, snaps, hist

def endgame(R, d, caps, states, L, ec, es, budget, log, max_states=None):
    """Exhaustive completion of the last L rows from states at row R-1-L:
    rows R-L..R-2 steered over the (ec, es) grid with capacity pruning,
    row R-1 by exact Bernstein decode.  Returns (full solution rows, msg)."""
    eopts = [None] + list(range(-es, es + 1))
    grid = [tg for tg in product(eopts, repeat=ec) if tg[0] in (1, -1)]
    cnt = [0]
    dlast = d[R - 1]

    def dfs(sig, j):
        if j == R - 1:
            cnt[0] += 1
            fin = final_decode(sig, dlast)
            return None if fin is None else [fin]
        cap_j, rem_j, dj = caps[j], R - 1 - j, d[j]
        for tg in grid:
            if cnt[0] >= budget:
                return "BUDGET"
            cnt[0] += 1
            r = step(sig, dj, tg)
            if r is None:
                continue
            de, ns = r
            if not ns or (ns[0] != 1 and ns[0] != -1):
                continue
            if not sig_ok(ns, cap_j, rem_j):
                continue
            sub = dfs(ns, j + 1)
            if sub == "BUDGET":
                return "BUDGET"
            if sub is not None:
                return [de] + sub
        return None

    use = states if max_states is None else states[:max_states]
    t0 = time.time()
    for si, (acc, sig) in enumerate(use):
        if si % 25 == 0:
            log(f"    endgame L={L} ({ec},{es}): state {si}/{len(use)} "
                f"nodes={cnt[0]} t={time.time()-t0:6.1f}s")
        res = dfs(sig, R - L)
        if res == "BUDGET":
            log(f"    endgame L={L}: NODE BUDGET {budget} exhausted at state "
                f"{si}/{len(use)}")
            return None, f"budget@state{si}"
        if res is not None:
            log(f"    endgame L={L}: SOLVED from state {si} "
                f"(nodes={cnt[0]})")
            return acc + res, "SOLVED"
    log(f"    endgame L={L}: exhausted {len(use)} states (nodes={cnt[0]})")
    return None, "exhausted"

ENDGAME_PLANS = [(2, 3, 6, 3_000_000, None),
                 (3, 3, 6, 1_200_000, 400),
                 (4, 2, 5, 1_200_000, 250),
                 (5, 2, 4, 1_200_000, 150)]

def hunt_phase(R, d, caps, beam, ctrl, span, dedup, log, ckpt_path, tag):
    snap_rows = {R - 6, R - 5, R - 4, R - 3}
    log(f"PHASE {tag}: R={R} beam={beam} ctrl={ctrl} span={span} rank=l1deg "
        f"dedup={dedup} (deterministic)")
    states, snaps, hist = run_beam(R, d, caps, beam, ctrl, span, dedup,
                                   snap_rows, tag, ckpt_path, log)
    if states is None:
        return None, hist
    for (L, ec, es, budget, mx) in ENDGAME_PLANS:
        t = R - 1 - L
        st = snaps.get(t)
        if st is None:
            continue
        log(f"  endgame plan L={L} from row {t} snapshot "
            f"({len(st)} states, grid ctrl={ec} span={es})")
        sol, msg = endgame(R, d, caps, st, L, ec, es, budget, log,
                           max_states=mx)
        if sol is not None:
            return sol, hist
    return None, hist

def finalize(R, d, sol, search_meta, log, out=WITNESS):
    a, b = verify_witness(R, sol, d)
    log(f"  EXACT VERIFY: admissible={a} identity={b}")
    assert a and b, "witness failed exact verification -- NOT saved"
    write_json(out, {"R": R, "profile": d, "blocks": sol, "verified": True,
                     "H": [1],
                     "source_label": "gamma* floor (direct l1deg beam, no lift)",
                     "search": search_meta})
    log(f"  WITNESS WRITTEN: {out}")

def hunt(args):
    R = args.R
    d = prof(R)
    caps = coeff_caps(R, d)
    log = lambda s: print(s, flush=True)
    log(f"HUNT R={R} gamma* floor profile, d[0]={d[0]} d[-1]={d[-1]} "
        f"eff={float(eff(R, d)):.6f}")
    phases = [int(x) for x in args.phases.split(",")]
    for beam in phases:
        tag = f"beam{beam}"
        sol, hist = hunt_phase(R, d, caps, beam, args.ctrl, args.span,
                               args.dedup, log, CKPT if R == 128 else None,
                               tag)
        if sol is not None:
            finalize(R, d, sol, {"beam": beam, "ctrl": args.ctrl,
                                 "span": args.span, "rank": "l1deg",
                                 "dedup": args.dedup, "prune": "caps+eval1",
                                 "phase": tag}, log)
            log("HUNT: CLOSED")
            return 0
        log(f"PHASE {tag}: no witness (search negative -- artefact, "
            f"never infeasibility evidence)")
    log("HUNT: all phases exhausted without witness")
    return 1

def ref64(args):
    """Engine validation: reproduce the R = 64 floor closure with the exact
    same recipe (beam 400) through this engine, verify, do NOT overwrite the
    saved witness.  Also logs the winning L1 trajectory for shape reference."""
    R = 64
    d = prof(R)
    caps = coeff_caps(R, d)
    log = lambda s: print(s, flush=True)
    log(f"REF64 engine validation: beam=400 ctrl=2 span=2 l1deg + caps prune")
    sol, hist = hunt_phase(R, d, caps, 400, 2, 2, 999, log, None, "ref64")
    if sol is None:
        log("REF64: FAILED -- engine adaptation broken (R=64 must close)")
        return 1
    a, b = verify_witness(R, sol, d)
    log(f"REF64: verify admissible={a} identity={b}")
    assert a and b
    log("REF64: PASS (engine reproduces the R=64 floor closure)")
    return 0

def referee():
    R = 128
    d = prof(R)
    print("AMM 12592 / R = 128 gamma* floor-profile closure referee")
    if not os.path.exists(WITNESS):
        print("no witness JSON present; run --hunt")
        return 1
    with open(WITNESS) as f:
        w = json.load(f)
    assert w["R"] == R and w["profile"] == d
    a, b = verify_witness(R, w["blocks"], d)
    print(f"witness blocks: {len(w['blocks'])}; admissible={a} identity={b}")
    e = eff(R, d)
    print(f"effective rate max_i d_i/(128+i) = {float(e):.6f} "
          f"(gamma* = 0.597987)")
    assert a and b
    print("ALL R128 FLOOR-WITNESS CHECKS PASSED")
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hunt", action="store_true")
    ap.add_argument("--ref64", action="store_true")
    ap.add_argument("--R", type=int, default=128)
    ap.add_argument("--phases", default="400,1000,2000")
    ap.add_argument("--ctrl", type=int, default=2)
    ap.add_argument("--span", type=int, default=2)
    ap.add_argument("--dedup", type=int, default=999)
    args = ap.parse_args()
    if args.ref64:
        sys.exit(ref64(args))
    if args.hunt:
        sys.exit(hunt(args))
    sys.exit(referee())

if __name__ == "__main__":
    main()
