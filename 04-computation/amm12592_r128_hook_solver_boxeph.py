"""AMM 12592 R=128: ATTRACTOR-HOOK beam solver.

Anatomy of the R=64 slim winner (amm12592_r64_attractor_anatomy_boxeph.py):
  * late residuals ride the attractor sigma_i = E_m := -1 + x + ... + x^m,
    m = R-2-i, with ride blocks Delta_j = -1 + 2x = p - q (ballot backbone)
    and final row Delta_{R-1} = -1 (negative full box);
  * the snap (row 54 at R=64) is Delta_i = sigma_{i-1} - p E_{R-2-i} being
    Bernstein-d_i admissible even though sigma_{i-1} != E_{R-1-i}.

KEY EXACT FACT: every residual satisfies sigma_{i-1} == 1 + x + ... +
x^{R-1-i} (mod 2) coefficientwise (each admissible block == 1 mod 2, since
sum_k C(d,k) B_{d,k} = 1).  The Bernstein-d decode matrix is triangular with
unit diagonal, so decode(junk) == 0 (mod 2) cellwise iff junk == 0 (mod 2)
coefficientwise.  Writing sigma_{i-1} = E_{R-1-i} + junk, the parity of junk
is ALWAYS even: the parity condition never blocks a snap.  The hook
    Delta_i := sigma_{i-1} - p (sgn E_{R-2-i}),  decode at d_i,
succeeds iff the CAPS hold; mid-coefficient junk ~1e14 faces caps ~1e40, so
hooks can fire long before the residual L1 collapses.  This removes the need
for the L1-collapse phase that the plain l1deg beam cannot reproduce at
R = 128 (its plateau ~1e14 = square of R=64's ~1e7 plateau).

Search layer only; every produced witness is verified exactly and
independently (identity in Z[x] + per-row admissibility) before saving.
"""
import argparse, json, os, sys, time
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng

WITNESS = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")
CKPT = os.path.join(HERE, "amm12592_r128_hook_checkpoint.json")

def E(m):
    """E_m = -1 + x + ... + x^m (E_0 = [-1])."""
    return [-1] + [1] * m

def ride_blocks(R, d, sgn):
    """Ride decode per row: rows j<R-1 use sgn*(2x-1), row R-1 uses sgn*(-1).
    Returns list rb with rb[j] = cells or None (inadmissible)."""
    rb = [None] * R
    for j in range(R):
        poly = [-sgn, 2 * sgn] if j < R - 1 else [-sgn]
        rb[j] = eng.final_decode(list(poly), d[j])
    return rb

def hook(sig, di, m_next, sgn):
    """Try Delta = sig - p*(sgn*E_{m_next}) as an admissible d_i block."""
    ptau = [0] + [sgn * c for c in E(m_next)]
    n = max(len(sig), len(ptau))
    cand = [(sig[t] if t < len(sig) else 0) - (ptau[t] if t < len(ptau) else 0)
            for t in range(n)]
    cand = eng.trim(cand)
    if len(cand) - 1 > di:
        return None
    return eng.final_decode(cand, di)

def soft_decode_badness(res0, d):
    """Clamped Bernstein-d decode: count cells where cap/parity clamping was
    needed plus a nonzero-leftover flag.  0 == exactly hook-decodable."""
    res = list(res0)
    if len(res) < d + 1:
        res += [0] * (d + 1 - len(res))
    caps = eng.cb(d)
    viol = 0
    for k in range(d, -1, -1):
        cap = caps[k]
        idx = d - k
        want = res[idx] if idx < len(res) else 0
        vv = max(-cap, min(cap, want))
        if (vv - cap) % 2:
            vv = vv - 1 if vv - 1 >= -cap else vv + 1
        if vv != want:
            viol += 1
        if vv:
            t = eng.qk(k)
            off = d - k
            for s in range(k + 1):
                res[off + s] -= vv * t[s]
    res[:d + 1] = [0] * min(len(res), d + 1)
    left = eng.trim(res)
    return viol + (1 if left else 0)

def hook_badness(sig, di, m_next):
    """min over sign of the soft-decode violation count of the hook block;
    0 means the hook fires exactly."""
    best = None
    for sgn in (1, -1):
        ptau = [0] + [sgn * c for c in E(m_next)]
        n = max(len(sig), len(ptau))
        cand = [(sig[t] if t < len(sig) else 0) -
                (ptau[t] if t < len(ptau) else 0) for t in range(n)]
        cand = eng.trim(cand)
        if len(cand) - 1 > di:
            b = 10 ** 6 + (len(cand) - 1 - di)
        else:
            b = soft_decode_badness(cand, di)
        if best is None or b < best:
            best = b
    return best

def assemble(R, d, acc, hook_de, i, sgn, rb):
    """Full block list: rows 0..i-1 = acc, row i = hook, rows i+1..R-1 ride."""
    blocks = list(acc) + [hook_de]
    for j in range(i + 1, R):
        assert rb[j] is not None
        blocks.append(rb[j])
    return blocks

def solve_hook(R, d, caps, beam, ctrl, span, dedup, log, ckpt_path, tag,
               hook_min_row=0, rank="l1deg"):
    opts = [None] + list(range(-span, span + 1))
    grids = [tg for tg in product(opts, repeat=ctrl) if tg[0] in (1, -1)]
    rb = {1: ride_blocks(R, d, 1), -1: ride_blocks(R, d, -1)}
    # rideok[j] = all rows j..R-1 rideable (suffix AND), per sign
    rideok = {}
    for sgn in (1, -1):
        ok = [False] * (R + 1)
        ok[R] = True
        for j in range(R - 1, -1, -1):
            ok[j] = ok[j + 1] and (rb[sgn][j] is not None)
        rideok[sgn] = ok
    first_ride = {s: next((j for j in range(R + 1) if rideok[s][j]), R)
                  for s in (1, -1)}
    log(f"  ride feasibility: sign +1 rideable from row {first_ride[1]}, "
        f"sign -1 from row {first_ride[-1]} (of {R})")
    states = [([], eng.qpow(R - 1))]
    t0 = time.time()
    hook_stats = []
    for i in range(R - 1):
        di = d[i]
        m_next = R - 2 - i
        # --- hook attempt on every parent state ---
        if i >= hook_min_row:
            for si, (acc, sig) in enumerate(states):
                for sgn in (1, -1):
                    if not rideok[sgn][i + 1]:
                        continue
                    de = hook(sig, di, m_next, sgn)
                    if de is not None:
                        log(f"  row {i:3d}: HOOK FIRED (state {si}, sign "
                            f"{sgn:+d}) t={time.time()-t0:.1f}s")
                        blocks = assemble(R, d, acc, de, i, sgn, rb[sgn])
                        return blocks, f"HOOKED@row{i}state{si}sign{sgn}"
        if i == R - 2:
            break
        # --- normal l1deg beam expansion with exact pruning ---
        cap_i, rem_i = caps[i], R - 1 - i
        nxt = []
        for acc, sig in states:
            for tg in grids:
                r = eng.step(sig, di, tg)
                if r is None:
                    continue
                de, ns = r
                if not ns or (ns[0] != 1 and ns[0] != -1):
                    continue
                if not eng.sig_ok(ns, cap_i, rem_i):
                    continue
                nxt.append((acc + [de], ns))
        if not nxt:
            log(f"  row {i:3d}: BEAM DIED")
            return None, f"died@row{i}"
        if rank == "hookdist" and i + 1 <= R - 2:
            dnx, mnx = d[i + 1], R - 3 - i
            nxt.sort(key=lambda st: (hook_badness(st[1], dnx, mnx),
                                     len(st[1]), eng.l1(st[1])))
        else:
            nxt.sort(key=lambda st: (len(st[1]), eng.l1(st[1])))
        seen, uniq = set(), []
        for a, sg in nxt:
            key = tuple(sg[:dedup])
            if key in seen:
                continue
            seen.add(key)
            uniq.append((a, sg))
        states = uniq[:beam]
        b = states[0][1]
        if i % 5 == 0 or i > R - 12:
            extra = ""
            if rank == "hookdist" and i + 1 <= R - 2:
                extra = f" bestbad={hook_badness(b, d[i+1], R-3-i)}"
            log(f"  row {i:3d}: kept={len(states):5d} bestdeg={len(b)-1:3d} "
                f"bestL1={eng.l1(b)}{extra} t={time.time()-t0:6.1f}s")
        if ckpt_path and i % 30 == 29:
            eng.write_json(ckpt_path, {
                "R": R, "tag": tag, "row": i,
                "states": [{"acc": a, "sigma": sg} for a, sg in states[:40]]})
    # states are sigma_{R-3} after the loop (i stopped at R-2 pre-expansion);
    # dump the full population for offline repair, then run the exact banded
    # 2-row repair as the endgame.
    if ckpt_path:
        dump = ckpt_path.replace(".json", "_final_states.json")
        eng.write_json(dump, {
            "R": R, "tag": tag, "row": R - 3,
            "states": [{"acc": a, "sigma": sg} for a, sg in states]})
        log(f"  [final population dumped: {dump} ({len(states)} states)]")
    import amm12592_r128_repair2_boxeph as rp2
    for si, (acc, sig) in enumerate(states):
        for sgn in (1, -1):
            r = rp2.repair2(sig, d[R - 2], d[R - 1], sgn)
            if r is not None:
                log(f"  repair2 endgame: state {si} sign {sgn}")
                return acc + [r[0], r[1]], f"REPAIR2@state{si}sign{sgn}"
    return None, "exhausted"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beam", type=int, default=400)
    ap.add_argument("--ctrl", type=int, default=2)
    ap.add_argument("--span", type=int, default=2)
    ap.add_argument("--dedup", type=int, default=999)
    ap.add_argument("--R", type=int, default=128)
    ap.add_argument("--rank", default="l1deg", choices=["l1deg", "hookdist"])
    ap.add_argument("--out", default=WITNESS)
    args = ap.parse_args()
    R = args.R
    d = eng.prof(R)
    caps = eng.coeff_caps(R, d)
    log = lambda s: print(s, flush=True)
    log(f"HOOK SOLVER R={R} beam={args.beam} ctrl={args.ctrl} "
        f"span={args.span} rank={args.rank} "
        f"(exact caps prune + attractor hook)")
    sol, msg = solve_hook(R, d, caps, args.beam, args.ctrl, args.span,
                          args.dedup, log, CKPT,
                          f"hook-{args.rank}-beam{args.beam}",
                          rank=args.rank)
    log(f"  -> {msg}")
    if sol is None:
        return 1
    a, b = eng.verify_witness(R, sol, d)
    log(f"  EXACT VERIFY: admissible={a} identity={b}")
    assert a and b, "hook witness failed exact verification -- NOT saved"
    eng.write_json(args.out, {
        "R": R, "profile": d, "blocks": sol, "verified": True, "H": [1],
        "source_label": "gamma* floor (direct l1deg beam + attractor hook)",
        "search": {"beam": args.beam, "ctrl": args.ctrl, "span": args.span,
                   "rank": "l1deg", "dedup": args.dedup,
                   "prune": "caps+eval1", "mode": "attractor-hook",
                   "provenance": msg}})
    log(f"  WITNESS WRITTEN: {args.out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
