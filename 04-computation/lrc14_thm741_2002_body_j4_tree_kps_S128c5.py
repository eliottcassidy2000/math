#!/usr/bin/env python3
"""
lrc14_thm741_2002_body_j4_tree_kps_S128c5.py
============================================
kind-pasteur-2026-07-13-S128 (cont.5).  THM-741: the j=4 rung (overnight run).

THEOREM (target): every 13-speed family with AT LEAST 9 speeds in {1..14} satisfies LRC(14).
For every 9-element body E subset {1..14} (all C(14,9)=2002) and all v1<v2<v3<v4 not in E,
{E,v1,v2,v3,v4} is lonely.

Tree per body (THM-735 Bonferroni legs; sqrt2 -> 99/70 safe; MISTAKE-141 exact thresholds):
  J4: all four >= V1(E)          [4/V1 < 3 m_E/(s2 r_E)]
  J3: per-v1 exact body E1; all three >= V2(E1)   [3/V2 < 4 m1/(s2 r1)]
  J2: per-(v1,v2) exact E2; both >= V3(E2)        [2/V3 < 5 m2/(s2 r2)]
  J1: per-(v1,v2,v3): v4 > v0                     [tail, THM-732(iii)]
  BOTTOM: v4 <= v0-bound: covering iff v4 covers Qb = the q in Q(E) not covered by v1,v2,v3:
          enumerate multiples of lcm(Qb) (or ALL v4 when Qb empty); exact-Q sweep;
          non-covering v4 -> THM-366.

KEY OPTIMIZATION (sound; P1: r' <= v m + (15/7) r, P2: m' >= (6/7) m - 8r/(49v), one level off
exact E2): v0(E3) <= v0u := s2*r3u/(6*m3l) with r3u, m3l the P1/P2 bounds, so "v4 > v0u => tail
fires" and the bottom candidate set (covering v4 in (v3, v0u]) is computable WITHOUT the exact E3
body; the exact level-3 subtract happens ONLY at nodes with candidates (sweeping covering v4 up
to v0u >= v0 is redundant-but-sound). Removes the E3-subtract bottleneck (first smoke test:
758,921 E3 subtracts / 278.6 s on one mid body; the naive per-level "lemma-skips" fired 0 times
-- algebraically vacuous inside the loop ranges -- and were removed).

REGRESSION: subtree (E={1..9}, v1=10) must reproduce THM-738 body-{1..10}: V=154, 143 E2-bodies,
7537 v3-nodes, and >= 27 sweeps (the v0-upper bound sweeps a superset of THM-738's 27).

Modes:
  --probe            run PROBE_BODIES serially with timing + regression, no state written
  (default)          full run: multiprocessing (heavy-first), resume-safe JSONL, summary at end

State: JSONL at SCRATCH/thm741_results.jsonl (one line per completed body; restart skips those).
On completion: writes 05-knowledge/results/lrc14_thm741_2002_body_j4_tree_kps_S128c5.out and
copies the JSONL beside it. The job never runs git.
"""
import sys, os, json, time
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd

ONE = F(1)
S2 = F(99, 70)                      # > sqrt2
SCRATCH = r"C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-math\f631d0eb-9f23-494b-bb86-e0501bc456e9\scratchpad"
STATE = os.path.join(SCRATCH, "thm741_results.jsonl")
REPO = r"C:\Users\Eliott\Documents\GitHub\math"
OUT = os.path.join(REPO, "05-knowledge", "results", "lrc14_thm741_2002_body_j4_tree_kps_S128c5.out")

def _build_bad(w):
    r = F(1, 14 * w)
    out = []
    for k in range(w):
        c = F(k, w)
        lo, hi = c - r, c + r
        if lo < 0:
            out.append((F(0), hi)); out.append((lo + 1, ONE))
        elif hi > 1:
            out.append((lo, ONE)); out.append((F(0), hi - 1))
        else:
            out.append((lo, hi))
    return tuple(sorted(out))

_bad_cached = lru_cache(maxsize=2048)(_build_bad)

def bad_pieces(w):
    # two-tier: cache small w only (bounds per-worker memory to ~tens of MB across 12 workers);
    # large w (sparse lcm-multiple candidates / capped flood ranges) built transiently.
    return _bad_cached(w) if w <= 600 else _build_bad(w)

def good_norm(speeds):
    pieces = []
    for w in speeds:
        pieces += bad_pieces(w)
    pieces = sorted(pieces)
    comps = []
    for lo, hi in pieces:
        if comps and lo <= comps[-1][1]:
            if hi > comps[-1][1]:
                comps[-1][1] = hi
        else:
            comps.append([lo, hi])
    out = []
    for i in range(len(comps)):
        a = comps[i][1]
        j = (i + 1) % len(comps)
        b = comps[j][0] + (ONE if j == 0 else 0)
        if a < b:
            if b <= 1:
                out.append((a, b))
            else:
                out.append((a, ONE))
                if b - 1 > 0:
                    out.append((F(0), b - 1))
    out.sort()
    return out, len(out), sum(b - a for a, b in out)

def subtract_sparse(G, w):
    """Measure of G minus D_w, building ONLY the D_w arcs that meet G (cheap for large w).
    Extended-k trick: k in [floor(w*(a-rad))-1, ceil(w*(b+rad))+1] without mod — wrapped arcs'
    unwrapped copies (k=0 at center 0, k=w at center 1) both appear, exactly covering [0,1)."""
    rad = F(1, 14 * w)
    tot = F(0)
    for a, b in G:
        flo = w * (a - rad)
        fhi = w * (b + rad)
        klo = flo.numerator // flo.denominator - 1
        khi = -((-fhi.numerator) // fhi.denominator) + 1
        cur = a
        seg = F(0)
        for k in range(klo, khi + 1):
            c = F(k, w)
            lo, hi = c - rad, c + rad
            if hi <= cur or lo >= b:
                continue
            if lo > cur:
                seg += min(lo, b) - cur
            cur = max(cur, hi)
            if cur >= b:
                break
        if cur < b:
            seg += b - cur
        tot += seg
    return tot

def subtract(G, w):
    B = bad_pieces(w)
    out = []
    j = 0
    for a, b in G:
        cur = a
        while j > 0 and B[j - 1][1] > cur:
            j -= 1
        k = j
        while k < len(B) and B[k][0] < b:
            lo, hi = B[k]
            if hi <= cur:
                k += 1; continue
            if lo > cur:
                out.append((cur, min(lo, b)))
            cur = max(cur, hi)
            if cur >= b:
                break
            k += 1
        if cur < b:
            out.append((cur, b))
        j = k
    return len(out), sum(b - a for a, b in out), out

def lcm(xs):
    L = 1
    for x in xs:
        L = L * x // gcd(L, x)
    return L

def minV(k, num, den):
    """minimal integer V with k/V < num/den."""
    V = (k * den) // num + 1
    while F(k, V) >= F(num, den):
        V += 1
    return V

def body_work(E):
    t0 = time.time()
    Eset = set(E)
    QE = [q for q in range(2, 15) if not any(w % q == 0 for w in E)]
    GE, rE, mE = good_norm(E)
    assert mE > 0
    thr = 3 * mE / (S2 * rE)
    V1 = minV(4, thr.numerator, thr.denominator)
    nE1 = nE2 = nv3 = nE3 = nsw = 0
    tights = []
    fails = []
    reg = None
    for v1 in range(1, V1):
        if v1 in Eset:
            continue
        r1_, m1, G1 = subtract(GE, v1)
        assert m1 > 0
        nE1 += 1
        thr1 = 4 * m1 / (S2 * r1_)
        V2 = minV(3, thr1.numerator, thr1.denominator)
        sub_nE2 = sub_nv3 = sub_nsw = 0
        for v2 in range(v1 + 1, V2):
            if v2 in Eset:
                continue
            r2_, m2, G2 = subtract(G1, v2)
            assert m2 > 0
            nE2 += 1; sub_nE2 += 1
            thr2 = 5 * m2 / (S2 * r2_)
            V3 = minV(2, thr2.numerator, thr2.denominator)
            for v3 in range(v2 + 1, V3):
                if v3 in Eset:
                    continue
                nv3 += 1; sub_nv3 += 1
                # v0 UPPER bound off exact E2 (P1/P2, screen only): r3<=r3u, m3>=m3l>0 => v0<=v0u.
                # Nodes surviving the screen compute the exact E3 and REBUILD the candidate set
                # with the TRUE v0 (~10x tighter -- the probe showed v0u-wide sweeping is 15x
                # slower); sweeps use subtract_sparse (only the D_v4 arcs meeting G3, cheap for
                # large v4). The v0u>2000 guard also avoids the giant-arc-tuple GC crash.
                r3u = v3 * m2 + F(15, 7) * r2_
                m3l = F(6, 7) * m2 - F(8 * r2_, 49 * v3)
                G3 = None
                v0 = S2 * r3u / (6 * m3l) if m3l > 0 else None
                if v0 is None or v0 > 2000:
                    r3_, m3, G3 = subtract(G2, v3)
                    assert m3 > 0
                    nE3 += 1
                    v0 = S2 * F(r3_) / (6 * m3)
                Qb = [q for q in QE if v1 % q != 0 and v2 % q != 0 and v3 % q != 0]
                L0 = lcm(Qb) if Qb else 1
                for _pass in (0, 1):
                    bmax = v0.numerator // v0.denominator + 1
                    while F(bmax) > v0:
                        bmax -= 1
                    if bmax <= v3:
                        cands = []
                        break
                    if Qb:
                        bs = range(((v3 // L0) + 1) * L0, bmax + 1, L0)
                    else:
                        bs = range(v3 + 1, bmax + 1)
                    cands = [b for b in bs if b not in Eset]
                    if not cands or G3 is not None:
                        break
                    # screen passed on v0u: now pay for the exact E3 and shrink to the true v0
                    r3_, m3, G3 = subtract(G2, v3)
                    assert m3 > 0
                    nE3 += 1
                    v0 = S2 * F(r3_) / (6 * m3)
                if not cands:
                    continue
                for v4 in cands:
                    nsw += 1; sub_nsw += 1
                    L = subtract_sparse(G3, v4)
                    if L == 0:
                        fam = tuple(sorted(E + (v1, v2, v3, v4)))
                        mq = [q for q in range(2, 15) if not any(w % q == 0 for w in fam)]
                        tights.append((fam, mq))
                        if not mq:
                            fails.append(fam)
        if E == tuple(range(1, 10)) and v1 == 10:
            reg = (V2, sub_nE2, sub_nv3, sub_nsw)
    return dict(E=list(E), Q=QE, r=rE, m=str(mE), V1=V1,
                E1=nE1, E2=nE2, v3nodes=nv3, E3=nE3, sweeps=nsw,
                tights=[[list(f), q] for f, q in tights], fails=[list(f) for f in fails],
                reg=reg, secs=round(time.time() - t0, 2))

PROBE_BODIES = [
    (2, 3, 4, 5, 6, 8, 9, 10, 12),            # smoke body (278.6s original / 4123s v0u-wide / target: ~2-4 min)
    (1, 2, 8, 9, 10, 11, 12, 13, 14),         # TRUE flood body (Q empty) -- the makespan anchor
]

def main():
    if "--probe" in sys.argv:
        print("=" * 100, flush=True)
        print("THM-741 PROBE (serial, %d bodies)" % len(PROBE_BODIES), flush=True)
        print("=" * 100, flush=True)
        tot = 0.0
        for E in PROBE_BODIES:
            res = body_work(E)
            tot += res["secs"]
            print("  %-30s Q=%-12s V1=%4d exact=%d/%d/%d v3nodes=%d sweeps=%7d fails=%d  [%.1fs]"
                  % ("{" + ",".join(map(str, E)) + "}", res["Q"], res["V1"],
                     res["E1"], res["E2"], res["E3"], res["v3nodes"], res["sweeps"],
                     len(res["fails"]), res["secs"]), flush=True)
            if res["reg"] is not None:
                rg = tuple(res["reg"])
                ok = rg[:3] == (154, 143, 7537) and rg[3] >= 27
                print("    REGRESSION subtree (E={1..9},v1=10): %s  expect (154,143,7537,>=27) -> %s"
                      % (rg, "OK" if ok else "**MISMATCH**"), flush=True)
                assert ok
        print("probe total %.1fs ; naive extrapolation: 2002 bodies ~ %.1f h serial"
              % (tot, 2002 * (tot / len(PROBE_BODIES)) / 3600), flush=True)
        return

    import multiprocessing as mp
    # keep the machine awake for the duration of the run (process-scoped, auto-released on exit;
    # the first overnight attempt died to system sleep with zero completed bodies)
    try:
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)  # ES_CONTINUOUS|ES_SYSTEM_REQUIRED
        print("sleep inhibitor armed (ES_CONTINUOUS|ES_SYSTEM_REQUIRED)", flush=True)
    except Exception as e:
        print("sleep inhibitor unavailable: %r" % (e,), flush=True)
    done = set()
    if os.path.exists(STATE):
        with open(STATE, encoding="utf-8") as f:
            for line in f:
                try:
                    done.add(tuple(json.loads(line)["E"]))
                except Exception:
                    pass
    bodies = [E for E in combinations(range(1, 15), 9) if E not in done]
    # HEAVY-FIRST (minimize makespan): flood bodies (Q empty) first, then descending V1.
    order = []
    for E in bodies:
        QE = [q for q in range(2, 15) if not any(w % q == 0 for w in E)]
        _, rE, mE = good_norm(E)
        thr = 3 * mE / (S2 * rE)
        V1 = minV(4, thr.numerator, thr.denominator)
        order.append((0 if not QE else 1, -V1, E))
    order.sort()
    bodies = [E for _, _, E in order]
    workers = max(2, (os.cpu_count() or 8) - 2)
    print("THM-741 full run: %d bodies remaining (of 2002), %d workers, heavy-first, state=%s"
          % (len(bodies), workers, STATE), flush=True)
    t0 = time.time()
    ndone = len(done)
    with mp.Pool(workers) as pool, open(STATE, "a", encoding="utf-8") as st:
        for res in pool.imap_unordered(body_work, bodies, chunksize=1):
            st.write(json.dumps(res) + "\n")
            st.flush()
            ndone += 1
            print("[%5d/2002] %-30s V1=%4d exact=%d/%d/%d v3n=%d sweeps=%7d fails=%d %.1fs (elapsed %.0f min)"
                  % (ndone, "{" + ",".join(map(str, res["E"])) + "}", res["V1"],
                     res["E1"], res["E2"], res["E3"], res["v3nodes"], res["sweeps"],
                     len(res["fails"]), res["secs"], (time.time() - t0) / 60), flush=True)
    rows = [json.loads(l) for l in open(STATE, encoding="utf-8")]
    assert len({tuple(r["E"]) for r in rows}) == 2002, "incomplete state -- rerun to resume"
    fails = [f for r in rows for f in r["fails"]]
    tights = {tuple(f): q for r in rows for f, q in r["tights"]}
    regs = [tuple(r["reg"]) for r in rows if r.get("reg")]
    reg_ok = len(regs) == 1 and regs[0][:3] == (154, 143, 7537) and regs[0][3] >= 27
    tot_sw = sum(r["sweeps"] for r in rows)
    tot_v3 = sum(r["v3nodes"] for r in rows)
    tot_E3 = sum(r["E3"] for r in rows)
    maxV1 = max((r["V1"], tuple(r["E"])) for r in rows)
    worst = max((r["secs"], tuple(r["E"])) for r in rows)
    lines = []
    lines.append("THM-741 AGGREGATE: 2002 bodies; v3-nodes %d ; exact E3 bodies %d ; bottom exact sweeps %d"
                 % (tot_v3, tot_E3, tot_sw))
    lines.append("max V1 = %d at %s ; worst body %.1fs at %s" % (maxV1[0], maxV1[1], worst[0], worst[1]))
    lines.append("REGRESSION (E={1..9},v1=10): %s  expect (154,143,7537,>=27) -> %s"
                 % (regs, "OK" if reg_ok else "MISMATCH"))
    lines.append("TIGHTS among swept covering quadruples: %s"
                 % (", ".join("%s missing q=%s" % (f, q) for f, q in tights.items()) if tights else "NONE"))
    lines.append("COVERING L=0 (blockers): %s" % (fails if fails else "NONE"))
    if not fails and reg_ok:
        lines.append("THM-741 ESTABLISHED: every 13-speed family with >=9 speeds in {1..14} satisfies LRC(14).")
    else:
        lines.append("THM-741 NOT ESTABLISHED -- investigate (fails or regression).")
    report = "\n".join(lines)
    print(report, flush=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(report + "\n")
    import shutil
    shutil.copy(STATE, OUT.replace(".out", ".results.jsonl"))
    print("summary written to %s ; DONE." % OUT, flush=True)

if __name__ == "__main__":
    main()
