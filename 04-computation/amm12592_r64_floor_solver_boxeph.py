"""AMM 12592 / THM-3029 follow-on: close epoch R = 64 at the gamma* FLOOR
profile d_i = floor(gamma*(64+i)), D0 = 0  (gamma* = log(phi)/log(sqrt 5)).

Machinery (THM-2966 / THM-3002 / THM-3026 / THM-3029):
  Epoch closure at R with degree profile (d_i) is the representation problem
      (*)  q^{R-1} = sum_{i=0}^{R-1} p^i Delta_i(p),   q = 1-p,
  where Delta_i = sum_k delta_{i,k} B_{d_i,k},  B_{d,k}(x) = x^{d-k}(1-x)^k,
  and each block is ADMISSIBLE:  |delta_{i,k}| <= binom(d_i,k)  and
  delta_{i,k} == binom(d_i,k) (mod 2)   (Lucas box capacity + parity).

  Profile monotonicity (THM-3029 (M) / THM-3026 (L)): if ANY pointwise-smaller
  profile closes, the floor profile closes, by convolving each block with the
  constant-1 block [binom(d'-d,k)]_k.

This file is BOTH the search engine and the exact referee:
  * solve_prof(): a generalized-profile beam search in the residual recursion
        sigma_{-1} = q^{R-1},  p sigma_i = sigma_{i-1} - Delta_i,
    steering the first `ctrl` coefficients of each new residual through
    enumerated small targets (the committed amm12592_gamma35_beam_deathstar.py
    policy, generalized to arbitrary profiles + stratified randomized beam +
    endgame widening).  Beam NEGATIVES ARE SEARCH ARTEFACTS (THM-3029 sec 1)
    and are never treated as evidence.
  * verify_witness(): EXACT integer verification, independent of the search:
    every block admissible at its profile degree, and the epoch identity
    sum_i x^i Delta_i == (1-x)^{R-1} as an identity in Z[x].

Exact arithmetic only (int); floats never touch the verification path.

Usage:
  python3 amm12592_r64_floor_solver_boxeph.py            # referee: re-verify + (re)produce witness JSON
  python3 amm12592_r64_floor_solver_boxeph.py --hunt ... # search driver (see --help)
"""
import argparse, json, os, random, sys
from fractions import Fraction
from itertools import product
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))

# gamma* = log(phi)/log(sqrt 5) = 0.59798743566544014974502650...; the 16-digit
# truncation has the same floor profile as gamma* for all m <= 127 (checked in
# recon with sympy at 60 digits; nearest approach of gamma*m to an integer on
# [64,127] is far above the 1e-16 truncation error).
GS = (5979874356654402, 10**16)

def prof(R, g1, g2, D0):
    return [(g1 * (R + i)) // g2 + D0 for i in range(R)]

# ---------------------------------------------------------------------------
# cached (1-x)^k coefficient tails:  B_{d,k} = x^{d-k} (1-x)^k, so the dense
# coefficient vector of B_{d,k} is QK[k] placed at offset d-k.
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

def step(sigma, d, targets):
    """One residual step at degree d.  targets[j] steers sigma_next[j]
    (j = 0..len(targets)-1); None = greedy clamp.  Returns (deltas, sigma_next)
    or None."""
    if not sigma or abs(sigma[0]) != 1:
        return None
    res = list(sigma) + [0] * (d + 2 - len(sigma) if len(sigma) < d + 2 else 1)
    if len(res) < d + 1:
        res += [0] * (d + 1 - len(res))
    deltas = [0] * (d + 1)
    v = sigma[0]
    deltas[d] = v                      # forced: B_{d,d} = (1-x)^d, offset 0
    t = qk(d)
    for j in range(min(len(t), len(res))):
        res[j] -= v * t[j]
    for j, tgt in enumerate(targets):
        k = d - 1 - j
        if k < 0:
            return None
        cap = comb(d, k)
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
            for s in range(k + 1):
                res[off + s] -= want * t[s]
    for k in range(d - 1 - len(targets), -1, -1):
        cap = comb(d, k)
        want = res[d - k]
        vv = max(-cap, min(cap, want))
        if (vv - cap) % 2:
            vv = vv - 1 if vv - 1 >= -cap else vv + 1
        deltas[k] = vv
        if vv:
            t = qk(k)
            off = d - k
            for s in range(k + 1):
                res[off + s] -= vv * t[s]
    if res[0] != 0:
        return None
    return deltas, trim(res[1:])

def final_decode(sigma, d):
    """Exact Bernstein-d decode of sigma; returns deltas iff admissible and
    exactly representable (i.e. Delta = sigma), else None."""
    if len(sigma) - 1 > d:
        return None
    res = list(sigma) + [0] * (d + 1 - len(sigma))
    de = [0] * (d + 1)
    for k in range(d, -1, -1):
        cap = comb(d, k)
        want = res[d - k]
        if abs(want) > cap or (want - cap) % 2:
            return None
        de[k] = want
        if want:
            t = qk(k)
            off = d - k
            for s in range(k + 1):
                res[off + s] -= want * t[s]
    return de if not trim(res) else None

RANKS = {
    "deghead": lambda sg: (len(sg), sum(abs(x) for x in sg[:6])),
    "l1": lambda sg: sum(abs(x) for x in sg),
    "l1deg": lambda sg: (len(sg), sum(abs(x) for x in sg)),
    "l1top": lambda sg: (sum(abs(x) for x in sg[-8:]), len(sg),
                         sum(abs(x) for x in sg)),
}

def solve_prof(d, beam=250, ctrl=2, span=2, seed=None, rand_frac=0.0,
               dedup=10, end_ctrl=None, end_span=None,
               rank="deghead", verbose=False, rhs=None):
    """Beam search for sum_i x^i Delta_i == rhs at profile d (len R), each
    Delta_i admissible at d[i].  rhs defaults to (1-x)^{R-1}, i.e. THM-3002 (*).
    A non-default rhs is the H != 1 widening of THM-3002 sec 1-2: any common
    factor H with rhs = (1-x)^{R-1} H(x) closes the block, and H(x) = H(1-x)
    keeps the mirror ansatz for the 1-side.  Returns (sol, msg).

    Structure: rows 0..R-3 by beam over the residual recursion (rank keys in
    RANKS; 'l1deg' = (degree, full L1) is what closes R = 64 at the floor
    profile -- the winning trajectories shed residual magnitude monotonically,
    which the original (degree, |first 6|) key cannot see).  Then an exhaustive
    WIDE 2-ROW COMPLETION: row R-2 steered over the (end_ctrl, end_span) target
    grid, row R-1 by exact Bernstein decode.  First success wins."""
    R = len(d)
    rankf = RANKS[rank]
    rng = random.Random(seed)
    states = [([], list(rhs) if rhs is not None else qpow(R - 1))]
    for i in range(R - 2):
        opts = [None] + list(range(-span, span + 1))
        nxt = []
        for acc, sig in states:
            for tg in product(opts, repeat=ctrl):
                if tg[0] not in (1, -1):
                    continue
                r = step(sig, d[i], tg)
                if r is None:
                    continue
                de, ns = r
                if not ns or abs(ns[0]) != 1:
                    continue
                nxt.append((acc + [de], ns))
        if not nxt:
            return None, f"died at row {i}"
        nxt.sort(key=lambda st: rankf(st[1]))
        seen, uniq = set(), []
        for a, sg in nxt:
            key = tuple(sg[:dedup])
            if key in seen:
                continue
            seen.add(key)
            uniq.append((a, sg))
        n_top = beam if rand_frac <= 0 else max(1, int(beam * (1 - rand_frac)))
        keep = uniq[:n_top]
        if rand_frac > 0 and len(uniq) > n_top:
            pool = uniq[n_top:]
            extra = rng.sample(pool, min(beam - len(keep), len(pool)))
            keep = keep + extra
        states = keep[:beam]
        if verbose:
            best = states[0][1]
            print(f"  row {i:3d}: states={len(states):5d} "
                  f"best deg={len(best)-1:4d} bestL1={sum(abs(x) for x in best)}",
                  flush=True)
    # wide 2-row completion: row R-2 steered, row R-1 exact decode
    ec = end_ctrl if end_ctrl is not None else ctrl
    es = end_span if end_span is not None else span
    da, db = d[R - 2], d[R - 1]
    eopts = [None] + list(range(-es, es + 1))
    for acc, sig in states:
        for tg in product(eopts, repeat=ec):
            if tg[0] not in (1, -1):
                continue
            r = step(sig, da, tg)
            if r is None:
                continue
            de, ns = r
            if not ns or abs(ns[0]) != 1:
                continue
            fin = final_decode(ns, db)
            if fin is not None:
                return acc + [de, fin], "SOLVED"
    return None, f"beam {beam} exhausted at 2-row completion"

# ---------------------------------------------------------------------------
# exact verification + lifting (THM-3026 (L)), independent of the search
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

def lift_block(delta, dd, dp):
    if dp == dd:
        return list(delta)
    ones = [comb(dp - dd, k) for k in range(dp - dd + 1)]
    r = [0] * (len(delta) + len(ones) - 1)
    for i, u in enumerate(delta):
        if u:
            for j, v in enumerate(ones):
                r[i + j] += u * v
    return r

def verify_witness(R, blocks, profile):
    ok_adm = all(admissible(blocks[i], profile[i]) for i in range(R))
    ok_id = epoch_identity(R, blocks, profile)
    return ok_adm, ok_id

def eff(R, d):
    return max(Fraction(d[i], R + i) for i in range(R))

# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------
WITNESS = os.path.join(HERE, "amm12592_floor_witness_R64.json")

def hunt(args):
    R = args.R
    tgt = prof(R, *GS, 0)
    if args.profile == "floor":
        src = tgt
        label = "gamma* floor"
    else:
        g1, g2, D0 = (int(x) for x in args.profile.split(","))
        src = prof(R, g1, g2, D0)
        label = f"gamma={g1}/{g2}, D0={D0}"
        assert all(src[i] <= tgt[i] for i in range(R)), \
            "source profile not pointwise <= floor profile"
    print(f"HUNT R={R} at {label}: beam={args.beam} ctrl={args.ctrl} "
          f"span={args.span} seed={args.seed} rand_frac={args.rand_frac} "
          f"dedup={args.dedup} rank={args.rank} "
          f"end=({args.end_ctrl},{args.end_span})",
          flush=True)
    sol, msg = solve_prof(src, beam=args.beam, ctrl=args.ctrl, span=args.span,
                          seed=args.seed, rand_frac=args.rand_frac,
                          dedup=args.dedup,
                          end_ctrl=args.end_ctrl, end_span=args.end_span,
                          rank=args.rank, verbose=args.verbose)
    print(f"  -> {msg}", flush=True)
    if sol is None:
        return 1
    a, b = verify_witness(R, sol, src)
    print(f"  source verify: admissible={a} identity={b}", flush=True)
    assert a and b
    lifted = [lift_block(sol[i], src[i], tgt[i]) for i in range(R)]
    a2, b2 = verify_witness(R, lifted, tgt)
    print(f"  lifted to floor profile: admissible={a2} identity={b2}", flush=True)
    assert a2 and b2
    out = args.out or WITNESS
    with open(out, "w") as f:
        json.dump({"R": R, "profile": tgt, "blocks": lifted, "verified": True,
                   "source_label": label, "source_profile": src,
                   "source_blocks": sol,
                   "search": {"beam": args.beam, "ctrl": args.ctrl,
                              "span": args.span, "rank": args.rank,
                              "seed": args.seed,
                              "rand_frac": args.rand_frac, "dedup": args.dedup,
                              "end_ctrl": args.end_ctrl,
                              "end_span": args.end_span}}, f)
    print(f"  WITNESS WRITTEN: {out}", flush=True)
    return 0

WINNING = dict(beam=400, ctrl=2, span=2, seed=None, rand_frac=0.0, dedup=999,
               end_ctrl=3, end_span=6, rank="l1deg")

def reproduce():
    """Deterministic from-scratch reproduction of the R = 64 floor witness
    (rand_frac = 0, so the seed is never consulted; the beam trajectory and
    the completion scan are fully deterministic)."""
    R = 64
    tgt = prof(R, *GS, 0)
    sol, msg = solve_prof(tgt, **WINNING)
    assert sol is not None, f"reproduction failed unexpectedly: {msg}"
    a, b = verify_witness(R, sol, tgt)
    assert a and b, "reproduced witness failed exact verification"
    with open(WITNESS, "w") as f:
        json.dump({"R": R, "profile": tgt, "blocks": sol, "verified": True,
                   "H": [1], "source_label": "gamma* floor (direct, no lift)",
                   "search": dict(WINNING, note="deterministic; seed unused")},
                  f)
    return sol

def referee():
    """Exact re-verification of the saved witness (search-independent).
    Regenerates the witness deterministically if the JSON is absent, so the
    reproduction chain cannot be broken by a lost data file."""
    R = 64
    tgt = prof(R, *GS, 0)
    print("AMM 12592 / R = 64 gamma* floor-profile closure (extends THM-3029 (A))")
    print(f"floor profile d_i = floor(gamma*(64+i)), i=0..63, gamma* = "
          f"{GS[0]}/{GS[1]} (= log(phi)/log(sqrt 5) to 16 digits):")
    print(f"  d = {tgt}")
    if not os.path.exists(WITNESS):
        print("witness JSON absent -- regenerating deterministically "
              f"(beam={WINNING['beam']}, ctrl={WINNING['ctrl']}, "
              f"span={WINNING['span']}, rank={WINNING['rank']}, "
              f"dedup=full, end=({WINNING['end_ctrl']},{WINNING['end_span']}))")
        reproduce()
    with open(WITNESS) as f:
        w = json.load(f)
    assert w["R"] == R and w["profile"] == tgt
    blocks = w["blocks"]
    a, b = verify_witness(R, blocks, tgt)
    print(f"witness blocks: {len(blocks)}; all admissible at profile degree: {a}")
    print(f"epoch identity sum_i x^i Delta_i == (1-x)^63 exact in Z[x]: {b}")
    e = eff(R, tgt)
    print(f"effective rate max_i d_i/(64+i) = {e} = {float(e):.6f} "
          f"(3/5 profile: 0.600000; gamma* = 0.597987)")
    print(f"provenance: {w.get('source_label', 'direct')}; H = {w.get('H', [1])} "
          f"(THM-3002 (*) normal form, no H != 1 widening needed, no lift needed)")
    print(f"sign word delta_(i,0): "
          + "".join("+" if de[0] > 0 else "-" for de in blocks))
    assert a and b
    print()
    print("search provenance (NOT part of the proof; negatives are artefacts,")
    print("THM-3029 sec 1): deterministic beam, rows 0..61 ranked by")
    print("(residual degree, residual L1) -- 'l1deg', beam 400, ctrl 2, span 2,")
    print("full-residual dedup -- then exhaustive 2-row completion, row 62")
    print("steered over the (ctrl 3, span 6) target grid, row 63 by exact")
    print("Bernstein decode.  The L1 rank is decisive: winning trajectories")
    print("shed residual L1 monotonically (2 by row 61), while the")
    print("(degree, |first 6|) key of the committed 3/5 solver keeps states")
    print("whose mid coefficients overshoot final-decode capacity ~1000x.")
    print()
    print("CONCLUSION: the gamma* floor profile CLOSES at R = 64, directly,")
    print("with D0 = 0 and H = 1.  With THM-3029 (A) (R = 8, 16, 32) every")
    print("dyadic epoch through [64,127] closes at the floor profile:")
    print("C = 1 + gamma* = log_5(5 phi^2) = 1.5979874356654402 is ATTAINED")
    print("for all critical values n <= 127, matching the archimedean floor")
    print("(THM-3009/THM-3017/THM-3024) exactly on that range.  R = 128 and")
    print("the all-R statement (which would give C* = log_5(5 phi^2) exactly)")
    print("remain OPEN.")
    print("ALL R64 FLOOR-WITNESS CHECKS PASSED")
    return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hunt", action="store_true")
    ap.add_argument("--R", type=int, default=64)
    ap.add_argument("--profile", default="floor",
                    help="'floor' or 'g1,g2,D0' (pointwise <= floor required)")
    ap.add_argument("--beam", type=int, default=400)
    ap.add_argument("--ctrl", type=int, default=2)
    ap.add_argument("--span", type=int, default=2)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--rand-frac", type=float, default=0.0, dest="rand_frac")
    ap.add_argument("--dedup", type=int, default=999)
    ap.add_argument("--end-ctrl", type=int, default=3, dest="end_ctrl")
    ap.add_argument("--end-span", type=int, default=6, dest="end_span")
    ap.add_argument("--rank", default="l1deg", choices=sorted(RANKS))
    ap.add_argument("--out", default=None)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    if args.hunt:
        sys.exit(hunt(args))
    sys.exit(referee())

if __name__ == "__main__":
    main()
