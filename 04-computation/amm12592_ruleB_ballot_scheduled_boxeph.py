#!/usr/bin/env python3
"""AMM 12592 -- ANGLE B1: ballot-scheduled attractor rule ("rule B"), boxeph.

Fixes rule A's transient by steering the forced unit cell delta_{i,0} = +-1
(the ONLY contributor to Delta_i(1)) so that the ballot walk
    sigma_i(1) = sigma_{i-1}(1) - delta_{i,0},   sigma_{-1}(1) = 0
tracks a corridor-legal schedule instead of drifting.

Schedules (policy):
  none : delta_{i,0} taken from the clamped ideal (EXACTLY rule A; baseline +
         instrumented anatomy).
  zero : walk sigma(1) toward 0 (task's literal suggestion).
  tent : walk sigma(1) toward min(i+1, R-3-i) -- the fastest lawful approach
         to the E_m attractor's ballot value E_{R-2-i}(1) = R-3-i.  In the
         attractor this is a no-op (ideal already picks +1), so tent only
         changes the transient.
Optional --k1: also steer the k=1 cell (the only other cell in Delta'(1) =
  d*delta_0 - delta_1) so tau_i = sigma_i'(1) tracks the attractor value
  E_m'(1) = m(m+1)/2 (m = R-2-i).  In the attractor the needed value is
  exactly ballot(d)[1] = d-2 (checked in selftest), so again a transient-only
  modification.
Clamp modes:
  cellwise : rule A's clamp (exact Bernstein coordinates of the truncated
             ideal, clamped cell-by-cell, parity fixed).
  seq      : triangular lowest-coefficient-first forced decode WITH clamping
             (junk left in place per coefficient instead of spread by the
             cell map).

Ballot cut (hard, necessary): |sigma_i(1)| <= R-1-i, parity automatic.
Both steering policies choose among corridor-legal signs only.

Exact integer arithmetic ONLY (no floats in any decision; the float in
floor_gamma_star is a seed corrected by exact Lucas/Fibonacci tests).
Search negatives are never infeasibility evidence.

Usage:
  --selftest
  --run --R 128 --D0 0 --clamp cellwise --policy tent [--k1] [--pf tozero]
"""
import argparse, json, os, sys, time

WT = "/tmp/math-wt-boxeph-multifront"
COMP = os.path.join(WT, "04-computation")
RES = os.path.join(WT, "05-knowledge", "results")

# ---------------------------------------------------------------- pascal cache
_PC = [[1]]
def crow(n):
    while len(_PC) <= n:
        r = _PC[-1]
        _PC.append([1] + [r[i] + r[i + 1] for i in range(len(r) - 1)] + [1])
    return _PC[n]

def qpow(m):
    row = crow(m)
    return [row[k] if k % 2 == 0 else -row[k] for k in range(m + 1)]

def Em(m):
    return [-1] + [1] * m if m >= 0 else []

def ptrim(a):
    n = len(a)
    while n and a[n - 1] == 0:
        n -= 1
    return a[:n]

def psub(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i, v in enumerate(a):
        r[i] += v
    for i, v in enumerate(b):
        r[i] -= v
    return ptrim(r)

# ------------------------------------------- exact floor (Lucas/Fib, referee)
def _fib_pair(n):
    if n == 0:
        return (0, 1)
    a, b = _fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)

def five_pow_le_phi2m(d, m):
    if d < 0:
        return True
    F, F1 = _fib_pair(2 * m)
    L = 2 * F1 - F
    A = 2 * 5 ** d - L
    if A <= 0:
        return True
    return A * A < 5 * F * F

def floor_gamma_star(m):
    d = int(m * 0.5979874356654402)
    while five_pow_le_phi2m(d + 1, m):
        d += 1
    while not five_pow_le_phi2m(d, m):
        d -= 1
    return d

# ----------------------------------------------------------------- bern algebra
def bern_to_poly(delta, d):
    acc = [0] * (d + 1)
    for k, v in enumerate(delta):
        if v:
            rowk = crow(k)
            off = d - k
            for s in range(k + 1):
                acc[off + s] += v * (rowk[s] if s % 2 == 0 else -rowk[s])
    return ptrim(acc)

def poly_to_bern(P, d):
    cells = [0] * (d + 1)
    for j, c in enumerate(P):
        if c:
            row = crow(d - j)
            for t in range(d - j + 1):
                cells[t] += c * row[t]
    return cells

def _parity_fix(w, v, cap, pf):
    """w already box-clamped; force w == cap (mod 2), stay in box."""
    if (w - cap) % 2 == 0:
        return w
    if pf == 'toward':
        if v > w and w + 1 <= cap:
            return w + 1
        if v < w and w - 1 >= -cap:
            return w - 1
    return w - 1 if w > 0 else w + 1

def clamp_cellwise(T, d, pf):
    cells = poly_to_bern(T, d)
    caps = crow(d)
    out = [0] * (d + 1)
    junk = 0
    for k in range(d + 1):
        v = cells[k]
        cap = caps[k]
        w = v if -cap <= v <= cap else (cap if v > cap else -cap)
        w = _parity_fix(w, v, cap, pf)
        junk += abs(v - w)
        out[k] = w
    return out, junk

def clamp_seq(T, d, pf):
    r = list(T) + [0] * (d + 1 - len(T)) if len(T) < d + 1 else list(T)
    caps = crow(d)
    out = [0] * (d + 1)
    for j in range(d + 1):
        k = d - j
        a = r[j]
        cap = caps[k]
        w = a if -cap <= a <= cap else (cap if a > cap else -cap)
        w = _parity_fix(w, a, cap, pf)
        out[k] = w
        if w:
            rowk = crow(k)
            r[j] -= w
            for s in range(1, k + 1):
                r[j + s] -= w * (rowk[s] if s % 2 == 0 else -rowk[s])
    junk = sum(abs(x) for x in r)
    return out, junk

def admissible(delta, d):
    if len(delta) != d + 1:
        return False
    caps = crow(d)
    return all(abs(c) <= caps[k] and (c - caps[k]) % 2 == 0
               for k, c in enumerate(delta))

# ----------------------------------------------------------------- schedules
def sig1_target(policy, i, R):
    if policy == 'zero':
        return 0
    # tent: fastest lawful approach to the attractor ballot value R-3-i
    if i >= R - 1:
        return 0
    return min(i + 1, R - 3 - i)

def tau_target(policy, i, R):
    if policy == 'zero':
        return 0
    m = R - 2 - i
    return m * (m + 1) // 2 if m >= 0 else 0

# ----------------------------------------------------------------- the rule
def run_rule(R, D0, clamp='cellwise', policy='none', k1=False, pf='tozero',
             log=print, trace_every=32):
    prof = [floor_gamma_star(R + i) + D0 for i in range(R)]
    clampf = clamp_cellwise if clamp == 'cellwise' else clamp_seq
    sig = qpow(R - 1)
    sig1 = 0            # sigma_{i-1}(1)
    tau = 0             # sigma_{i-1}'(1)
    blocks = []
    hist = []           # per-row health record
    t0 = time.time()
    for i in range(R):
        d = prof[i]
        m = R - 2 - i
        # ideal = sig - x * E_m
        ideal = list(sig)
        if m >= 0:
            need = m + 2
            if len(ideal) < need:
                ideal += [0] * (need - len(ideal))
            e = Em(m)
            for j, c in enumerate(e):
                ideal[j + 1] -= c
        T = ideal[:d + 1]
        delta, junk = clampf(T, d, pf)
        override0 = 0
        if policy != 'none':
            budget = R - 1 - i          # |sigma_i(1)| <= budget
            targ = sig1_target(policy, i, R)
            opts = []
            for s in (1, -1):
                v = sig1 - s
                if abs(v) <= budget:
                    opts.append((abs(v - targ), s))
            if opts:
                opts.sort()
                if len(opts) == 2 and opts[0][0] == opts[1][0]:
                    s = delta[0] if delta[0] in (1, -1) else opts[0][1]
                else:
                    s = opts[0][1]
                if delta[0] != s:
                    override0 = 1
                delta[0] = s
        override1 = 0
        if k1 and policy != 'none' and d >= 2:
            s0 = delta[0]
            s1new = sig1 - s0
            need1 = tau_target(policy, i, R) - tau + d * s0 + s1new
            cap1 = d
            w = need1 if -cap1 <= need1 <= cap1 else (cap1 if need1 > cap1 else -cap1)
            if (w - cap1) % 2:
                cands = [c for c in (w - 1, w + 1) if -cap1 <= c <= cap1]
                w = min(cands, key=lambda c: abs(c - need1))
            if delta[1] != w:
                override1 = 1
            delta[1] = w
        D = bern_to_poly(delta, d)
        t = psub(sig, D)
        if t and t[0] != 0:
            rec = dict(die_row=i, const_bits=abs(t[0]).bit_length(),
                       sig1=sig1, targ=sig1_target(policy, i, R) if policy != 'none' else None,
                       head_bits=[abs(c).bit_length() for c in sig[:12]],
                       junk_bits=junk.bit_length(),
                       elapsed=round(time.time() - t0, 1))
            return None, rec, blocks, hist
        sig = t[1:] if t else []
        blocks.append(delta)
        sig1 -= delta[0]
        tau -= (d * delta[0] - (delta[1] if d >= 1 else 0)) + sig1
        # exact invariants (cheap)
        assert sig1 == sum(sig), (i, sig1)
        assert tau == sum(j * c for j, c in enumerate(sig)), (i, tau)
        mb = max((abs(c) for c in sig), default=0).bit_length()
        hb = max((abs(c) for c in sig[:8]), default=0).bit_length()
        hist.append((i, sig1, hb, mb, junk.bit_length(), override0, override1))
        if trace_every and (i % trace_every == 0 or i >= R - 3):
            log(f"  i={i:4d} d={d:4d} sig1={sig1:5d} headbits={hb:4d} "
                f"maxbits={mb:4d} junkbits={junk.bit_length():4d} "
                f"ov0={override0} ov1={override1} t={time.time()-t0:7.1f}s")
    if sig:
        rec = dict(die_row=R, const_bits=-1, residual_L1=sum(abs(c) for c in sig),
                   elapsed=round(time.time() - t0, 1))
        return None, rec, blocks, hist
    # closed: exact verification (admissibility + identity + profile)
    adm = all(admissible(blocks[i], prof[i]) for i in range(R))
    S = []
    for i in range(R):
        pl = bern_to_poly(blocks[i], prof[i])
        L = i + len(pl)
        if len(S) < L:
            S += [0] * (L - len(S))
        for j, c in enumerate(pl):
            S[i + j] += c
    idt = ptrim(S) == qpow(R - 1)
    fl = [floor_gamma_star(R + i) for i in range(R)]
    pok = all(prof[i] == fl[i] + D0 for i in range(R))
    rec = dict(closed=True, adm=adm, identity=idt, profile_ok=pok,
               elapsed=round(time.time() - t0, 1))
    return (blocks if (adm and idt and pok) else None), rec, blocks, hist

# ----------------------------------------------------------------- selftest
def selftest():
    print("== rule B selftest ==")
    # 1. bern round trip + seq==cellwise when no clamping binds
    import random
    rng = random.Random(12592)
    ok = True
    for d in (5, 9, 14):
        caps = crow(d)
        for _ in range(25):
            delta = [rng.choice(range(-caps[k], caps[k] + 1, 2)) for k in range(d + 1)]
            P = bern_to_poly(delta, d)
            P = list(P) + [0] * (d + 1 - len(P))
            ok &= poly_to_bern(P, d) == delta
            cw, j1 = clamp_cellwise(P, d, 'tozero')
            sq, j2 = clamp_seq(P, d, 'tozero')
            ok &= cw == delta == sq and j1 == 0 == j2
    print(f"round trip + exact-decode agreement (cellwise==seq, junk 0): {'PASS' if ok else 'FAIL'}")
    assert ok
    # 2. full box represents -1; ballot cells represent 2x-1
    ok = True
    for d in range(1, 60):
        caps = crow(d)
        ok &= bern_to_poly([-caps[k] for k in range(d + 1)], d) == [-1]
        bal = [crow(d - 1)[k] - (crow(d - 1)[k - 1] if 1 <= k <= d - 1 else 0)
               if k <= d - 1 else -crow(d - 1)[d - 1] for k in range(d + 1)]
        ok &= bern_to_poly(bal, d) == [-1, 2]
    print(f"full box == -1, ballot == 2x-1 (d<=59): {'PASS' if ok else 'FAIL'}")
    assert ok
    # 3. attractor compatibility of the steering targets:
    #    sigma_{i-1} = E_{m+1} => ideal = 2x-1; tent picks delta0 = +1 = ballot[0];
    #    k1 need = d-2 = ballot[1]
    R = 64
    for i in (40, 50, 60):
        m = R - 2 - i
        d = floor_gamma_star(R + i)
        sig_att = Em(m + 1)
        sig1_att = sum(sig_att)
        tau_att = sum(j * c for j, c in enumerate(sig_att))
        targ = sig1_target('tent', i, R)
        # descend branch: sig1_att - 1 == R-3-i == targ
        assert sig1_att - 1 == targ == R - 3 - i, (i, sig1_att, targ)
        need1 = tau_target('tent', i, R) - tau_att + d * 1 + (sig1_att - 1)
        assert need1 == d - 2, (i, need1, d)
    print("attractor compatibility: tent+k1 steering is a no-op on E_m tail: PASS")
    # 4. policy none == rule A (original module) at R=32, D0=0
    sys.path.insert(0, COMP)
    from amm12592_allR_greedy_attractor_rule_boxeph import run_rule as runA
    solA, msgA, _ = runA(32, 'tozero')
    solB, rec, blocks, _h = run_rule(32, 0, 'cellwise', 'none', False, 'tozero',
                                     log=lambda s: None, trace_every=0)
    same = solA is not None and solB is not None and \
        all(list(a) == list(b) for a, b in zip(solA, solB))
    print(f"policy none reproduces rule A exactly at R=32 (blocks identical): "
          f"{'PASS' if same else 'FAIL'} [{msgA} | {rec}]")
    assert same
    # 5. (S3) parity sanity: 2G_R = q^{R-1} - E_{R-1} even coefficientwise iff R dyadic
    def twoG(R):
        return psub(qpow(R - 1), Em(R - 1))
    for R2, dy in ((8, True), (16, True), (32, True), (12, False), (100, False), (256, True)):
        ev = all(c % 2 == 0 for c in twoG(R2))
        assert ev == dy, (R2, ev)
        g = twoG(R2)
        assert g[0] == 2 and all(g[j] == (-1) ** j * crow(R2 - 1)[j] - 1
                                 for j in range(1, R2)), R2
    print("(S3) recomputed: [x^0]2G=2, [x^j]2G=(-1)^j C(R-1,j)-1, even iff R=2^t: PASS")
    print("== selftest PASSED ==")

# ----------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--R", type=int, default=128)
    ap.add_argument("--D0", type=int, default=0)
    ap.add_argument("--clamp", default='cellwise', choices=['cellwise', 'seq'])
    ap.add_argument("--policy", default='tent', choices=['none', 'zero', 'tent'])
    ap.add_argument("--k1", action="store_true")
    ap.add_argument("--pf", default='tozero', choices=['tozero', 'toward'])
    ap.add_argument("--trace-every", type=int, default=32, dest="trace_every")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return 0
    if not args.run:
        return 0
    tag = args.tag or f"{args.clamp}-{args.policy}{'-k1' if args.k1 else ''}" \
                      f"{'-toward' if args.pf == 'toward' else ''}"
    label = f"R={args.R} D0={args.D0} [{tag}]"
    print(f"RULE-B {label} start", flush=True)
    sol, rec, blocks, hist = run_rule(args.R, args.D0, args.clamp, args.policy,
                                      args.k1, args.pf,
                                      log=lambda s: print(s, flush=True),
                                      trace_every=args.trace_every)
    # always dump the trace (die anatomy or closure health)
    trace_path = os.path.join(
        RES, f"amm12592_ruleB_trace_R{args.R}_D0_{args.D0}_{tag}_boxeph.json")
    with open(trace_path, "w") as f:
        json.dump({"R": args.R, "D0": args.D0, "variant": tag, "rec": rec,
                   "hist_cols": ["i", "sig1", "headbits", "maxbits",
                                 "junkbits", "ov0", "ov1"],
                   "hist": hist}, f)
    if sol is None:
        print(f"RULE-B {label} NEGATIVE: {rec}", flush=True)
        print(f"trace: {trace_path}", flush=True)
        return 1
    print(f"RULE-B {label} CLOSED VERIFIED: {rec}", flush=True)
    wpath = os.path.join(COMP, f"amm12592_witness_R{args.R}_ruleB_D0_{args.D0}_boxeph.json")
    prof = [floor_gamma_star(args.R + i) + args.D0 for i in range(args.R)]
    with open(wpath, "w") as f:
        json.dump({"R": args.R, "D0": args.D0, "profile": prof,
                   "floor_profile": [d - args.D0 for d in prof],
                   "blocks": sol, "verified": True,
                   "ballot_trajectory": [h[1] for h in hist],
                   "source": f"rule B ballot-scheduled ({tag}), "
                             f"d_i = exact_floor(gamma*(R+i)) + {args.D0}"}, f)
    print(f"witness: {wpath}", flush=True)
    print(f"trace: {trace_path}", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
