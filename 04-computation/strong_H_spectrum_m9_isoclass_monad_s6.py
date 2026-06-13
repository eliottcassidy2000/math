"""strong-min(9) via iso-class generation (monad-compute-2026-06-06-S6).

Follow-up to S5, which exhaustively established strong-min(m)=3,5,9,15,25,45 for m=3..8
(generator validated against A000568, strong spectra matched HYP-2180 canon). The Busch
recurrence guess p(n)=p(n-1)+p(n-2)+1 (=>41 at n=8) is REFUTED by S5 (actual 45). This
script computes strong-min(9) and the m=9 strong H-spectrum to characterize the true law,
by generating all A000568(9)=191536 non-isomorphic tournaments on 9 vertices.

Optimized canonical form: colour refinement (1-WL) with discreteness early-exit; most
tournaments discretize to singletons => one encode. Reuses the validated S5 building blocks.
"""
import sys, time
from itertools import permutations, product

sys.stdout.reconfigure(line_buffering=True)

def Hcount(n, adj):
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            av = adj[v] & ~mask
            while av:
                w = (av & -av).bit_length() - 1
                dp[mask | (1 << w)][w] += c
                av &= av - 1
    full = size - 1
    return sum(dp[full][v] for v in range(n))

def is_strong(n, adj):
    FULL = (1 << n) - 1
    seen = 1; frontier = 1
    while frontier:
        nf = 0; mm = frontier
        while mm:
            v = (mm & -mm).bit_length() - 1; mm &= mm - 1
            nf |= adj[v]
        nf &= ~seen
        if not nf:
            break
        seen |= nf; frontier = nf
    if seen != FULL:
        return False
    radj = [0] * n
    for v in range(n):
        av = adj[v]
        while av:
            w = (av & -av).bit_length() - 1
            radj[w] |= 1 << v
            av &= av - 1
    seen = 1; frontier = 1
    while frontier:
        nf = 0; mm = frontier
        while mm:
            v = (mm & -mm).bit_length() - 1; mm &= mm - 1
            nf |= radj[v]
        nf &= ~seen
        if not nf:
            break
        seen |= nf; frontier = nf
    return seen == FULL

def _encode(n, adj, perm):
    code = 0; b = 0
    for p in range(n):
        ap = adj[perm[p]]
        for q in range(p + 1, n):
            code |= ((ap >> perm[q]) & 1) << b
            b += 1
    return code

def canon(n, adj):
    inmask = [0] * n
    for v in range(n):
        av = adj[v]
        while av:
            w = (av & -av).bit_length() - 1
            inmask[w] |= 1 << v
            av &= av - 1
    colour = [bin(adj[v]).count("1") for v in range(n)]
    ndist = len(set(colour))
    while ndist < n:
        sig = []
        for v in range(n):
            oc = []; ic = []
            av = adj[v]
            while av:
                w = (av & -av).bit_length() - 1; oc.append(colour[w]); av &= av - 1
            iv = inmask[v]
            while iv:
                w = (iv & -iv).bit_length() - 1; ic.append(colour[w]); iv &= iv - 1
            oc.sort(); ic.sort()
            sig.append((colour[v], tuple(oc), tuple(ic)))
        rank = {s: i for i, s in enumerate(sorted(set(sig)))}
        newc = [rank[sig[v]] for v in range(n)]
        nd = len(rank)
        if newc == colour or nd == ndist:   # stabilized (discrete or not)
            colour = newc
            ndist = nd
            break
        colour = newc; ndist = nd
    cells_by_colour = {}
    for v in range(n):
        cells_by_colour.setdefault(colour[v], []).append(v)
    cells = [cells_by_colour[c] for c in sorted(cells_by_colour)]
    if all(len(c) == 1 for c in cells):           # discrete: single perm
        perm = [c[0] for c in cells]
        return _encode(n, adj, perm)
    best = None
    for combo in product(*[list(permutations(c)) for c in cells]):
        perm = [v for grp in combo for v in grp]
        code = _encode(n, adj, perm)
        if best is None or code < best:
            best = code
    return best

def extend(adj_prev, nprev, ext):
    adj = list(adj_prev) + [0]
    new = nprev
    for i in range(nprev):
        if (ext >> i) & 1:
            adj[new] |= 1 << i
        else:
            adj[i] |= 1 << new
    return adj

def gen_level(prev_reps, nprev):
    """All non-iso tournaments on nprev+1 vertices from reps on nprev."""
    n = nprev + 1
    seen = {}
    for R in prev_reps:
        for ext in range(1 << nprev):
            adj = extend(R, nprev, ext)
            c = canon(n, adj)
            if c not in seen:
                seen[c] = adj
    return list(seen.values())

def decode(code, n):
    adj = [0] * n; b = 0
    for p in range(n):
        for q in range(p + 1, n):
            if (code >> b) & 1:
                adj[p] |= 1 << q
            else:
                adj[q] |= 1 << p
            b += 1
    return adj

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=9)
    ap.add_argument("--probe", action="store_true", help="time partial n=9 gen then exit")
    args = ap.parse_args()

    A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536}
    print("=" * 78)
    print(f"strong-min via iso-class generation up to n={args.nmax}  (monad-compute-S6)")
    print("=" * 78)

    reps = [[0]]
    nprev = 1
    t0 = time.time()
    while nprev < 8:
        reps = gen_level(reps, nprev); nprev += 1
        print(f"  n={nprev}: {len(reps)} classes (A000568={A000568[nprev]}, "
              f"ok={len(reps)==A000568[nprev]})  t={time.time()-t0:.1f}s")
    if args.probe:
        # time canon throughput on first 400 reps' worth of extensions
        tp = time.time(); cnt = 0; seen = set()
        for R in reps[:400]:
            for ext in range(1 << 8):
                canon(9, extend(R, 8, ext)); cnt += 1
        rate = cnt / (time.time() - tp)
        total = 6880 * 256
        print(f"\n  PROBE: {cnt} canon@n=9 in {time.time()-tp:.1f}s => {rate:.0f}/s")
        print(f"  projected full n=9 generation: {total/rate:.0f}s for {total} extensions")
        return

    if args.nmax >= 9:
        reps = gen_level(reps, 8); nprev = 9
        print(f"  n=9: {len(reps)} classes (A000568={A000568[9]}, "
              f"ok={len(reps)==A000568[9]})  t={time.time()-t0:.1f}s")

    m = nprev
    ts = time.time(); sv = {}; nstrong = 0
    for adj in reps:
        if not is_strong(m, adj):
            continue
        nstrong += 1
        h = Hcount(m, adj)
        sv[h] = sv.get(h, 0) + 1
    print("\n" + "-" * 78)
    print(f"  m={m}: #strong-classes={nstrong}/{len(reps)}  strong-min={min(sv)}  "
          f"|spectrum|={len(sv)}  dt={time.time()-ts:.1f}s")
    print(f"  strong-min(m) sequence so far: 3,5,9,15,25,45 (m=3..8, from S5); "
          f"strong-min({m})={min(sv)}")
    full = sorted(sv)
    print(f"  smallest 25 strong H-values at m={m}: {full[:25]}")
    for q in (7, 21, 35, 49, 63, 189):
        print(f"    is {q:3d} a strong H-value at m={m}? {q in sv}")
    # recurrence checks against 3,5,9,15,25,45,X
    seq = [3, 5, 9, 15, 25, 45, min(sv)]
    print(f"\n  candidate sequence m=3..9: {seq}")
    print(f"    p(n)=p(n-1)+p(n-2)+1 predicts: {[3,5,9,15,25,41,67][:7]}  (BROKE at n=8)")
    a, b = seq[-2], seq[-1]
    print(f"    last ratio {b}/{a} = {b/a:.4f}")
    print("DONE.")

if __name__ == "__main__":
    main()
