"""AMM 12592 R=128 exact multi-row (L-row) banded attractor repair.

State sigma_t at row t (from a checkpoint), rows t+1..127 to fill, m = 126-t.
Ansatz: rows t+1..126 = ride(sgn at d_r) + even low-cell corrections a^{(r)}_j;
row 127 = exact remainder decode (must be box-admissible).

Exact facts used:
  * pure rides telescope: sum_{r=t+1}^{126} p^{r-t-1} ride_r
      = sgn E_m + sgn p^m,  so  Delta_127 = J/p^m - sgn,  J = sigma_t - sgn E_m,
    requiring p^m | J (state head must match the attractor through m coeffs).
  * a correction a at row r cell j shifts the final decode cells by
      a * C(G_r, k-j) at cell k,  G_r = (d_127 - d_r) + (127 - r)
    (cascade of the per-step degree-gap bands (x+q)^{g}).
  * margins: |ride_j + a_j| <= C(d_r, j), a_j even (parity preserved);
    final |de_k - T_k| <= C(152, k) (parity automatic).

Search: per final cell k = 0..K, the reachable T_k from the new unknowns
{a^{(r)}_k}_r is an interval; intersect with the needed window and distribute
greedily (several row orderings tried).  Assembly is EXACT (no band formula
on the trusted path): remainder polynomial division + decode + independent
witness verification.

Usage: python3 amm12592_r128_repairL_boxeph.py <ckpt.json> [K] [J]
"""
import json, os, sys, time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng
from amm12592_r128_repair2_boxeph import decode_cells, box_ok, poly_of_cells

def E(m):
    return [-1] + [1] * m

def repairL(sigma, t, d, sgn, K=10, J=5, verbose=False):
    """Returns full rows t+1..127 (list of cell lists) or None."""
    R = 128
    m = 126 - t
    rows = list(range(t + 1, 127))          # correction rows
    d127 = d[127]
    # J-poly and divisibility
    Em = E(m)
    n = max(len(sigma), m + 1)
    Jp = [(sigma[i] if i < len(sigma) else 0) -
          (sgn * Em[i] if i < len(Em) else 0) for i in range(n)]
    if any(Jp[i] != 0 for i in range(min(m, len(Jp)))):
        return None
    Z = Jp[m:] if len(Jp) > m else [0]
    if not Z:
        Z = [0]
    Z[0] -= sgn
    Z = eng.trim(Z) or [0]
    de = decode_cells(Z, d127)
    if de is None:
        return None
    # parity screen on the small cells
    for k in range(min(K + 8, d127 + 1)):
        if (de[k] - comb(d127, k)) % 2:
            return None
    rides = {}
    for r in rows:
        rd = eng.final_decode([-sgn, 2 * sgn], d[r])
        if rd is None:
            return None
        rides[r] = rd
    G = {r: (d127 - d[r]) + (127 - r) for r in rows}
    c2 = [comb(d127, k) for k in range(K + 2)]
    # margins: a^{(r)}_j in [loM, hiM] even
    marg = {}
    for r in rows:
        for j in range(J + 1):
            c = comb(d[r], j)
            rk = rides[r][j]
            lo, hi = -c - rk, c - rk
            lo += lo % 2
            hi -= hi % 2
            marg[(r, j)] = (lo, hi)

    orderings = [sorted(rows, key=lambda r: -G[r]), sorted(rows),
                 sorted(rows, reverse=True)]
    budget = [200000]

    def distribute(order, k, rem):
        """distribute an even total `rem` over rows' cell-k margins; returns
        dict r->v or None."""
        out = {}
        for idx, r in enumerate(order):
            lo_r, hi_r = marg[(r, k)]
            rest_lo = sum(marg[(rr, k)][0] for rr in order[idx + 1:])
            rest_hi = sum(marg[(rr, k)][1] for rr in order[idx + 1:])
            v = max(lo_r, min(hi_r, rem - rest_lo))
            v = max(v, rem - rest_hi)
            v = max(lo_r, min(hi_r, v))
            v -= v % 2
            out[r] = v
            rem -= v
        return out if rem == 0 else None

    def dfs(order, k, carry, a):
        if budget[0] <= 0:
            return None
        budget[0] -= 1
        if k > K:
            return a
        lo_need = de[k] - c2[k] - carry[k]
        hi_need = de[k] + c2[k] - carry[k]
        if k <= J:
            lo_sum = sum(marg[(r, k)][0] for r in rows)
            hi_sum = sum(marg[(r, k)][1] for r in rows)
        else:
            lo_sum = hi_sum = 0
        lo_t = max(lo_need, lo_sum)
        hi_t = min(hi_need, hi_sum)
        lo_t += lo_t % 2
        hi_t -= hi_t % 2
        if lo_t > hi_t:
            return None
        # candidate totals: exact-correct (near de_k - carry), extremes, zero
        want = de[k] - carry[k]
        want -= want % 2
        cands = []
        for T in (max(lo_t, min(hi_t, want)), lo_t, hi_t,
                  0 if lo_t <= 0 <= hi_t else None,
                  (lo_t + hi_t) // 2 - ((lo_t + hi_t) // 2) % 2):
            if T is None or T < lo_t or T > hi_t or T in cands:
                continue
            cands.append(T)
        for T in cands[:4]:
            dist = distribute(order, k, T) if k <= J else ({} if T == 0
                                                           else None)
            if dist is None:
                continue
            nc = list(carry)
            for r, v in dist.items():
                if v:
                    for i in range(1, min(G[r], K + J + 7 - k) + 1):
                        nc[k + i] += v * comb(G[r], i)
            na = dict(a)
            for r, v in dist.items():
                na[(r, k)] = v
            res = dfs(order, k + 1, nc, na)
            if res is not None:
                return res
        return None

    for order in orderings:
        budget[0] = 200000
        a = dfs(order, 0, [0] * (K + J + 8), {})
        if a is None:
            continue
        # EXACT assembly and verification
        blocks = []
        acc_poly = None
        for r in rows:
            cells = list(rides[r])
            for j in range(J + 1):
                v = a.get((r, j), 0)
                if v:
                    cells[j] = rides[r][j] + v
            blocks.append(cells)
        # remainder: sigma_t - sum p^{r-t-1} poly(block_r), then / p^{126-t}
        tot = [0] * (200 + d127)
        for i, c in enumerate(sigma):
            tot[i] += c
        for ri, r in enumerate(rows):
            poly = poly_of_cells(blocks[ri], d[r])
            off = r - t - 1
            for i, c in enumerate(poly):
                tot[off + i] -= c
        sh = 126 - t
        if any(tot[i] != 0 for i in range(sh)):
            continue
        Zx = eng.trim(tot[sh:]) or [0]
        de127 = decode_cells(Zx, d127)
        if de127 is None or not box_ok(de127, d127):
            continue
        return blocks + [de127]
    return None

def main():
    R = 128
    d = eng.prof(R)
    out = os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")
    ck_path = sys.argv[1]
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    Jc = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    with open(ck_path) as f:
        ck = json.load(f)
    t = ck["row"]
    print(f"repairL on {os.path.basename(ck_path)} (row {t}, "
          f"{len(ck['states'])} states, K={K}, J={Jc})", flush=True)
    t0 = time.time()
    diag = []
    for si, st in enumerate(ck["states"]):
        sig = st["sigma"]
        for sgn in (1, -1):
            r = repairL(sig, t, d, sgn, K=K, J=Jc)
            if r is None:
                continue
            blocks = list(st["acc"]) + r
            a, b = eng.verify_witness(R, blocks, d)
            print(f"  repairL state {si} sign {sgn}: verify adm={a} id={b}",
                  flush=True)
            if a and b:
                eng.write_json(out, {
                    "R": R, "profile": d, "blocks": blocks, "verified": True,
                    "H": [1],
                    "source_label": f"gamma* floor (beam ckpt row {t} + "
                                    f"exact {127-t}-row banded repair)",
                    "search": {"mode": "repairL", "row": t, "state": si,
                               "sign": sgn, "K": K, "J": Jc}})
                print(f"  WITNESS WRITTEN: {out}")
                print(f"CLOSED ({time.time()-t0:.1f}s)")
                return 0
    print(f"repairL: no closure over {len(ck['states'])} states "
          f"({time.time()-t0:.1f}s)")
    return 1

if __name__ == "__main__":
    sys.exit(main())
