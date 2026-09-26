#!/usr/bin/env python3
"""Orchestrator audit of the cube-distance lane (HYP-9138), written from the
lane note's statements only (the lane's code was not read).

Parts (run with the part name as argument; default runs A-F):
  A  entropy constants
  B  Bad_k by direct simulation vs a ballot DP; Bad_k in 3 mod 4; -1, -5 in Bad_k
  C  sigma_k = flip exactly Bad_k: class (i) by exact Karp (k<=11) and by an
     integer potential iteration at threshold F_k (k<=16); attainment of F_k
     by the explicit first-descent periodic orbit
  D  integer dynamics of T_(sigma_k): descent within k steps above n_0(k),
     n_0(k) by first-descent-word DFS, every n <= N reaches 1
  E  necklace counts N_k (enumeration and Burnside), N_k >= 2^(hk)/(3k^2) to
     k=200, and the necklace rotations are closed walks of Collatz's parity graph
  F  delta_k for k<=6 by exhaustive search over flip sets of increasing size
  G  delta_7, delta_8 (3n+1) by an implicit hitting set with OR-tools CP-SAT
  H  5n+-1: no class-(i) strategy at k<=5 (exhaustive) and k=6 (IHS, UNSAT);
     a class-(i) strategy exists at k=7; 5n+1's exact distance at k=7
  K  the stored k=10 certificate: every no-good re-derived; own seeds; HiGHS
     proves no hitting set of size <= 39 (CP-SAT returned UNKNOWN after 3000 s)
"""
import sys, math, itertools, gzip, json, time
from fractions import Fraction
import numpy as np

C3 = math.log(2) / math.log(3)


def H2(x):
    return -x * math.log2(x) - (1 - x) * math.log2(1 - x)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


# ---------------------------------------------------------------- graph basics
def succ_table(k, mul, sig):
    """sig[r//2] in {+1,-1} for odd r mod 2^k.  Returns (t, H): node s -> t[s], t[s]+H."""
    M = 1 << k
    H = M >> 1
    t = np.empty(M, dtype=np.int64)
    for s in range(M):
        if s % 2 == 0:
            t[s] = (s // 2) % H
        else:
            t[s] = ((mul * s + sig[s // 2]) // 2) % H
    return t, H


def karp_max_density(k, mul, sig):
    """Exact max odd-density over cycles of G_sigma (Karp, super source)."""
    M = 1 << k
    t, H = succ_table(k, mul, sig)
    src = np.concatenate([np.arange(M), np.arange(M)])
    dst = np.concatenate([t, t + H])
    w = (np.arange(M) % 2).astype(np.int64)
    NEG = -(1 << 40)
    D = np.full((M + 1, M), NEG, dtype=np.int64)
    D[0, :] = 0
    for j in range(1, M + 1):
        prev = D[j - 1]
        vals = prev[src] + w[src]
        vals[prev[src] <= NEG // 2] = NEG
        row = np.full(M, NEG, dtype=np.int64)
        np.maximum.at(row, dst, vals)
        D[j] = row
    best = None
    Dn = D[M]
    js = np.arange(M)
    for v in range(M):
        if Dn[v] <= NEG // 2:
            continue
        col = D[:M, v]
        ok = col > NEG // 2
        ratios = (Dn[v] - col[ok]) / (M - js[ok])
        i = int(np.argmin(ratios))
        j = js[ok][i]
        val = Fraction(int(Dn[v] - col[ok][i]), int(M - j))
        if best is None or val > best:
            best = val
    return best


def potential_ok(k, mul, sig, F):
    """True iff every cycle of G_sigma has odd density <= F (integer potential
    iteration psi = max(0, w_F + max psi(succ)); converges iff no positive cycle)."""
    q, r = F.numerator, F.denominator
    M = 1 << k
    t, H = succ_table(k, mul, sig)
    wF = np.where(np.arange(M) % 2 == 1, r - q, -q).astype(np.int64)
    psi = np.zeros(M, dtype=np.int64)
    for it in range(M + 2):
        new = np.maximum(0, wF + np.maximum(psi[t], psi[t + H]))
        if np.array_equal(new, psi):
            # verify certificate edge by edge
            ok1 = np.all(psi[t] <= psi - wF)
            ok2 = np.all(psi[t + H] <= psi - wF)
            return bool(ok1 and ok2), int(psi.max())
        psi = new
    return False, None


# ---------------------------------------------------------------- Bad_k
def bad_mask(k):
    """Bad_k by simulating Collatz on representatives 0..2^k-1."""
    M = 1 << k
    x = np.arange(M, dtype=object) if k > 40 else np.arange(M, dtype=np.int64)
    a = np.zeros(M, dtype=np.int64)
    bad = np.ones(M, dtype=bool)
    for j in range(1, k + 1):
        odd = (x % 2 == 1)
        a = a + odd
        x = np.where(odd, (3 * x + 1) // 2, x // 2)
        # M_j = 3^a / 2^j > 1  <=>  a > j * log_3 2 (never equal for j >= 1)
        bad &= (a * math.log(3) > j * math.log(2))
    return bad


def ballot_count(k):
    """number of words of length k with 3^(a_j) > 2^j for all j = 1..k."""
    cur = {0: 1}
    for j in range(1, k + 1):
        nxt = {}
        for a, c in cur.items():
            for b in (0, 1):
                a2 = a + b
                if 3 ** a2 > 2 ** j:
                    nxt[a2] = nxt.get(a2, 0) + c
        cur = nxt
    return sum(cur.values())


def sigma_from_mask(k, mask):
    M = 1 << k
    return np.array([-1 if mask[r] else 1 for r in range(1, M, 2)], dtype=np.int64)


def F_k(k):
    best = Fraction(0)
    for d in range(1, k + 1):
        for a in range(0, d + 1):
            if 3 ** a < 2 ** d and Fraction(a, d) > best:
                best = Fraction(a, d)
    return best


# ---------------------------------------------------------------- parts
def part_A():
    print("A. constants")
    h = H2(C3)
    print(f"  log_3 2 = {C3:.10f}, h = H2(log_3 2) = {h:.10f}, 1-h = {1-h:.10f}")
    check(abs(h - 0.9499555) < 5e-7, "h = 0.9499555")
    check(abs((1 - h) - 0.0500445) < 5e-7, "1-h = 0.0500445")


BAD_TABLE = {2: 1, 3: 2, 4: 3, 5: 4, 6: 8, 7: 13, 8: 19, 9: 38, 10: 64, 11: 128,
             12: 226, 13: 367, 14: 734, 16: 2114, 20: 27328}


def part_B():
    print("B. Bad_k")
    for k in range(2, 21):
        m = bad_mask(k)
        nb = int(m.sum())
        bc = ballot_count(k)
        res = np.nonzero(m)[0]
        assert nb == bc, (k, nb, bc)
        assert np.all(res % 4 == 3), k
        M = 1 << k
        assert m[M - 1], ("-1 not bad", k)
        if k >= 3:
            assert m[(M - 5) % M], ("-5 not bad", k)
        assert not m[1]
        if k in BAD_TABLE:
            assert nb == BAD_TABLE[k], (k, nb, BAD_TABLE[k])
    check(True, "|Bad_k| simulation = ballot DP for k=2..20; equals the note's table; Bad_k in 3 mod 4; -1, -5 in Bad_k, 1 not")
    big = {40: 6.40e9, 100: 3.03e26, 200: 4.92e54}
    h = H2(C3)
    for k, v in big.items():
        bc = ballot_count(k)
        assert abs(bc / v - 1) < 0.01, (k, bc, v)
    for k in range(2, 201):
        bc = ballot_count(k)
        assert math.log2(bc) <= h * k + 1e-9, k
    check(True, "ballot DP matches |Bad_k| at k=40,100,200 (3 digits); |Bad_k| <= 2^(hk) for k=2..200")


def part_C():
    print("C. sigma_k class (i)")
    for k in range(2, 17):
        m = bad_mask(k)
        sig = sigma_from_mask(k, m)
        F = F_k(k)
        ok, mx = potential_ok(k, 3, sig, F)
        assert ok, k
        assert F < Fraction(C3).limit_denominator(10 ** 12) and float(F) < C3
        line = f"  k={k:2d} F_k={str(F):>5} potential certificate ok (max psi {mx})"
        if k <= 11:
            rho = karp_max_density(k, 3, sig)
            assert rho == F, (k, rho, F)
            line += f"; exact Karp rho_max = {rho}"
        # attainment: explicit periodic orbit of the first-descent word of density F
        a, d = F.numerator, F.denominator
        word = []
        prev = 0
        for j in range(1, d + 1):
            aj = a if j == d else math.floor(C3 * j) + 1
            word.append(aj - prev)
            prev = aj
        assert all(b in (0, 1) for b in word), (k, word)
        # periodic point x = c/(2^d - 3^a)
        c = 0
        for j, b in enumerate(word):
            if b:
                c = 3 * c + 2 ** j
        x = Fraction(c, 2 ** d - 3 ** a)
        M = 1 << k
        t, H = succ_table(k, 3, sig)
        orbit = []
        y = x
        for j in range(d):
            num, den = y.numerator, y.denominator
            res = (num * pow(den, -1, M)) % M
            orbit.append(res)
            assert (num % 2 == 1) == bool(word[j]), (k, j)
            y = (3 * y + 1) / 2 if num % 2 else y / 2
        assert y == x
        for j in range(d):
            s, s2 = orbit[j], orbit[(j + 1) % d]
            assert s2 in (t[s], t[s] + H), (k, j)
            if s % 2:
                assert sig[s // 2] == 1
        line += f"; attained by the period-{d} orbit of c/(2^d-3^a) = {x}"
        print(line)
    check(True, "sigma_k is class (i) with rho_max <= F_k (k=2..16), = F_k exactly (Karp, k<=11), attained (k<=16)")


def first_descent_n0(k):
    """n_0(k) = max floor(c_w/(2^d-3^a)) over first-descent words of length d <= k."""
    best = 0
    stack = [(0, 0, 0)]  # (j, a, c)
    while stack:
        j, a, c = stack.pop()
        for b in (0, 1):
            j2, a2 = j + 1, a + b
            c2 = 3 * c + 2 ** j if b else c
            if 3 ** a2 < 2 ** j2:
                best = max(best, c2 // (2 ** j2 - 3 ** a2))
            elif j2 < k:
                stack.append((j2, a2, c2))
    return best


N0_TABLE = [(2, 1), (5, 4), (8, 24)]


def part_D(N=200000):
    print(f"D. integer dynamics of T_(sigma_k), n <= {N}")
    for k in range(2, 21):
        m = bad_mask(k)
        M = 1 << k
        n0 = first_descent_n0(k)
        expect = [v for (kk, v) in N0_TABLE if k >= kk][-1]
        assert n0 == expect, (k, n0, expect)
        n = np.arange(2, N + 1, dtype=np.int64)
        x = n.copy()
        desc = np.zeros(len(n), dtype=bool)
        for j in range(1, k + 1):
            odd = x % 2 == 1
            sgn = np.where(m[x % M], -1, 1)
            x = np.where(odd, (3 * x + sgn) // 2, x // 2)
            desc |= x < n
        fails = n[~desc]
        assert len(fails) == 0 or fails.max() <= n0, (k, fails[:10])
        # every n <= n0 (hence, by induction, every n) reaches 1
        for s in range(2, max(n0, 2) + 1):
            y, steps = s, 0
            while y != 1:
                y = (3 * y + (-1 if m[y % M] else 1)) // 2 if y % 2 else y // 2
                steps += 1
                assert steps < 10 ** 5, (k, s)
        print(f"  k={k:2d} n_0={n0:3d}  non-descending n in [2,{N}]: {len(fails)} (max {fails.max() if len(fails) else '-'})")
    check(True, "n_0(k) by DFS = 1, 4, 24 at k>=2,5,8 (k<=20); every n in (n_0, N] descends within k; every n <= n_0 reaches 1")


def necklaces_by_enum(k):
    c = C3
    seen = set()
    cnt = 0
    for v in range(1 << k):
        if bin(v).count("1") <= c * k:
            continue
        rots = [((v << i) | (v >> (k - i))) & ((1 << k) - 1) for i in range(k)]
        m = min(rots)
        if m not in seen:
            seen.add(m)
            cnt += 1
    return cnt


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def necklaces_burnside(k, pred):
    tot = Fraction(0)
    for mm in range(k + 1):
        if not pred(mm):
            continue
        g = math.gcd(k, mm) if mm else k
        s = 0
        for d in range(1, g + 1):
            if g % d == 0 and k % d == 0 and mm % d == 0:
                s += phi(d) * math.comb(k // d, mm // d)
        tot += Fraction(s, k)
    assert tot.denominator == 1
    return int(tot)


NK_TABLE = {2: 1, 3: 2, 4: 2, 5: 2, 6: 5, 7: 5, 8: 6, 9: 16, 10: 19, 11: 52, 12: 70,
            13: 85, 14: 251, 16: 434, 20: 6910}


def terras_residue_of_word(k):
    """map parity word (bit i = parity of T^i) -> residue mod 2^k"""
    M = 1 << k
    table = {}
    for r in range(M):
        x, w = r, 0
        for i in range(k):
            if x % 2:
                w |= 1 << i
                x = (3 * x + 1) // 2
            else:
                x //= 2
        table[w] = r
    assert len(table) == M
    return table


def part_E():
    print("E. necklaces")
    h = H2(C3)
    for k in range(2, 21):
        ne = necklaces_by_enum(k)
        nb = necklaces_burnside(k, lambda mm, k=k: mm > C3 * k)
        assert ne == nb, (k, ne, nb)
        if k in NK_TABLE:
            assert ne == NK_TABLE[k], (k, ne, NK_TABLE[k])
    for k in range(2, 201):
        nb = necklaces_burnside(k, lambda mm, k=k: mm > C3 * k)
        assert math.log2(nb) >= h * k - math.log2(3 * k * k), k
        if k in (40, 100, 200):
            print(f"  N_{k} = {nb:.3e}")
    check(True, "N_k: enumeration = Burnside (k<=20), equals the note's table; N_k >= 2^(hk)/(3k^2) for k=2..200")
    for k in range(2, 13):
        tab = terras_residue_of_word(k)
        M = 1 << k
        sig = np.ones(M // 2, dtype=np.int64)
        t, H = succ_table(k, 3, sig)
        mask = (1 << k) - 1
        for v in range(M):
            r = tab[v]
            # rotation: drop the first letter (bit 0), append it at the end (bit k-1)
            v2 = (v >> 1) | ((v & 1) << (k - 1))
            r2 = tab[v2]
            assert r2 in (t[r], t[r] + H), (k, v)
            assert (r % 2) == (v & 1)
    check(True, "rotation of parity words is an edge of Collatz's parity graph G_0 (k<=12): necklaces are closed walks")


def class_i_batch(k, mul, signs):
    """signs: (B, 2^(k-1)) array of +-1.  Returns bool (B,): class (i)?"""
    M = 1 << k
    H = M >> 1
    B = signs.shape[0]
    wodd, wev = math.log(mul / 2), -math.log(2)
    w = np.array([wodd if s % 2 else wev for s in range(M)])
    D = np.zeros((B, M))
    rows = np.arange(B)
    changed = np.ones(B, dtype=bool)
    for it in range(M + 1):
        new = D.copy()
        for s in range(M):
            vals = D[:, s] + w[s]
            if s % 2 == 0:
                tt = np.full(B, (s // 2) % H)
            else:
                tt = ((mul * s + signs[:, s // 2]) // 2) % H
            for tgt in (tt, tt + H):
                cur = new[rows, tgt]
                new[rows, tgt] = np.maximum(cur, vals)
        changed = np.any(new > D + 1e-9, axis=1)
        D = new
        if not changed.any():
            break
    return ~changed


def part_F():
    print("F. delta_k exhaustive for k<=6")
    expect = {2: 1, 3: 2, 4: 2, 5: 4, 6: 5}
    for k in range(2, 7):
        H = 1 << (k - 1)
        found = None
        for size in range(0, H + 1):
            combos = list(itertools.combinations(range(H), size))
            for start in range(0, len(combos), 20000):
                chunk = combos[start:start + 20000]
                signs = np.ones((len(chunk), H), dtype=np.int64)
                for i, cmb in enumerate(chunk):
                    for c in cmb:
                        signs[i, c] = -1
                ok = class_i_batch(k, 3, signs)
                if ok.any():
                    found = (size, int(ok.sum()), [2 * c + 1 for c in chunk[int(np.argmax(ok))]])
                    break
            if found:
                break
        print(f"  k={k}: delta_k = {found[0]} (first optimum found: flips at residues {found[2]})")
        assert found[0] == expect[k], (k, found)
    check(True, "delta_k = 1,2,2,4,5 for k=2..6 by exhaustive search over flip sets of increasing size")


# ---------------------------------------------------------------- IHS (CP-SAT)
def find_expanding_cycles(k, mul, sig):
    """Bellman-Ford longest path with predecessors; returns list of positive cycles
    (node lists in forward order); empty list iff class (i)."""
    M = 1 << k
    t, H = succ_table(k, mul, sig)
    lw = [math.log(mul / 2) if s % 2 else -math.log(2) for s in range(M)]
    D = [0.0] * M
    pred = [-1] * M
    cycles = []
    for it in range(1, 4 * M + 1):
        changed = False
        for s in range(M):
            val = D[s] + lw[s]
            for tg in (int(t[s]), int(t[s]) + H):
                if val > D[tg] + 1e-12:
                    D[tg] = val
                    pred[tg] = s
                    changed = True
        if not changed:
            return []
        if it % 8 == 0 or it >= M:
            # look for cycles in the predecessor graph
            color = [0] * M
            for v0 in range(M):
                if color[v0]:
                    continue
                path, v = [], v0
                while v != -1 and color[v] == 0:
                    color[v] = 1
                    path.append(v)
                    v = pred[v]
                if v != -1 and color[v] == 1:
                    cyc = path[path.index(v):]
                    cyc = cyc[::-1]  # forward order: pred[x] -> x
                    a = sum(1 for x in cyc if x % 2)
                    if a * math.log(mul) > len(cyc) * math.log(2):
                        cycles.append(cyc)
                for x in path:
                    color[x] = 2
            if cycles:
                return cycles
    raise RuntimeError("no cycle extracted")


def verify_nogood(k, mul, cyc, minus):
    """cyc: node list forming a closed walk when odd nodes in `minus` take sign -1, others +1."""
    M = 1 << k
    H = M >> 1
    a = 0
    for i, s in enumerate(cyc):
        nx = cyc[(i + 1) % len(cyc)]
        if s % 2:
            a += 1
            sg = -1 if s in minus else 1
            tt = ((mul * s + sg) // 2) % H
        else:
            tt = (s // 2) % H
        if nx not in (tt, tt + H):
            return False
    return mul ** a > 2 ** len(cyc)


def ihs(k, mul, base_sign=1, objective=True, time_limit=600, seeds=(), bound=None, log_every=50):
    from ortools.sat.python import cp_model
    H = 1 << (k - 1)
    nogoods = [ng for ng in seeds]
    t0 = time.time()
    it = 0
    while True:
        it += 1
        mdl = cp_model.CpModel()
        x = [mdl.NewBoolVar(f"x{i}") for i in range(H)]  # x=1: sign differs from base
        for cyc, minus in nogoods:
            lits = []
            for s in set(cyc):
                if s % 2:
                    sg = -1 if s in minus else 1
                    lits.append(x[s // 2] if sg == base_sign else x[s // 2].Not())
            mdl.AddBoolOr(lits)
        if bound is not None:
            mdl.Add(sum(x) <= bound)
        if objective:
            mdl.Minimize(sum(x))
        sol = cp_model.CpSolver()
        sol.parameters.num_search_workers = 2
        sol.parameters.max_time_in_seconds = time_limit
        st = sol.Solve(mdl)
        if st == cp_model.INFEASIBLE:
            return {"status": "INFEASIBLE", "iterations": it, "nogoods": len(nogoods), "seconds": time.time() - t0}
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return {"status": sol.StatusName(st), "iterations": it, "nogoods": len(nogoods)}
        if objective:
            assert st == cp_model.OPTIMAL
        xs = [sol.Value(v) for v in x]
        sig = np.array([(-base_sign if xs[i] else base_sign) for i in range(H)], dtype=np.int64)
        cyc = find_expanding_cycles(k, mul, sig)
        if not cyc:
            return {"status": "CLASS_I", "iterations": it, "nogoods": len(nogoods), "size": sum(xs),
                    "flips": [2 * i + 1 for i in range(H) if xs[i]], "sig": sig, "seconds": time.time() - t0}
        for c in cyc:
            minus = {s for s in c if s % 2 and sig[s // 2] == -1}
            assert verify_nogood(k, mul, c, minus)
            nogoods.append((c, minus))
        if it % log_every == 0:
            print(f"    ... it {it}: {len(nogoods)} no-goods, LB {sum(xs) if objective else '-'}, {time.time()-t0:.0f}s", flush=True)


def simple_expanding_cycles_G0(k, mul, P):
    """all simple cycles of length <= P of the base graph (all + signs) that are expanding."""
    M = 1 << k
    sig = np.ones(M // 2, dtype=np.int64)
    t, H = succ_table(k, mul, sig)
    succ = [(int(t[s]), int(t[s]) + H) for s in range(M)]
    out = []
    for s0 in range(M):
        stack = [(s0, [s0])]
        while stack:
            v, path = stack.pop()
            for u in succ[v]:
                if u == s0:
                    a = sum(1 for z in path if z % 2)
                    if mul ** a > 2 ** len(path):
                        out.append((list(path), set()))
                elif u > s0 and u not in path and len(path) < P:
                    stack.append((u, path + [u]))
    return out


def part_G():
    print("G. delta_7, delta_8 (3n+1) by IHS with CP-SAT (independent solver)")
    for k, expect in ((7, 9), (8, 14)):
        seeds = simple_expanding_cycles_G0(k, 3, k + 3)
        r = ihs(k, 3, seeds=seeds, time_limit=600)
        print(f"  k={k}: {r['status']} size {r.get('size')} after {r['iterations']} iterations, {r['nogoods']} no-goods ({len(seeds)} seeds), {r.get('seconds',0):.0f}s")
        assert r["status"] == "CLASS_I" and r["size"] == expect, (k, r)
        rho = karp_max_density(k, 3, r["sig"])
        print(f"        optimum flips {r['flips']}, exact rho_max = {rho}")
        assert rho < Fraction(C3).limit_denominator(10 ** 9)
    check(True, "delta_7 = 9, delta_8 = 14 re-derived with a third solver (CP-SAT IHS)")


def part_H():
    print("H. 5n+-1 cube")
    c5 = math.log(2) / math.log(5)
    for k in range(2, 6):
        Hh = 1 << (k - 1)
        allsig = np.array(list(itertools.product((1, -1), repeat=Hh)), dtype=np.int64)
        ok = np.zeros(len(allsig), dtype=bool)
        for st in range(0, len(allsig), 16384):
            ok[st:st + 16384] = class_i_batch(k, 5, allsig[st:st + 16384])
        print(f"  k={k}: {len(allsig)} strategies, class (i): {int(ok.sum())}")
        assert ok.sum() == 0
    r = ihs(6, 5, objective=False, time_limit=600,
            seeds=simple_expanding_cycles_G0(6, 5, 9))
    print(f"  k=6: IHS feasibility -> {r['status']} after {r['iterations']} iterations, {r['nogoods']} no-goods")
    assert r["status"] == "INFEASIBLE"
    r = ihs(7, 5, objective=False, time_limit=600, seeds=simple_expanding_cycles_G0(7, 5, 10))
    assert r["status"] == "CLASS_I"
    rho = karp_max_density(7, 5, r["sig"])
    print(f"  k=7: class-(i) strategy found ({r['size']} flips from 5n+1), exact rho_max = {rho} < log_5 2 = {c5:.6f}")
    assert float(rho) < c5
    check(True, "5n+-1: provable class empty at k=2..6 (exhaustive k<=5, IHS UNSAT k=6), nonempty at k=7")
    t0 = time.time()
    r = ihs(7, 5, objective=True, time_limit=1200, seeds=simple_expanding_cycles_G0(7, 5, 10), log_every=100)
    print(f"  k=7: 5n+1 exact distance = {r.get('size')} of 64 ({r['status']}, {r['iterations']} it, {time.time()-t0:.0f}s)")
    assert r["status"] == "CLASS_I" and r["size"] == 29
    check(True, "5n+1's exact distance to class (i) at k=7 is 29/64")


def part_K(time_limit=1500):
    print("K. k=10 certificate")
    path = "05-knowledge/results/procgen_cubedist_20260925_k10_certificate.json.gz"
    d = json.load(gzip.open(path, "rt"))
    k = d["k"]
    assert k == 10 and d["mul"] == 3 and d["LB"] == 40
    ngs = []
    for cyc, minus in d["nogoods"]:
        assert verify_nogood(k, 3, cyc, set(minus)), cyc
        ngs.append((cyc, set(minus)))
    check(True, f"all {len(ngs)} stored no-goods re-derived as expanding closed walks")
    best = d["best"]
    H = 1 << (k - 1)
    sig = np.ones(H, dtype=np.int64)
    for r in best:
        sig[r // 2] = -1
    ok, mx = potential_ok(k, 3, sig, Fraction(5, 8))
    rho = karp_max_density(k, 3, sig)
    print(f"  upper-bound set: {len(best)} flips, potential at 5/8 ok={ok}, exact rho_max={rho}")
    assert ok and len(best) == 44 and rho == Fraction(5, 8)
    seeds = simple_expanding_cycles_G0(k, 3, k + 3)
    print(f"  own seeds: {len(seeds)} simple expanding cycles of G_0 of length <= 13")
    # The hitting-set bound is solved with HiGHS (MIP) on this script's own model:
    # the stored no-goods (each re-derived above) plus this script's own seeds.
    # (A CP-SAT attempt with 2 workers returned UNKNOWN after 3000 s.)
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 2)
    h.setOptionValue("time_limit", float(time_limit))
    inf = highspy.kHighsInf
    for i in range(H):
        h.addVar(0, 1)
        h.changeColIntegrality(i, highspy.HighsVarType.kInteger)
        h.changeColCost(i, 1.0)
    for cyc, minus in ngs + seeds:
        pos = [s // 2 for s in set(cyc) if s % 2 and s not in minus]
        neg = [s // 2 for s in set(cyc) if s % 2 and s in minus]
        idx = pos + neg
        val = [1.0] * len(pos) + [-1.0] * len(neg)
        h.addRow(1.0 - len(neg), inf, len(idx), np.array(idx, dtype=np.int32), np.array(val))
    t0 = time.time()
    h.run()
    info = h.getInfo()
    status = h.modelStatusToString(h.getModelStatus())
    print(f"  HiGHS on {len(ngs)+len(seeds)} no-goods: {status}, optimum {info.objective_function_value:.6f}, "
          f"dual bound {info.mip_dual_bound:.6f} ({time.time()-t0:.0f}s)")
    check(status == "Optimal" and info.mip_dual_bound > 39.5,
          "min hitting set of the re-derived no-goods is 40: delta_10 >= 40 (independent model and seeds)")


if __name__ == "__main__":
    parts = sys.argv[1:] or ["A", "B", "C", "D", "E", "F"]
    for p in parts:
        t0 = time.time()
        globals()["part_" + p]()
        print(f"  [part {p}: {time.time()-t0:.1f}s]", flush=True)
