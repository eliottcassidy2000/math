"""aud_followup.py - follow-up audit of the post-checkpoint commits of opus-2026-11-08-S86
(ebe3096762 PREREGISTRATION-2, 1a1606cc66 coverage_products, 1d57e6b23c ternary_attainment).

Independent code (does not import or copy coverage_products.py / ternary_attainment.py); reads the session's .out
files only as data.

Construction sets (every piece is an explicit strong tournament generated and checked by the audit's own C code):
  CEN : every strong class of order <= 9 (aud_fewcycle census mode on kind-pasteur's classgen lists), full alpha vectors
  FC  : every strong class with alpha_1 <= 80 (aud_fewcycle 80 0 ...), full alpha vectors (alpha_1..alpha_7)
  IC  : every strong intercyclic class with alpha_1 <= 8191 (aud_fewcycle 8191 1 ..., complete: order 17 EMPTY)
  APEX: one-apex values (reduced words of length <= L) by an incremental DFS, plus Session 16's two rays, used only
        to mirror the session's universe for counts beyond alpha_1 = 8191
Closure: a value x = 1 (mod b) is attained iff x is a piece value or x = s * l with s, l > 1, one of them a piece
value and the other attained (divisor-pair DP over increasing x; a different algorithm from the session's
forward-multiplication closure).
usage: python aud_followup.py WORKDIR SESSIONDIR > aud_followup.out
"""
import os
import re
import sys
from math import isqrt


def load_census(work):
    """all classes of orders 3..9: (n, strong, alpha tuple a1..a7, H)"""
    rows = []
    for n in range(1, 10):
        for line in open(os.path.join(work, f"census_{n}.txt")):
            t = line.split()
            rows.append((int(t[0]), int(t[1]), tuple(int(x) for x in t[2:9]), int(t[9])))
    return rows


def load_inv(path):
    out = []
    for line in open(path):
        t = line.split()
        out.append((int(t[0]), tuple(int(x) for x in t[1:8]), int(t[8])))
    return out


def hb(a, b):
    return 1 + sum(c * b ** (k + 1) for k, c in enumerate(a))


def one_apex_values(L, cap):
    """values sum_{i<j, w_i=1, w_j=0} E[j-i] over words starting with 1 (E = 1, 1, 2, 4, ...), recorded at every
    terminal 0; incremental state: alpha, s = amount added by a 0 at the next position, last letter"""
    vals = set()
    stack = [(1, 0, 1, 1)]           # (length, alpha, s, last) after the first letter '1': s = E[1] = 1
    while stack:
        ln, alpha, s, last = stack.pop()
        if ln == L:
            continue
        # next position p = ln; contributions for a 0 placed there = s
        for x in (0, 1):
            na = alpha + (s if x == 0 else 0)
            if na > cap:
                continue
            if x == 0:
                vals.add(na)
            # s for the position after p: every earlier 1 moves one step further away: 2s - [last letter = 1];
            # plus E[1] = 1 if the new letter is 1
            ns = 2 * s - (1 if last == 1 else 0) + (1 if x == 1 else 0)
            stack.append((ln + 1, na, ns, x))
    return vals


def closure_divisor(piece_values, b, M):
    """attained alpha-sums m <= M in base b, with a witness factor for each attained value"""
    V = 1 + b * M
    atoms = {v for v in piece_values if 1 < v <= V}
    reach = {1: None}
    for m in range(1, M + 1):
        x = 1 + b * m
        if x in atoms:
            reach[x] = (x, 1)
            continue
        r = isqrt(x)
        for s in range(1 + b, r + 1, b):
            if x % s:
                continue
            l = x // s
            if (s in atoms and l in reach) or (l in atoms and s in reach):
                reach[x] = (s, l) if s in atoms else (l, s)
                break
    att = {(x - 1) // b for x in reach}
    return att, reach, atoms


def factor(x, reach):
    f = []
    while x > 1:
        a, rest = reach[x]
        f.append(a)
        x = rest
    return f


def main():
    work, sdir = sys.argv[1], sys.argv[2]
    cen = load_census(work)
    fc = load_inv(os.path.join(work, "fc80.txt"))
    ic = load_inv(os.path.join(work, "ic8191.txt"))
    CEN = {(n, a) for (n, st, a, H) in cen if st and n >= 3}
    FC = {(n, a) for (n, a, H) in fc}
    IC = {(n, a) for (n, a, H) in ic}
    ncen = sum(1 for (n, st, a, H) in cen if st and n >= 3)
    print("1. Pieces (audit's own explicit strong tournaments):")
    print(f"   CEN strong classes of orders 3-9: {ncen} classes; {len(CEN)} distinct (order, packing vector) pairs, "
          f"{len({a for n, a in CEN})} distinct vectors")
    print(f"   FC  strong classes with alpha_1 <= 80: {len(fc)} classes; {len(FC)} distinct (order, vector) pairs, "
          f"{len({a for n, a in FC})} distinct vectors")
    print(f"   IC  strong intercyclic classes with alpha_1 <= 8191: {len(ic)} classes; {len(IC)} distinct (order, vector) "
          f"pairs, {len({a[0] for n, a in IC})} distinct alpha_1")

    # ---- (ii) independent attainment, bases 2 and 3
    print()
    print("2. Attainment by products of CEN u FC u IC (independent construction, divisor-pair closure):")
    pieces = {}
    for src, S in (("CEN", CEN), ("FC", FC), ("IC", IC)):
        for n, a in S:
            pieces.setdefault(a, (src, n))
    res = {}
    for b, M in ((2, 1800), (3, 2200)):
        vals = {hb(a, b) for a in pieces}
        att, reach, atoms = closure_divisor(vals, b, M)
        unc = [m for m in range(M + 1) if m not in att]
        first = next(m for m in unc if m not in (3, 10))
        res[b] = (att, reach, atoms, first)
        print(f"   base {b}: uncovered alpha-sums m <= {M}: first {unc[:12]}; first uncovered m other than 3, 10: "
              f"{first} (value {1 + b * first}); every value 1 + {b}m with m < {first}, m not 3 or 10, is attained")
    # restricted re-check: CEN only (all strong classes of order <= 9) for values up to 1000
    for b in (2, 3):
        vals = {hb(a, b) for n, a in CEN}
        M = 1000 // b + 1
        att, _, _ = closure_divisor(vals, b, M)
        unc = [m for m in range(M + 1) if m not in att]
        print(f"   check with CEN alone (strong classes of order <= 9), base {b}: uncovered alpha-sums m <= {M}: {unc}")
    # direct check: every class of order <= 9 (no products at all)
    for b in (2, 3):
        direct = {hb(a, b) for (n, st, a, H) in cen}
        gaps = (7, 21) if b == 2 else (10, 31)
        miss = [x for x in range(1, max(direct) + 1) if (x - 1) % b == 0 and x not in direct and x not in gaps]
        print(f"   direct H_{b} values of all classes of order <= 9 (no products): every value = 1 (mod {b}) up to "
              f"{miss[0] - b} except {gaps} is attained; first further misses {miss[:5]}; max {max(direct)}")

    # ---- the session's construction universe (a) + (b) + (c), with correct vectors, counts to m = 10000
    print()
    print("3. The session's universe (a) strong order <= 9, (b) alpha_1 <= 80, (c) one-apex length <= 19 + rays,")
    print("   with CORRECT full vectors, counts to m = 10000 (ternary_attainment.out reports 1183 and 1698):")
    apex19 = one_apex_values(19, 70000)
    rays = {3 * 2 ** s for s in range(1, 20)} | {9 * 2 ** s - 2 for s in range(0, 20)}
    inter = {a for a in apex19 | rays if 1 <= a <= 70000}
    icvals = {a[0] for n, a in IC}
    print(f"   intercyclic values (one-apex length <= 19 and rays, <= 70000): {len(inter)}; "
          f"agree with the exhaustive IC below 8192: {({a for a in inter if a < 8192}) == ({a for a in icvals if a < 8192})}")
    vecs_ok = {a for n, a in CEN} | {a for n, a in FC}
    vecs_trunc = {a for n, a in CEN} | {a[:4] + (0, 0, 0) for n, a in FC}       # the session's v[:4] truncation
    fake = vecs_trunc - vecs_ok
    print(f"   distinct vectors: correct {len(vecs_ok)}; with the session's v[:4] truncation {len(vecs_trunc)} "
          f"({len(fake)} truncated vectors that are not packing vectors of their classes)")
    for b in (2, 3):
        print(f"   base {b}: smallest truncated (unrealised) piece value {min(hb(a, b) for a in fake)}")
    keep = {}
    for label, vecs in (("correct", vecs_ok), ("truncated", vecs_trunc)):
        out = []
        for b in (2, 3):
            vals = {hb(a, b) for a in vecs} | {1 + b * a for a in inter}
            att, reach, atoms = closure_divisor(vals, b, 10000)
            keep[(label, b)] = (att, reach)
            unc = [m for m in range(10001) if m not in att]
            first = next(m for m in unc if m not in (3, 10))
            out.append(f"base {b}: {len(unc)} uncovered m <= 10000, first beyond 3, 10: {first}")
        print(f"   {label:9s}: " + "; ".join(out))
    for b in (2, 3):
        diff = sorted(keep[("truncated", b)][0] - keep[("correct", b)][0])
        print(f"   base {b}: alpha-sums <= 10000 'covered' only through truncated (unrealised) piece values: {diff} "
              f"{[factor(1 + b * m, keep[('truncated', b)][1]) for m in diff]}")
        lost = sorted(keep[("correct", b)][0] - keep[("truncated", b)][0])
        print(f"   base {b}: alpha-sums <= 10000 covered by the genuine alpha_5 > 0 pieces but lost under truncation: "
              f"{lost} {[factor(1 + b * m, keep[('correct', b)][1]) for m in lost]}")

    # ---- (iii) the witnesses printed in ternary_attainment.out
    print()
    print("4. Witness factorisations printed in ternary_attainment.out, each factor matched to an explicit piece:")
    txt = open(os.path.join(sdir, "ternary_attainment.out")).read()
    base = None
    for line in txt.splitlines():
        mb = re.match(r"base (\d+):", line)
        if mb:
            base = int(mb.group(1))
        mw = re.search(r"witness m = (\d+): 1 \+ (\d+)\*\d+ = (\d+) = product of piece values \[([0-9, ]+)\]", line)
        if mw:
            m, b, v = int(mw.group(1)), int(mw.group(2)), int(mw.group(3))
            fs = [int(x) for x in mw.group(4).split(",")]
            prod = 1
            for x in fs:
                prod *= x
            desc = []
            for x in fs:
                hit = None
                for (src, S) in (("CEN", CEN), ("FC", FC), ("IC", IC)):
                    cand = sorted((n, a) for n, a in S if hb(a, b) == x)
                    if cand:
                        n, a = cand[0]
                        hit = f"{x} = H_{b} of a strong {src} class of order {n}, alpha = {a[:4]}"
                        break
                desc.append(hit or f"{x}: NO REALISING PIECE FOUND")
            print(f"   base {b}, m = {m}: product {prod} {'=' if prod == v else '!='} {v}")
            for d in desc:
                print(f"      {d}")
    # pure-Python re-derivation of the two large single-piece base-3 witnesses from the class codes
    for target in (3001, 5398):
        hit = None
        for line in open(os.path.join(work, "census_9.txt")):
            t = line.split()
            if t[1] == "1" and hb(tuple(int(x) for x in t[2:9]), 3) == target:
                hit = t
                break
        hx = hit[15]
        w = [int(hx[16 * (4 - i):16 * (4 - i) + 16], 16) for i in range(5)]
        n = 9
        adj = [[False] * n for _ in range(n)]
        bit = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (w[bit >> 6] >> (bit & 63)) & 1:
                    adj[i][j] = True
                else:
                    adj[j][i] = True
                bit += 1
        cyc = []                                   # odd cycles as vertex sets, each found once from its least vertex

        def dfs(start, v, seen, length):
            for u in range(start + 1, n):
                if adj[v][u] and not (seen >> u) & 1:
                    dfs(start, u, seen | (1 << u), length + 1)
            if adj[v][start] and length >= 3 and length % 2 == 1:
                cyc.append(seen)
        for s in range(n):
            dfs(s, s, 1 << s, 1)
        a1 = len(cyc)
        a2 = sum(1 for i in range(a1) for j in range(i + 1, a1) if not cyc[i] & cyc[j])
        a3 = sum(1 for i in range(a1) for j in range(i + 1, a1) if not cyc[i] & cyc[j]
                 for k in range(j + 1, a1) if not (cyc[k] & (cyc[i] | cyc[j])))
        c_alpha = tuple(int(x) for x in hit[2:5])
        print(f"   pure-Python re-derivation of an order-9 piece with H_3 = {target}: alpha = ({a1}, {a2}, {a3}) "
              f"(C census: {c_alpha}, {'agree' if c_alpha == (a1, a2, a3) else 'DISAGREE'}), "
              f"H_3 = {1 + 3 * a1 + 9 * a2 + 27 * a3}, H = {1 + 2 * a1 + 4 * a2 + 8 * a3}, out-masks "
              f"{[sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)]}")
    claims_gap = [ln for ln in txt.splitlines() if re.search(r"\bgap|never|impossible|not attain", ln, re.I)]
    print("   lines of ternary_attainment.out that assert a gap / non-attainment:", claims_gap if claims_gap else "none")

    # ---- (iv) the numbers quoted in the PREREGISTRATION-2 outcome
    print()
    print("5. PREREGISTRATION-2 universe: (i) alpha_1 <= 80 strong pieces (full vectors), (ii) one-apex reduced length")
    print("   <= 21 (script) / <= 23 (pre-registration text), (iii) the rays; uncovered alpha-sums by base:")
    apex21 = one_apex_values(21, 100000)
    apex23 = one_apex_values(23, 100000)
    print(f"   one-apex values <= 100000 from reduced length <= 21 and <= 23 coincide: {apex21 == apex23}")
    inter2 = {a for a in apex21 | rays if 1 <= a <= 100000}
    fcvecs = {a for n, a in FC}
    plan = {2: 100000, 3: 100000, 4: 20000, 5: 20000, 6: 20000, 7: 20000, 8: 20000}
    for b, M in plan.items():
        vals = {hb(a, b) for a in fcvecs} | {1 + b * a for a in inter2}
        att, _, _ = closure_divisor(vals, b, M)
        low = [m for m in range(81) if m not in att]
        unc = [m for m in range(81, M + 1) if m not in att]
        print(f"   b = {b}: m <= 80 uncovered {low}; uncovered in (80, {M}]: {len(unc)}, first {unc[:3]}")
    print("   the unrun in-sample analogue (strong pieces with alpha_1 <= K plus the intercyclic families):")
    for K in (20, 40):
        sub = {a for a in fcvecs if a[0] <= K}
        for b in (2, 3):
            vals = {hb(a, b) for a in sub} | {1 + b * a for a in inter2}
            att, _, _ = closure_divisor(vals, b, 80)
            unc = [m for m in range(K + 1, 81) if m not in att]
            print(f"     K = {K}, b = {b}: uncovered alpha-sums in ({K}, 80]: {len(unc)} {unc}")


if __name__ == "__main__":
    main()
