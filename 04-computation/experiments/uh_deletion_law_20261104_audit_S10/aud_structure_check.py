"""aud_structure_check.py -- independent audit code (uh_deletion_law_20261104_audit_S10). Shares no code with the session.

Reads the DUMP lines of aud_uhdel_census.c (every sole-exit-free class with hc = k in {3,4,5}, orders 6..11, with the
actual hc(T - v) for every v) and, for every class, independently:
  1. finds all one-arc-leaving sets and the minimal ones (checks uniqueness), the arc u -> l, R = V - D;
  2. checks every structural claim of THM-4874's "Setting": e_u(D) = k; s_l(R) = 1; if |R| >= 2 the contraction T/D is
     uniquely Hamiltonian with z a sole exit; l beats D - u; every other prefix vertex beats all of D; D - u has
     condensation (trivial chain) => W with W strong, |W| >= 3; A = N-(u) & W, u beats W - A; u -> w_1 if p >= 1;
     p = 0 only if A != W; (W, A) satisfies the block conditions and has type k;
  3. predicts hc(T - v) for every v from the table (1) and compares with the census value;
  4. identifies the block up to isomorphism among the 8 blocks of aud_blocks.out;
  5. tallies the joint (k1, k2, |Aut|) per order and per block and compares with HYP-9467's formulas;
  6. checks (4): for hc = 3, |Aut| > 1 iff k2 = 0;
  7. (optionally, n <= NRECHECK) recomputes every hc(T - v) with this file's own DP.
Usage: python aud_structure_check.py dumpfile [NRECHECK]"""
import sys
from itertools import permutations
from collections import defaultdict, Counter


def fib(m):
    a, b = 0, 1
    if m < 0:
        return (-1) ** (m + 1) * fib(-m)
    for _ in range(m):
        a, b = b, a + b
    return a


def sub_out(rows, verts, converse=False):
    """out-neighbour bitmasks of T[verts] (or of its converse), in the order of verts."""
    k = len(verts)
    out = [0] * k
    for i, v in enumerate(verts):
        m = 0
        for j, w in enumerate(verts):
            if i != j:
                arc = ((rows[v] >> w) & 1) if not converse else ((rows[w] >> v) & 1)
                if arc:
                    m |= 1 << j
        out[i] = m
    return out


def paths_from(out, s):
    """list c[t] = #Hamiltonian paths of the tournament 'out' (bitmask rows) starting at s and ending at t."""
    k = len(out)
    full = (1 << k) - 1
    if k == 1:
        return [1]
    dp = [[0] * k for _ in range(1 << k)]
    dp[1 << s][s] = 1
    for S in range(1, full + 1):
        if not (S >> s) & 1:
            continue
        row = dp[S]
        for v in range(k):
            c = row[v]
            if not c:
                continue
            nb = out[v] & ~S & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                dp[S | (1 << w)][w] += c
    return dp[full]


def ham_paths_matrix(rows, verts):
    """P[(s,t)] = #Hamiltonian paths of T[verts] from s to t."""
    verts = list(verts)
    out = sub_out(rows, verts)
    P = {}
    for i, s in enumerate(verts):
        c = paths_from(out, i)
        for j, t in enumerate(verts):
            P[(s, t)] = c[j]
    return P


def end_count(rows, verts, x):
    verts = list(verts)
    out = sub_out(rows, verts, converse=True)
    return sum(paths_from(out, verts.index(x)))


def start_count(rows, verts, x):
    verts = list(verts)
    out = sub_out(rows, verts)
    return sum(paths_from(out, verts.index(x)))


def hcycles(rows, verts):
    verts = list(verts)
    if len(verts) < 3:
        return 0
    out = sub_out(rows, verts)
    c = paths_from(out, 0)
    return sum(c[j] for j in range(1, len(verts)) if (out[j] >> 0) & 1)


def leaving(rows, n, Dmask):
    return sum(bin(rows[v] & ~Dmask & ((1 << n) - 1)).count("1") for v in range(n) if (Dmask >> v) & 1)


def strong_components_in_order(rows, verts):
    """strong components of T[verts] in the order C_1 => C_2 => ... (tournament: components are linearly ordered)."""
    verts = list(verts)
    vs = set(verts)
    reach = {}
    for s in verts:
        seen = {s}
        st = [s]
        while st:
            v = st.pop()
            for w in vs:
                if w not in seen and (rows[v] >> w) & 1:
                    seen.add(w)
                    st.append(w)
        reach[s] = seen
    comps = []
    done = set()
    for s in verts:
        if s in done:
            continue
        comp = {t for t in reach[s] if s in reach[t]}
        comps.append(comp)
        done |= comp
    # order: C before C' iff C reaches C'; sort by size of reach set (descending)
    comps.sort(key=lambda C: -len(reach[next(iter(C))]))
    return comps


def canon_pair(rows, Wl, A):
    """canonical key of the pair (T[W], A) by brute force over permutations."""
    k = len(Wl)
    best = None
    for p in permutations(range(k)):
        code = []
        for i in range(k):
            # vertex at new position i is Wl[p[i]]
            v = Wl[p[i]]
            code.append(tuple(1 if (rows[v] >> Wl[p[j]]) & 1 else 0 for j in range(k)))
        a = tuple(1 if Wl[p[i]] in A else 0 for i in range(k))
        key = (k, tuple(code), a)
        if best is None or key < best:
            best = key
    return best


def block_conditions(rows, Wl, A):
    Wset = set(Wl)
    k = len(Wl)
    if len(A) < 2:
        return False
    for v in Wl:
        if v not in A and sum(1 for w in Wl if (rows[v] >> w) & 1) < 2:
            return False
    for zm in range(1, 1 << k):
        Z = {Wl[i] for i in range(k) if (zm >> i) & 1}
        lv = sum(1 for v in Z for w in Wset - Z if (rows[v] >> w) & 1)
        if lv + len(Z & A) == 1:
            return False
    # strong
    comps = strong_components_in_order(rows, Wl)
    return len(comps) == 1


# the 8 blocks of aud_blocks.out (out-neighbour lists, A), named as in THM-4874
def rows_from_lists(lists):
    return [sum(1 << w for w in L) for L in lists]


W4L = [[3], [0], [0, 1], [1, 2]]   # aud_blocks labels: 0 = d, 1 = c, 2 = b, 3 = a  (W4: a->b->c->d->a, a->c, b->d)
BLOCKS = [
    ("k3 (C3, all)", rows_from_lists([[2], [0], [1]]), {0, 1, 2}, 3),
    ("k3 (W4, {c,d})", rows_from_lists(W4L), {0, 1}, 3),
    ("k4 (W4, {b,c,d})", rows_from_lists(W4L), {0, 1, 2}, 4),
    ("k4 (W4, {a,c,d})", rows_from_lists(W4L), {0, 1, 3}, 4),
    ("k5 (W4, all)", rows_from_lists(W4L), {0, 1, 2, 3}, 5),
    ("k5 (AT5, {4,5})", rows_from_lists([[4], [0], [0, 1], [0, 1, 2], [1, 2, 3]]), {0, 1}, 5),
    ("k5 (W5, A1={0,1})", rows_from_lists([[3, 4], [0], [0, 1], [1, 2], [1, 2, 3]]), {0, 1}, 5),
    ("k5 (W5, A2={1,4})", rows_from_lists([[3, 4], [0], [0, 1], [1, 2], [1, 2, 3]]), {1, 4}, 5),
]
BLOCK_KEYS = {}
BLOCK_DATA = {}
for name, r, A, k in BLOCKS:
    Wl = list(range(len(r)))
    BLOCK_KEYS[canon_pair(r, Wl, A)] = name
    # delta and h0 from this file's own DP
    deltas = []
    for v in Wl:
        rest = [x for x in Wl if x != v]
        deltas.append(sum(end_count(r, rest, x) for x in A if x != v))
    P = ham_paths_matrix(r, Wl)
    h0 = sum(P[(s, t)] for s in Wl for t in Wl if s != t and s not in A and t in A)
    BLOCK_DATA[name] = (deltas.count(1), deltas.count(2), h0, len(Wl), len(A) == len(Wl), k)

# HYP-9467 predicted joint table (n >= 7): block -> (k1, k2, aut, F-index offset as function of n)
PRED = {
    "k3 (C3, all)": ((3, 0, 3), lambda n: fib(2 * n - 10)),
    "k3 (W4, {c,d})": ((3, 1, 1), lambda n: fib(2 * n - 11)),
    "k4 (W4, {a,c,d})": ((2, 1, 1), lambda n: fib(2 * n - 11)),
    "k4 (W4, {b,c,d})": ((2, 2, 1), lambda n: fib(2 * n - 11)),
    "k5 (W4, all)": ((2, 0, 1), lambda n: fib(2 * n - 12)),
    "k5 (AT5, {4,5})": ((3, 0, 1), lambda n: fib(2 * n - 13)),
    "k5 (W5, A1={0,1})": ((2, 1, 1), lambda n: fib(2 * n - 13)),
    "k5 (W5, A2={1,4})": ((2, 2, 1), lambda n: fib(2 * n - 13)),
}


def main(path, nrecheck=0):
    classes = []
    with open(path) as f:
        for line in f:
            if not line.startswith("DUMP"):
                continue
            left, right = line.split("|")
            parts = left.split()
            mode, n, h, aut = parts[1], int(parts[2]), int(parts[3]), int(parts[4])
            rows = [int(x) for x in parts[5:5 + n]]
            hv = [int(x) for x in right.split()]
            classes.append((mode, n, h, aut, rows, hv))
    print(f"read {len(classes)} sole-exit-free classes with hc in {{3,4,5}} from {path}")
    print("block data recomputed here (#delta=1, #delta=2, h0, |W|, A=W, k):")
    for name in BLOCK_DATA:
        print(f"  {name}: {BLOCK_DATA[name]}")
    nfail = 0
    nvert = 0
    per_n_block = defaultdict(Counter)
    joint = defaultdict(Counter)
    aut_k2_viol = 0
    setting_viol = Counter()
    coverage = Counter()
    nrechecked = [0]
    for mode, n, k, aut, rows, hv in classes:
        V = list(range(n))
        full = (1 << n) - 1
        if nrecheck and n <= nrecheck:
            nrechecked[0] += n + 1
            for v in V:
                hv2 = hcycles(rows, [x for x in V if x != v])
                if hv2 != hv[v]:
                    print("RECHECK MISMATCH", n, rows, v, hv[v], hv2)
                    nfail += 1
            if hcycles(rows, V) != k:
                print("RECHECK hc MISMATCH", n, rows)
                nfail += 1
        one_arc = [D for D in range(1, full) if leaving(rows, n, D) == 1]
        minimal = [D for D in one_arc if not any((E & D) == E and E != D for E in one_arc)]
        if len(minimal) != 1:
            in_range = (k == 3 and n >= 6) or (k in (4, 5) and n >= 7)
            if in_range:
                setting_viol["unique minimal one-arc-leaving set"] += 1
            print(f"no unique minimal one-arc-leaving set ({len(minimal)} minimal; arc-strong connectivity >= 2 if 0):"
                  f" n={n} hc={k} |Aut|={aut} k1={hv.count(1)} k2={hv.count(2)} rows={rows}"
                  f" {'IN THE THEOREM RANGE -- VIOLATION' if in_range else '(outside the theorem range)'}")
            continue
        Dm = minimal[0]
        D = [v for v in V if (Dm >> v) & 1]
        R = [v for v in V if not (Dm >> v) & 1]
        (u, l), = [(v, w) for v in D for w in R if (rows[v] >> w) & 1]
        if len(D) < 3:
            setting_viol["|D| >= 3"] += 1
        # every arc between R and D points into D except u -> l
        for x in R:
            for y in D:
                if (x, y) != (l, u) and not (rows[x] >> y) & 1:
                    setting_viol["R -> D arcs"] += 1
        # e_u(D) = k
        eu = end_count(rows, D, u)
        if eu != k:
            setting_viol["e_u(D) = k"] += 1
        sl = 1 if len(R) == 1 else start_count(rows, R, l)
        if sl != 1:
            setting_viol["s_l(R) = 1"] += 1
        if len(R) >= 2:
            # contraction Q = T/D: vertices R + z; z -> l only, R - l -> z
            z = n  # new label
            qrows = {x: rows[x] & ~Dm for x in R}
            for x in R:
                if x != l:
                    qrows[x] |= 1 << z
            qrows[z] = 1 << l
            Qrows = [0] * (n + 1)
            for x, m in qrows.items():
                Qrows[x] = m
            Qv = R + [z]
            if hcycles(Qrows, Qv) != 1:
                setting_viol["Q uniquely Hamiltonian"] += 1
        # condensation of D - u
        Du = [v for v in D if v != u]
        comps = strong_components_in_order(rows, Du)
        Wset = comps[-1]
        chain = [next(iter(C)) for C in comps[:-1]]
        if any(len(C) != 1 for C in comps[:-1]) or len(Wset) < 3:
            setting_viol["chain + block form"] += 1
            print("NOT chain+block", n, rows, [sorted(C) for C in comps])
            continue
        p = len(chain)
        Wl = sorted(Wset)
        A = {x for x in Wl if (rows[x] >> u) & 1}
        if any(not (rows[u] >> x) & 1 for x in Wl if x not in A):
            setting_viol["u beats W - A"] += 1
        if p >= 1 and not (rows[u] >> chain[0]) & 1:
            setting_viol["u -> w_1"] += 1
        if p == 0 and len(A) == len(Wl):
            setting_viol["p = 0 only if A != W"] += 1
        if not block_conditions(rows, Wl, A):
            setting_viol["block conditions"] += 1
        typ = sum(end_count(rows, Wl, x) for x in A)
        if typ != k:
            setting_viol["type = k"] += 1
        bkey = canon_pair(rows, Wl, A)
        bname = BLOCK_KEYS.get(bkey, "UNKNOWN BLOCK")
        if bname == "UNKNOWN BLOCK":
            setting_viol["block in the list of 8"] += 1
        # h0
        PW = ham_paths_matrix(rows, Wl)
        h0 = sum(PW[(s, t)] for s in Wl for t in Wl if s != t and s not in A and t in A)
        # predictions
        for v in V:
            if v in R and v != l:
                Rv = [x for x in R if x != v]
                pred = k * (1 if len(Rv) == 1 else start_count(rows, Rv, l))
                pos = "prefix"
            elif v == u:
                pred, pos = 0, "u"
            elif v == l:
                pred = 0 if len(R) >= 2 else (k if p >= 1 else h0)
                pos = "l"
            elif v in chain:
                pred, pos = k, "chain"
            else:
                rest = [x for x in Wl if x != v]
                pred = sum(end_count(rows, rest, x) for x in A if x != v)
                pos = "W"
            if pos == "l":
                pos = "l, |R|>=2" if len(R) >= 2 else ("l, R={l}, p>=1" if p >= 1 else "l, R={l}, p=0")
            if pos == "prefix":
                pos = "prefix v != l"
            coverage[pos] += 1
            nvert += 1
            if pred != hv[v]:
                nfail += 1
                print("MISMATCH", n, k, rows, "v", v, pos, "pred", pred, "actual", hv[v])
        k1 = hv.count(1)
        k2 = hv.count(2)
        per_n_block[(n, k)][(bname, k1, k2, aut)] += 1
        joint[(n, k)][(k1, k2, aut)] += 1
        if k == 3 and ((aut > 1) != (k2 == 0)):
            aut_k2_viol += 1
            print("(4) VIOLATION", n, rows, aut, k2)
    print(f"\nvertex deletions compared with table (1): {nvert}; mismatches: {nfail}")
    print(f"coverage of the cases of (1): {dict(coverage)}")
    print(f"hc(T) and hc(T - v) values recomputed in Python (orders <= {nrecheck}): {nrechecked[0]} (mismatches are counted above)")
    print(f"Setting-claim violations: {dict(setting_viol) if setting_viol else 'none'}")
    print(f"(4) hc = 3: |Aut| > 1 iff k2 = 0 -- violations: {aut_k2_viol}")
    print("\nper order and hc: (block, k1, k2, |Aut|): count   [HYP-9467 prediction for n >= 7]")
    allok = True
    for (n, k) in sorted(per_n_block):
        print(f" n={n} hc={k}:")
        for key, c in sorted(per_n_block[(n, k)].items()):
            bname, k1, k2, aut = key
            note = ""
            if n >= 7 and bname in PRED:
                (pk1, pk2, paut), f = PRED[bname]
                ok = (pk1, pk2, paut) == (k1, k2, aut) and f(n) == c
                note = f"  predicted ({pk1},{pk2},{paut}) x {f(n)}  {'OK' if ok else 'FAIL'}"
                allok &= ok
            print(f"   {key}: {c}{note}")
        # every predicted block present?
        if n >= 7:
            for bname, ((pk1, pk2, paut), f) in PRED.items():
                if int(bname[1]) == k and f(n) > 0 and not any(key[0] == bname for key in per_n_block[(n, k)]):
                    print(f"   MISSING block {bname}: predicted {f(n)}")
                    allok = False
    print(f"\nall HYP-9467 joint predictions at the orders >= 7 present here hold: {allok}")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 0)
