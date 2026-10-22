# Audit2 (C): the two-loop kernel table of THM-4757 (4), extracted purely combinatorially from the brute-force
# enumeration of excess-1 connected (run tuple, partition) pairs (b_partitions.partitions_pruned), with no use of the
# session's kernel code.  Each such pair has a quotient multigraph that is a subdivided figure-eight, theta or dumbbell;
# suppressing degree-2 blocks, Z_pi = sign * Z_core(lengths), sign = (-1)^(#arcs traversed against the core path's
# reference direction).  We accumulate C(class) = sum mu(pi) * sign / prod(mult!) over the canonical class
# (core type, sorted lengths) and compare with the table's prediction
#   figure-eight (a,b): K8 * #orderings / 8,   theta (a,b,c): Kth * #orderings / 12,   dumbbell (a,b;c): Kdb * #orderings(a,b) / 8,
#   K8 = 8(e_a + e_b) - 12 e_a e_b (a, b even), Kth = 4(e_a+e_b+e_c) - 8 e_a e_b e_c (all even), -4 e_c (a, b odd, c even),
#   Kdb = -32 e_c (a = b = 2 mod 4, c even), else 0;  e_a = (-1)^(a/2).
# Classes whose Z_core vanishes identically for skew matrices (a loop of odd length) are reported separately.
import sys, itertools, time
from fractions import Fraction
from math import factorial
from collections import Counter, defaultdict
sys.path.insert(0, '.')
from b_partitions import multisets, partitions_pruned

def eps(a): return (-1) ** (a // 2)

def K8(a, b): return 8 * (eps(a) + eps(b)) - 12 * eps(a) * eps(b) if a % 2 == 0 and b % 2 == 0 else 0
def Kth(a, b, c):
    if a % 2 == 0 and b % 2 == 0 and c % 2 == 0: return 4 * (eps(a) + eps(b) + eps(c)) - 8 * eps(a) * eps(b) * eps(c)
    odd = [x for x in (a, b, c) if x % 2]
    if len(odd) == 2: return -4 * eps([x for x in (a, b, c) if x % 2 == 0][0])
    return 0
def Kdb(a, b, c): return -32 * eps(c) if a % 4 == 2 and b % 4 == 2 and c % 2 == 0 else 0

def n_orderings(t): return factorial(len(t)) // __import__('math').prod(factorial(v) for v in Counter(t).values())

def classify(V, arcs):
    """arcs: list of (tail, head) block ids.  Returns (type, lengths, sign) or None."""
    deg = [0] * V; inc = defaultdict(list)
    for k, (a, b) in enumerate(arcs):
        deg[a] += 1; deg[b] += 1
        inc[a].append((k, 0)); inc[b].append((k, 1))      # (arc id, which end: 0 = tail, 1 = head)
    core = [v for v in range(V) if deg[v] >= 3]
    used = set(); paths = []
    for c0 in core:
        for (k, end) in inc[c0]:
            if k in used: continue
            # traverse from c0 along arc k
            length = 0; back = 0; v = c0; arc, e = k, end
            while True:
                used.add(arc); length += 1
                tail, head = arcs[arc]
                if e == 0: nxt = head            # leaving v through the tail end: forward traversal
                else: nxt = tail; back += 1      # leaving through the head end: against the arc
                v = nxt
                if deg[v] >= 3: break
                # degree-2: continue through the other incident arc-end
                others = [(kk, ee) for (kk, ee) in inc[v] if not (kk == arc and ee == (1 - e))]
                assert len(others) == 1, (others, arc, e)
                arc, e = others[0]
            paths.append((c0, v, length, back))
    if len(core) == 1:
        assert len(paths) == 2 and deg[core[0]] == 4
        ls = sorted(p[2] for p in paths)
        sign = (-1) ** sum(p[3] for p in paths)
        return ('figure-eight', tuple(ls), sign, any(l % 2 for l in ls))
    assert len(core) == 2 and len(paths) == 3
    u, w = core
    loops = [p for p in paths if p[0] == p[1]]
    if not loops:
        # theta: orient every path u -> w
        sign = 1; ls = []
        for (s, t, l, b) in paths:
            if s == u: sign *= (-1) ** b
            else: sign *= (-1) ** (l - b)
            ls.append(l)
        return ('theta', tuple(sorted(ls)), sign, False)
    assert len(loops) == 2
    bridge = [p for p in paths if p[0] != p[1]][0]
    lu = [p for p in loops if p[0] == u][0]; lw = [p for p in loops if p[0] == w][0]
    s, t, c, b = bridge
    sign = (-1) ** lu[3] * (-1) ** lw[3] * ((-1) ** b if s == u else (-1) ** (c - b))
    a1, b1 = lu[2], lw[2]
    return ('dumbbell', (min(a1, b1), max(a1, b1), c), sign, (a1 % 2 == 1) or (b1 % 2 == 1))

def run(Lmax, Lmin=4):
    C = defaultdict(Fraction); zero_classes = set(); counts = Counter()
    for L in range(Lmin, Lmax + 1, 2):
        t0 = time.time(); np_ = 0
        for tup in multisets(L, odd=False):
            wden = 1
            for c in Counter(tup).values(): wden *= factorial(c)
            pos = []; prev = []; deg = []
            for j, l in enumerate(tup):
                for i in range(l + 1):
                    prev.append(len(pos) - 1 if i > 0 else -1)
                    deg.append(1 if i in (0, l) else 2)
                    pos.append((j, i))
            arcs_pos = [(k - 1, k) for k in range(len(pos)) if prev[k] >= 0]
            for blocks in partitions_pruned(deg, prev, True, 1, L):
                V = len(blocks)
                if L - V != 1: continue
                bid = [0] * len(pos)
                for t, b in enumerate(blocks):
                    for q in b: bid[q] = t
                arcs = [(bid[a], bid[b]) for a, b in arcs_pos]
                par = list(range(V))
                def f(a):
                    while par[a] != a: par[a] = par[par[a]]; a = par[a]
                    return a
                for a, b in arcs: par[f(a)] = f(b)
                if len({f(a) for a in range(V)}) != 1: continue
                mu = 1
                for b in blocks: mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
                typ, ls, sign, zero = classify(V, arcs)
                C[(typ, ls)] += Fraction(mu * sign, wden)
                if zero: zero_classes.add((typ, ls))
                np_ += 1
        counts[L] = np_
        print(f"L={L}: {np_} excess-1 connected partitions ({time.time() - t0:.1f}s)", flush=True)
    return C, zero_classes, counts

if __name__ == '__main__':
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    C, zc, counts = run(Lmax)
    bad = 0
    print('\nclass: brute-force coefficient C vs table prediction')
    for (typ, ls) in sorted(C, key=lambda k: (sum(k[1]), k[0], k[1])):
        v = C[(typ, ls)]
        if (typ, ls) in zc:
            print(f"  {typ:12s} {str(ls):12s} C = {str(v):>8s}   (Z_core = 0 identically: odd loop)")
            continue
        if typ == 'figure-eight': pred = Fraction(K8(*ls) * n_orderings(ls), 8)
        elif typ == 'theta':
            # Kth is symmetric in the three lengths for the stated table
            pred = Fraction(Kth(*ls) * n_orderings(ls), 12)
        else:
            a, b, c = ls
            pred = Fraction(Kdb(a, b, c) * n_orderings((a, b)), 8)
        ok = (v == pred); bad += (not ok)
        print(f"  {typ:12s} {str(ls):12s} C = {str(v):>8s}   table: {str(pred):>8s}   {'ok' if ok else 'MISMATCH'}")
    # also report table classes that never occurred (should be classes with zero prediction)
    print('\nall classes with nonzero table prediction and total length <= Lmax occur:', end=' ')
    missing = []
    for L in range(4, Lmax + 1, 2):
        for a in range(2, L - 1, 2):
            b = L - a
            if b >= a and b % 2 == 0 and K8(a, b) != 0 and ('figure-eight', (a, b)) not in C: missing.append(('figure-eight', (a, b)))
        for a in range(1, L):
            for b in range(a, L):
                c = L - a - b
                if c < b: continue
                if Kth(a, b, c) != 0 and ('theta', (a, b, c)) not in C: missing.append(('theta', (a, b, c)))
        for a in range(2, L, 2):
            for b in range(a, L, 2):
                c = L - a - b
                if c >= 1 and Kdb(a, b, c) != 0 and ('dumbbell', (a, b, c)) not in C: missing.append(('dumbbell', (a, b, c)))
    print('yes' if not missing else f'NO: {missing}')
    print('KERNEL TABLE: ' + ('ALL MATCH' if bad == 0 and not missing else f'{bad} MISMATCHES'))
