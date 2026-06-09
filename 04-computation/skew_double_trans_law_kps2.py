#!/usr/bin/env python3
"""
skew_double_trans_law_kps2.py — kind-pasteur-2026-06-09-S2, Branch I part 4 (HYP-2357)

The doubling law for trans under the skew-Sylvester double D(T) (THM-447):
  D(T) dominance blocks [[M, M+I], [M-I, -M]]  (copy 1 = T, copy 2 = T^op,
  cross arcs i->j' iff (i->j in T or i=j), i'->j iff i->j in T).

Claims to test:
  L1 (lower):  trans(D(T)) >= trans(T) + 1   [sink-twin append: max transitive
       chain v1<...<vk in copy 1, then vk' is beaten by every vi (vi->vk' iff
       vi->vk or i=k), so append vk' as new sink]
  L2 (upper):  trans(D(T)) <= 2*trans(T)     [any transitive sub splits into
       S1 (copy 1, transitive in T) and S2 (copy 2, transitive in T^op), and
       trans(T^op) = trans(T)]
  L3 (the "+1 law"): is trans(D(T)) == trans(T) + 1 always? for strongly
       connected T? LITERATURE FORCES FAILURE SOMEWHERE: trans(QR7) = 3 but
       D(QR7) has 14 vertices and R(5) = 14 forces trans(D(QR7)) >= 5 = trans+2.
       Where does it first fail?

Tests:
  1. Exhaustive on iso classes n=3..5: trans(T), trans(D(T)), delta histogram,
     exact counterexamples to delta == 1, strong-connectivity refinement.
  2. Exhaustive on all 32768 labeled n=6 tournaments: delta histogram x strong
     connectivity; exact smallest counterexamples.
  3. Specific large: QR7, ST13 (circulant found in script A), Galois ST27,
     tower cores T7, T15, T31 (D on 14/30/62 vertices), T63 if budget allows.
  4. Tower-step growth (bordered doubling n -> 2n+1) is reported in script A;
     here we contrast D(T_m) (2m vertices) with T_{2m+1} (the actual tower step).

Output: 05-knowledge/results/skew_double_trans_law_kps2.out
"""
import sys, time
from collections import Counter
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import (D_skew, M_of, A_of, scores, iso_classes,
    normalize_first_row, core_tournament, is_skew_hadamard)
from erdos_moser_trans_tower_kps2 import (TransSolver, masks_of, trans_of, paley,
    galois_27, circulant, double_S, build_tower)

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\skew_double_trans_law_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def strongly_connected(A):
    n = A.shape[0]
    for M in (A, A.T):
        seen = {0}
        stack = [0]
        while stack:
            u = stack.pop()
            for v in range(n):
                if M[u, v] and v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(seen) != n:
            return False
    return True

def bits_of(A):
    n = A.shape[0]
    return ''.join(str(A[i, j]) for i in range(n) for j in range(i + 1, n))

def main():
    t_start = time.time()
    w('=== skew_double_trans_law_kps2 — trans under the skew-Sylvester double (HYP-2357) ===')
    w('')

    # ---------- PART 1: iso classes n=3..5 ----------
    w('--- PART 1: all iso classes n=3..5 ---')
    w(f'{"n":>2} {"idx":>3} {"scores":>16} {"SC?":>4} {"trans":>6} {"transD":>7} {"delta":>6}')
    viol_L1 = viol_L2 = 0
    delta_hist = Counter()
    for n in (3, 4, 5):
        for idx, A in enumerate(iso_classes(n)):
            tr, _ = trans_of(A)
            Ad, _ = D_skew(A)
            trD, _ = trans_of(Ad)
            d = trD - tr
            sc = strongly_connected(A)
            delta_hist[(n, d, sc)] += 1
            viol_L1 += (trD < tr + 1)
            viol_L2 += (trD > 2 * tr)
            w(f'{n:>2} {idx:>3} {str(scores(A)):>16} {"Y" if sc else "N":>4} '
              f'{tr:>6} {trD:>7} {d:>6}')
    w('')
    w(f'L1 (trans(D) >= trans+1) violations: {viol_L1}')
    w(f'L2 (trans(D) <= 2*trans) violations: {viol_L2}')
    w('delta histogram by (n, delta, strongly_connected): '
      f'{dict(sorted(delta_hist.items()))}')
    w('')

    # ---------- PART 2: exhaustive labeled n=6 ----------
    w('--- PART 2: exhaustive labeled n=6 (32768 tournaments) ---')
    t0 = time.time()
    pairs6 = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    hist6 = Counter()
    examples = {}
    v1 = v2 = 0
    for bits in range(1 << 15):
        A = np.zeros((6, 6), dtype=np.int64)
        for k, (i, j) in enumerate(pairs6):
            if (bits >> k) & 1:
                A[i, j] = 1
            else:
                A[j, i] = 1
        tr, _ = trans_of(A)
        Ad, _ = D_skew(A)
        trD, _ = trans_of(Ad)
        d = trD - tr
        sc = strongly_connected(A)
        key = (tr, d, sc)
        hist6[key] += 1
        if key not in examples:
            examples[key] = (bits, scores(A))
        v1 += (trD < tr + 1)
        v2 += (trD > 2 * tr)
    w(f'L1 violations: {v1}   L2 violations: {v2}   ({time.time()-t0:.1f}s)')
    w(f'{"trans":>6} {"delta":>6} {"strong":>7} {"count":>7}   example (upper-tri bits, scores)')
    for key in sorted(hist6):
        tr, d, sc = key
        bits, scs = examples[key]
        w(f'{tr:>6} {d:>6} {"Y" if sc else "N":>7} {hist6[key]:>7}   bits={bits} scores={scs}')
    w('')
    sc_deltas = sorted({d for (tr, d, sc) in hist6 if sc})
    w(f'deltas over STRONGLY CONNECTED n=6 tournaments: {sc_deltas}')
    w(f'deltas over all n=6: {sorted({d for (tr, d, sc) in hist6})}')
    w('')

    # ---------- PART 3: specific large tournaments ----------
    w('--- PART 3: specific tournaments ---')
    cores = build_tower()
    specials = []
    specials.append(('QR7 (=ST7, unique TT4-free)', paley(7)))
    try:
        ST13 = np.load(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\st13_adj_kps2.npy')
        specials.append(('ST13 (unique TT5-free, circulant from script A)', ST13))
    except FileNotFoundError:
        w('ST13 adjacency not found (script A may still be running) — rebuilding via Z13 scan')
        found = None
        pairs13 = [(d, 13 - d) for d in range(1, 7)]
        for b in range(64):
            Sset = {a if (b >> k) & 1 else bb for k, (a, bb) in enumerate(pairs13)}
            A = circulant(13, Sset)
            if trans_of(A)[0] == 4:
                found = A
                break
        if found is not None:
            specials.append(('ST13 (unique TT5-free, circulant)', found))
        else:
            w('  no TT5-free circulant on Z13 found — ST13 not circulant, skipped')
    specials.append(('Galois_27 (=ST27, unique TT6-free)', galois_27()))
    specials.append(('T7 (tower core)', cores[7]))
    specials.append(('T15 (tower core)', cores[15]))
    specials.append(('T31 (tower core)', cores[31]))

    w(f'{"name":>45} {"n":>4} {"trans":>6} {"transD":>7} {"delta":>6} {"2n":>4}')
    for name, A in specials:
        t0 = time.time()
        tr, _ = trans_of(A)
        Ad, _ = D_skew(A)
        trD, sv = trans_of(Ad)
        w(f'{name:>45} {A.shape[0]:>4} {tr:>6} {trD:>7} {trD-tr:>6} {2*A.shape[0]:>4}'
          f'   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
    w('')

    # T63 double: 126 vertices. Attempt exact trans with node budget via bracket:
    w('--- D(T63) on 126 vertices (bracketed) ---')
    T63 = cores[63]
    tr63, _ = trans_of(T63)
    Ad, _ = D_skew(T63)
    mk, n126 = masks_of(Ad)
    sv = TransSolver(mk)
    full = (1 << n126) - 1
    t0 = time.time()
    # lower bound is tr63 + 1 by L1; probe upward with decision queries
    k = tr63 + 1
    last_yes = None
    try:
        while True:
            res = sv.has_TT(full, k)
            w(f'has_TT(D(T63), {k}) = {res}   (nodes={sv.nodes}, {time.time()-t0:.1f}s)')
            if res:
                last_yes = k
                k += 1
            else:
                break
            if time.time() - t0 > 1200:
                w('  budget reached; stopping probe')
                break
    except (MemoryError, RecursionError) as e:
        w(f'  aborted: {type(e).__name__}')
    if last_yes is not None:
        w(f'trans(D(T63)) = {last_yes} exactly' if not res else
          f'trans(D(T63)) >= {last_yes} (upper probe incomplete)')
    w('')

    # ---------- PART 4: contrast with the tower step ----------
    w('--- PART 4: D(T_m) (2m vertices) vs tower step T_{2m+1} ---')
    w('tower step n -> 2n+1 is the BORDERED double (skew-Hadamard order doubling);')
    w('D(T) is the unbordered skew double on 2n vertices. Compare trans growth:')
    for m, nxt in ((3, 7), (7, 15), (15, 31), (31, 63)):
        trm, _ = trans_of(cores[m])
        trnxt, _ = trans_of(cores[nxt])
        AD, _ = D_skew(cores[m])
        trD, _ = trans_of(AD)
        w(f'  T{m}: trans={trm}   D(T{m}) [n={2*m}]: trans={trD} (delta {trD-trm})   '
          f'T{nxt} [n={nxt}]: trans={trnxt} (delta {trnxt-trm})')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
