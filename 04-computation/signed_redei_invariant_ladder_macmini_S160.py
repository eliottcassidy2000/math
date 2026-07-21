#!/usr/bin/env python3
"""
THM-1966: where the signed Redei count |R| sits on THM-1780's moment/co-spectral ladder.
  n<=5: |R| SPECTRAL (constant on co-spectral classes).
  n=6 : |R| leaves the ladder WITH H (same 3 co-spectral splits, perfectly coupled), but
        |R| = f(spectrum,H) -- adds nothing beyond spectrum+H.  (H,|R|) distinguishes 31/56
        iso classes vs H's 19 and the spectrum's 28.
  n=7 : |R| DECOUPLES -- co-spectral tournaments with SAME H but DIFFERENT |R| exist.
        Witness: two non-iso 7-tournaments, char poly x^7-9x^4-12x^3-16x^2-8x-1, both H=81,
        |R| = 1 vs 17.  => |R| is a genuinely independent invariant from n=7.
  max|R| at n=7 = 147 = QR(7) Paley (regular) tournament.
                                                                    mac-mini-2026-07-21-S160
"""
import numpy as np, random
from itertools import combinations, product, permutations
from collections import defaultdict

def hR(n, A):
    dph = [[0]*n for _ in range(1 << n)]; dpr = [[0]*n for _ in range(1 << n)]
    for v in range(n): dph[1 << v][v] = 1; dpr[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            ch, cr = dph[mask][last], dpr[mask][last]
            if not ch: continue
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if A[last][nxt]:
                    dph[mask | (1 << nxt)][nxt] += ch
                    invs = bin(mask >> (nxt+1)).count("1")
                    dpr[mask | (1 << nxt)][nxt] += cr * (-1 if invs & 1 else 1)
    full = (1 << n)-1
    return sum(dph[full][v] for v in range(n)), abs(sum(dpr[full][v] for v in range(n)))

def moments(n, A):
    Anp = np.array(A, dtype=object); tr = []; Cur = Anp.copy()
    for k in range(1, n+1): tr.append(int(np.trace(Cur))); Cur = Cur @ Anp
    return tuple(tr)

def all_tours(n):
    pairs = list(combinations(range(n), 2))
    for bits in product((0, 1), repeat=len(pairs)):
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def canon(n, A):
    best = None
    An = np.array(A)
    for p in permutations(range(n)):
        code = tuple(An[np.ix_(list(p), list(p))].flatten())
        if best is None or code < best: best = code
    return best

if __name__ == "__main__":
    # ladder placement: |R|-split of co-spectral classes (fast spectral pass, exhaustive n<=6).
    # iso-refinement counts at n=6 (from the full-canon pass, R_refine2.out) are documented constants
    # to keep this script fast; rerun with `iso[canon(n,A)]=...` for n=6 to reproduce them (~2 min).
    ISO_N6 = dict(H=19, R=8, HR=31, spectrum=28, spec_H=32, spec_H_R=32)  # verified exhaustively
    for n in (4, 5, 6):
        by_spec = defaultdict(lambda: {'H': set(), 'R': set()})
        iso = {}
        do_canon = (n <= 5)
        for A in all_tours(n):
            h, r = hR(n, A); sp = moments(n, A)
            by_spec[sp]['H'].add(h); by_spec[sp]['R'].add(r)
            if do_canon: iso[canon(n, A)] = (h, r, sp)
        Rsplit = sum(1 for d in by_spec.values() if len(d['R']) > 1)
        if do_canon:
            isos = list(iso.values()); nd = lambda f: len({f(v) for v in isos})
            extra = (f"iso by H={nd(lambda v:v[0])}, |R|={nd(lambda v:v[1])}, "
                     f"(H,|R|)={nd(lambda v:(v[0],v[1]))}, spectrum={nd(lambda v:v[2])}, "
                     f"(spec,H,|R|)={nd(lambda v:(v[2],v[0],v[1]))}")
        else:
            extra = (f"iso by H={ISO_N6['H']}, |R|={ISO_N6['R']}, (H,|R|)={ISO_N6['HR']}, "
                     f"spectrum={ISO_N6['spectrum']}, (spec,H)={ISO_N6['spec_H']}, "
                     f"(spec,H,|R|)={ISO_N6['spec_H_R']}  [documented; |R|=f(spec,H) at n=6]")
        print(f"n={n}: |R| splits {Rsplit} co-spectral classes; {extra}", flush=True)

    # n=7 independence witness
    def rand_tour(n, rng):
        A = [[0]*n for _ in range(n)]
        for i, j in combinations(range(n), 2):
            A[i][j], A[j][i] = (1, 0) if rng.random() < 0.5 else (0, 1)
        return A
    rng = random.Random(77); n = 7; buckets = defaultdict(list); found = None
    for _ in range(200000):
        A = rand_tour(n, rng); h, r = hR(n, A); key = (moments(n, A), h)
        if any(rr != r for rr, _ in buckets[key]):
            A0 = [x for x in buckets[key] if x[0] != r][0]; found = (A0[1], A0[0], A, r, key); break
        buckets[key].append((r, A))
    if found:
        A0, r0, A1, r1, (sp, h) = found
        print(f"n=7 WITNESS: co-spectral, H={h}, |R|={r0} vs {r1}, non-iso={canon(7,A0)!=canon(7,A1)} "
              f"=> |R| independent of (spectrum,H) from n=7", flush=True)
