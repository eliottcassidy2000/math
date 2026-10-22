"""Auditor: hostile structured families for HYP-9370 (a)-(c) (not in the lane's MCMC samples).
- all circulant regular q-tournaments (q = 13, 15, 17) and their one-vertex deletions (almost-regular n = q - 1);
- R[K_2] for every regular R on 5, 7 (3 classes) and a few on 9 (almost-regular n = 10, 14, 18);
- lexicographic products C_3[C_5], C_5[C_3] (regular 15) and their vertex deletions (n = 14).
Usage: python zeros_structured.py > zeros_structured.out"""
import os, sys, itertools
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn

HERE = os.path.dirname(os.path.abspath(__file__))
K2 = np.array([[0, 1], [0, 0]])


def circulants(q):
    half = list(range(1, (q - 1) // 2 + 1))
    seen = set(); out = []
    for signs in itertools.product((1, -1), repeat=len(half)):
        conn = sorted((s * x) % q for s, x in zip(signs, half))
        # multiplier equivalence is not used: keep every connection set (isomorphic duplicates are harmless)
        key = tuple(conn)
        if key in seen:
            continue
        seen.add(key); out.append(conn)
    return out


def t0_of(As, label):
    by_n = {}
    for i, A in enumerate(As):
        by_n.setdefault(A.shape[0], []).append(i)
    res = [None] * len(As)
    for n, idx in by_n.items():
        codes = [tourn.to_code(As[i]) for i in idx]
        nch = min(12, len(codes))
        chunks = [codes[k::nch] for k in range(nch)]
        with ThreadPoolExecutor(nch) as ex:
            outs = list(ex.map(lambda ch: tourn.forest_numbers(ch, n), chunks))
        flat = [None] * len(codes)
        for k in range(nch):
            for j, r in enumerate(outs[k]):
                flat[k + j * nch] = r
        for i, (H, c, _) in zip(idx, flat):
            N = tourn.N_from_c(c, n)
            res[i] = tourn.smallest_zero(N, n) + (H,)
    t = np.array([r[0] for r in res]); kinds = [r[2] for r in res]
    i = int(np.argmin(t))
    print(f'{label}: {len(As)} tournaments, t0 min {t.min():.4f} median {np.median(t):.4f} max {t.max():.4f}; '
          f'kinds {dict((k, kinds.count(k)) for k in set(kinds))}; argmin #{i} H={res[i][3]} w={res[i][1]}', flush=True)
    return res


def main():
    for q in (13, 15, 17):
        cs = circulants(q)
        regs = [tourn.circulant(q, c) for c in cs]
        t0_of(regs, f'circulant regular q={q} (all {len(cs)} connection sets)')
        dels = [tourn.delete(R, 0) for R in regs]
        t0_of(dels, f'circulant q={q} minus a vertex (almost-regular n={q - 1})')
    R5 = tourn.circulant(5, [1, 2])
    reg7 = [tourn.from_code(c) for c in ('000111010010010100010', '000111000110001000000', '000111000110100000010')]
    t0_of([tourn.lex(R5, K2)], 'C_5[K_2] (almost-regular 10)')
    t0_of([tourn.lex(R, K2) for R in reg7], 'R[K_2], R the 3 regular 7-tournaments (almost-regular 14)')
    C3 = tourn.circulant(3, [1])
    t0_of([tourn.lex(C3, R5), tourn.lex(R5, C3)], 'C_3[C_5], C_5[C_3] (regular 15)')
    t0_of([tourn.delete(tourn.lex(C3, R5), 0), tourn.delete(tourn.lex(R5, C3), 0)], 'C_3[C_5], C_5[C_3] minus a vertex (n=14)')
    t0_of([tourn.lex(C3, tourn.lex(C3, K2))], 'C_3[C_3[K_2]] (almost-regular 18)')
    t0_of([tourn.lex(tourn.circulant(9, [1, 2, 3, 4]), K2), tourn.lex(tourn.lex(C3, C3), K2)],
          'C_9(1,2,3,4)[K_2], C_3[C_3][K_2] (almost-regular 18)')


if __name__ == '__main__':
    main()
