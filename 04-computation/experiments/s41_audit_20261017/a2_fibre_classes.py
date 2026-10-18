#!/usr/bin/env python3
"""S41 audit A2 (independent referee code): THM-4667 statement 8 and its proof's class description, by brute force on
every primitive seed below LIM.

Checks
  * a(s) = nu_2(3s+1) is 2 exactly on s = 1 (8) and 1 on s = 3, 7 (8);
  * heads s = 1 (8): 2s+1 = 3 (16), a(2s+1) = 1, P(2s+1) = P(s);
  * the residue claims of the proof: x = 8k+3 gives (3x+1)/2 = 4k+5 (8); x = 7 (8) gives (3x+1)/2 = 3 or 7 (8);
    x = 11 (16) gives (3x+1)/2 = 1 (8); hence singletons x = 7 (8) or 11 (16) have P(x) = (3x+1)/2, x = (2P(x)-1)/3,
    P(x) = 2 (3); and the three classes partition the seeds;
  * every fibre P^(-1)(u), restricted to seeds below LIM, equals the explicit formula
    {(2^(a_j) B^j(u) - 1)/3 : B^j(u) != 0 (3)} below LIM, and has the shape [singleton (2u-1)/3 iff u = 2 (3)] followed by
    pairs (s_k, 2 s_k + 1) with s_(k+1) = 64 s_k + 49.
Usage: python a2_fibre_classes.py [LIM]
"""
import sys, time, json, os
sys.dont_write_bytecode = True
import numpy as np


def seed_vec(u):
    u = u.copy()
    m = (u % 8) == 5
    while m.any():
        u[m] = (u[m] - 1) // 4
        m = (u % 8) == 5
    return u


def main():
    LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 4_000_000
    t0 = time.time()
    odd = np.arange(1, LIM, 2, dtype=np.int64)
    S = odd[odd % 8 != 5]
    y = 3 * S + 1
    low = y & (-y)
    a = np.log2(low).round().astype(np.int64)
    assert np.all((a == 1) | (a == 2)), 'a in {1, 2} on seeds'
    assert np.array_equal(a == 2, S % 8 == 1), 'a = 2 exactly on s = 1 (8)'
    img_raw = y >> a
    Pimg = seed_vec(img_raw)

    # map seed -> position via searchsorted
    def pos(x):
        k = np.searchsorted(S, x)
        ok = (k < len(S))
        ok[ok] &= S[k[ok]] == x[ok]
        return k, ok

    # heads and tails
    H = S[S % 8 == 1]
    T = 2 * H + 1
    assert np.all(T % 16 == 3), 'tails are 3 (16)'
    k, ok = pos(T)
    kh, okh = pos(H)
    assert okh.all()
    assert np.all(a[k[ok]] == 1), 'a(2s+1) = 1'
    assert np.array_equal(Pimg[k[ok]], Pimg[kh[ok]]), 'P(2s+1) = P(s)'
    n_pairs = int(ok.sum())
    # residue claims
    x3 = S[S % 8 == 3]
    kk = (x3 - 3) // 8
    assert np.all(((3 * x3 + 1) // 2) % 8 == (4 * kk + 5) % 8), '(3x+1)/2 = 4k+5 (8) for x = 8k+3'
    x7 = S[S % 8 == 7]
    r7 = ((3 * x7 + 1) // 2) % 8
    assert np.all((r7 == 3) | (r7 == 7)), 'x = 7 (8): (3x+1)/2 = 3 or 7 (8)'
    x11 = S[S % 16 == 11]
    assert np.all(((3 * x11 + 1) // 2) % 8 == 1), 'x = 11 (16): (3x+1)/2 = 1 (8)'
    single_mask = (S % 8 == 7) | (S % 16 == 11)
    xs = S[single_mask]
    assert np.all(a[single_mask] == 1)
    assert np.array_equal(Pimg[single_mask], (3 * xs + 1) // 2), 'singletons: P(x) = (3x+1)/2 (no braid stripped)'
    assert np.all((2 * Pimg[single_mask] - 1) == 3 * xs) and np.all(Pimg[single_mask] % 3 == 2)
    head_mask = S % 8 == 1
    tail_mask = S % 16 == 3
    assert np.all(head_mask.astype(int) + tail_mask.astype(int) + single_mask.astype(int) == 1), 'partition of the seeds'
    # tails are never singletons in their fibre: their image needs a stripped braid (j >= 1)
    assert np.all(img_raw[tail_mask] != Pimg[tail_mask]), 'tails have j >= 1'

    # fibres: group by image
    order = np.lexsort((S, Pimg))
    Ps = Pimg[order]
    Ss = S[order]
    starts = np.flatnonzero(np.r_[True, Ps[1:] != Ps[:-1]])
    ends = np.r_[starts[1:], len(Ps)]
    nfib = 0
    nsingle = 0
    maxlen = 0
    for st, en in zip(starts, ends):
        u = int(Ps[st])
        elems = [int(v) for v in Ss[st:en]]
        # explicit formula below LIM
        form = []
        Bj = u
        j = 0
        while True:
            r = Bj % 3
            if r != 0:
                aj = 2 if r == 1 else 1
                s = ((1 << aj) * Bj - 1) // 3
                if s < LIM:
                    form.append(s)
            # stop when even the smallest possible next element exceeds LIM
            if (2 * Bj - 1) // 3 >= LIM:
                break
            Bj = 4 * Bj + 1
            j += 1
        assert sorted(form) == elems, ('fibre formula', u, form[:6], elems[:6])
        # shape
        rest = elems
        if u % 3 == 2:
            assert rest[0] == (2 * u - 1) // 3, ('singleton', u)
            nsingle += 1
            rest = rest[1:]
        else:
            assert all(e != (2 * u - 1) / 3 for e in rest)
        heads = rest[0::2]
        tails = rest[1::2]
        for i, h in enumerate(heads):
            assert h % 8 == 1, ('head residue', u, h)
            if i < len(tails):
                assert tails[i] == 2 * h + 1, ('tail', u, h)
            if i + 1 < len(heads):
                assert heads[i + 1] == 64 * h + 49, ('spacing 64', u, h, heads[i + 1])
        nfib += 1
        maxlen = max(maxlen, len(elems))
    res = dict(LIM=LIM, seeds=int(len(S)), heads=int(head_mask.sum()), pairs_checked=n_pairs,
               singletons=int(single_mask.sum()), fibres=nfib, fibres_with_singleton=nsingle, max_fibre_len=maxlen)
    print('A2: all checks passed below LIM = %d: %s [%.0f s]' % (LIM, res, time.time() - t0))
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a2_fibre_classes.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
