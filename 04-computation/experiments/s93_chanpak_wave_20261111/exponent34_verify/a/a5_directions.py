#!/usr/bin/env python3
"""Independent verifier, side checks on how (A) is used/quoted in the remarks and receipts of THM-5141 (draft).

R1  Remark 1 parenthetical: "at most C 6^(w/2) pairs for each direction".  For every primitive direction (k, k')
    (identified with its negative) the number of ordered pairs x != x' with (X, Xi) in Z (k, k') is computed from the
    signed strings (weight 2^z each), for all pairs and for odd pairs; the maximum is compared with 6^(w/2).
R2  Remark 2: pairs of height 1 (X = +-Xi).  Corollary of (A)1: D = D* iff eps is a palindrome, D = -D* iff eps is an
    antipalindrome.  Checked by enumeration, with the closed forms
        odd pairs, even w: 2 * 2 * 6^(w/2 - 1) = (2/3) 6^(w/2);   odd w: 12 * 6^((w-3)/2) = 0.8165 * 6^(w/2);
        all pairs, direction (1,1): palindromes minus the zero string: even w 6^(w/2) - 2^w, odd w 4 * 6^((w-1)/2) - 2^w.
Usage: python a5_directions.py [WMAX (default 14)]
"""
import sys
import time
import numpy as np


def strings(w):
    n = 3 ** w
    t = np.arange(n, dtype=np.int64)
    D = np.zeros(n, np.int64)
    Ds = np.zeros(n, np.int64)
    z = np.zeros(n, np.int64)
    pal = np.ones(n, bool)
    apal = np.ones(n, bool)
    digs = []
    for i in range(w):
        e = (t % 3) - 1
        t = t // 3
        digs.append(e.astype(np.int8))
        D += e << i
        Ds += e << (w - 1 - i)
        z += (e == 0)
    for i in range(w):
        pal &= (digs[i] == digs[w - 1 - i])
        apal &= (digs[i] == -digs[w - 1 - i])
    return D, Ds, z, pal, apal


def main():
    WMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    ok = True
    for w in range(2, WMAX + 1):
        t0 = time.time()
        D, Ds, z, pal, apal = strings(w)
        nzmask = (D != 0) | (Ds != 0)
        wt = (np.int64(1) << z)
        # R2: D = D* iff palindrome; D = -D* iff antipalindrome
        c_pal = bool(np.all((D == Ds) == pal)) and bool(np.all((D == -Ds) == apal))
        odd = ((D & 1) == 1) & ((Ds & 1) == 1)
        h1_odd = int(wt[odd & ((D == Ds) | (D == -Ds))].sum())
        if w % 2 == 0:
            pred_h1 = 4 * 6 ** (w // 2 - 1)
            pred_11 = 6 ** (w // 2) - 2 ** w
        else:
            pred_h1 = 12 * 6 ** ((w - 3) // 2) if w >= 3 else None
            pred_11 = 4 * 6 ** ((w - 1) // 2) - 2 ** w
        n11 = int(wt[nzmask & (D == Ds)].sum())
        c_h1 = (pred_h1 is None or h1_odd == pred_h1) and (n11 == pred_11)
        # R1: direction counts
        Dn = D[nzmask]
        Dsn = Ds[nzmask]
        wn = wt[nzmask]
        g = np.gcd(np.abs(Dn), np.abs(Dsn))
        k = Dn // g
        kp = Dsn // g
        sgn = np.where((k < 0) | ((k == 0) & (kp < 0)), -1, 1)
        k = k * sgn
        kp = kp * sgn
        size = 1 << w
        key = k * (4 * size) + (kp + size)
        uk, inv = np.unique(key, return_inverse=True)
        cnt = np.bincount(inv, weights=wn.astype(np.float64)).astype(np.int64)
        im = int(np.argmax(cnt))
        kb = int(uk[im] // (4 * size))
        kpb = int(uk[im] % (4 * size)) - size
        mx = int(cnt[im])
        oddn = odd[nzmask]
        uko, invo = np.unique(key[oddn], return_inverse=True)
        cnto = np.bincount(invo, weights=wn[oddn].astype(np.float64)).astype(np.int64)
        imo = int(np.argmax(cnto))
        kbo = int(uko[imo] // (4 * size))
        kpbo = int(uko[imo] % (4 * size)) - size
        mxo = int(cnto[imo])
        six = 6 ** (w / 2)
        s1 = (w + 1) // 2
        s2 = w // 2
        bound_odd = 6 ** s1 + 2 * 6 ** s2          # (A)3 + (D1): (1/2)(2 * 6^s1 + 4 * 6^s2)
        c_b = (mxo <= bound_odd)
        ok &= c_pal and c_h1 and c_b
        print(f"  w = {w:2d}: D = D* iff palindrome, D = -D* iff antipalindrome: {c_pal}; height-1 odd pairs = {h1_odd} "
              f"(closed form: {c_h1}; = {h1_odd / six:.4f} * 6^(w/2)); all pairs: max direction ({kb},{kpb}) with {mx} pairs = {mx / six:.4f} * 6^(w/2); "
              f"odd pairs: max direction ({kbo},{kpbo}) with {mxo} = {mxo / six:.4f} * 6^(w/2) (<= 6^s1 + 2*6^s2 = {bound_odd / six:.3f} * 6^(w/2): {c_b})   [{time.time()-t0:.1f}s]", flush=True)
    print("  limits: all pairs, odd w, direction (1,1): 4/sqrt(6) = %.4f;  even w, direction (1,2): 4/3 = %.4f" % (4 / 6 ** 0.5, 4 / 3))
    print("ALL CHECKS PASSED" if ok else "*** SOME CHECK FAILED ***")


if __name__ == "__main__":
    main()
