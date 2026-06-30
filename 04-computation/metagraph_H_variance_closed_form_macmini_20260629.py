#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THM-589: the metagraph H-variance closed form Var(H)=(n!/4^{n-1})(W(n)-n!), W(n)=THM-219's
no-unit-descent succession count = the owner's odd-composition sum, OGF kernel x(1+x^2)/(1-x^2);
exact rate CV(H)^2 ~ 2/n; Poisson-Euler = Poisson(1). (mac-mini-2026-06-29-S22)
"""
from __future__ import annotations
import functools, itertools, math
print = functools.partial(print, flush=True)


def W_succession(n):
    """W(n) = sum_{perm: no unit descent} 2^{#unit ascents} (THM-219)."""
    tot = 0
    for p in itertools.permutations(range(1, n + 1)):
        asc = 0; ok = True
        for i in range(n - 1):
            if p[i + 1] == p[i] + 1: asc += 1
            elif p[i + 1] == p[i] - 1: ok = False; break
        if ok: tot += 1 << asc
    return tot


def W_oddcomp(n):
    """W(n) = sum_{compositions of n into odd parts} k! * 2^{#parts>=3} (owner's form)."""
    def comps(m):
        if m == 0: yield (); return
        for first in range(1, m + 1, 2):
            for rest in comps(m - first): yield (first,) + rest
    return sum(math.factorial(len(c)) * (1 << sum(1 for p in c if p >= 3)) for c in comps(n))


def W_ogf(n):
    """W(n) = sum_k k! [x^n] W(x)^k, W(x)=x(1+x^2)/(1-x^2)=x+2x^3+2x^5+..."""
    w = [0] * (n + 1)
    for p in range(1, n + 1):
        if p % 2 == 1: w[p] = 1 if p == 1 else 2
    Wk = [0] * (n + 1); Wk[0] = 1; tot = 0
    for k in range(1, n + 1):
        new = [0] * (n + 1)
        for a in range(n + 1):
            if Wk[a]:
                for p in range(1, n + 1 - a):
                    if w[p]: new[a + p] += Wk[a] * w[p]
        Wk = new; tot += math.factorial(k) * Wk[n]
    return tot


def var_H_brute(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    perms = list(itertools.permutations(range(n)))
    s = 0; s2 = 0; c = 0
    for bits in range(1 << len(pairs)):
        arc = [[False] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1: arc[i][j] = True
            else: arc[j][i] = True
        H = sum(1 for pm in perms if all(arc[pm[k]][pm[k + 1]] for k in range(n - 1)))
        s += H; s2 += H * H; c += 1
    from fractions import Fraction as F
    EH = F(s, c); return EH, F(s2, c) - EH * EH


def main():
    print("THM-589: metagraph H-variance closed form (mac-mini-S22)")
    print("=" * 72)
    print(f"{'n':>2} {'W_succ':>7} {'W_comp':>7} {'W_ogf':>7} {'Var(H)closed':>13} "
          f"{'Var(H)brute':>13} {'CV^2':>8} {'n*CV^2':>7}")
    for n in range(3, 8):
        ws, wc, wo = W_succession(n), W_oddcomp(n), W_ogf(n)
        nf = math.factorial(n)
        var_closed = nf * (wc - nf) / 4**(n - 1)
        cv2 = wc / nf - 1
        chk = ""
        if n <= 6:
            EH, vb = var_H_brute(n)
            chk = f"{float(vb):.4f}" + ("" if abs(float(vb) - var_closed) < 1e-9 else " MISMATCH")
        agree = "OK" if ws == wc == wo else "MISMATCH"
        print(f"{n:>2} {ws:>7} {wc:>7} {wo:>7} {var_closed:>13.4f} {chk:>13} "
              f"{cv2:>8.4f} {n*cv2:>7.4f}  [{agree}]")
    print("\n  W(n) = 1,2,8,32,158,928,6350,49752 (=THM-219 W, =S90/S112 simplicial-Redei).")
    print("  Var(H) closed = brute (n<=6). CV(H)^2 = W(n)/n! - 1 ~ 2/n (n*CV^2 -> 2).")
    print("  Poisson-Euler (THM-219) = Poisson(1) (HYP-3560): NUD/n!->1/e, E[2^asc]->e, product->1.")


if __name__ == "__main__":
    main()
