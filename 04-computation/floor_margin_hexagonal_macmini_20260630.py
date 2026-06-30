#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The LRC floor margin is 1/(hexagonal number). mac-mini-2026-06-30-S48.
margin(n) = 2/(2n-1) - 1/n = 1/(n(2n-1)) = 1/H_n = 1/T_(2n-1) = 1/C(2n,2) = 1/dim so(2n) = 2B(2n-1,2).
Sum_n = 2 ln 2 = ln 4. Doubling bridge: Phi_6(2n) = 2*[denom margin(n)] + 1.
"""
from __future__ import annotations
import functools, math
print = functools.partial(print, flush=True)


def main():
    print("margin(n) = 2/(2n-1) - 1/n = 1/(n(2n-1)); the denominator's many faces:\n")
    print(f"  {'n':>3} {'n(2n-1)':>8} {'H_n':>6} {'T_(2n-1)':>9} {'C(2n,2)':>8} {'dim so(2n)':>11} {'Phi6(2n)=2d+1':>14}")
    ok = True
    for n in range(2, 16):
        d = n*(2*n-1)
        H = n*(2*n-1)
        T = (2*n-1)*(2*n)//2
        C = math.comb(2*n, 2)
        so = n*(2*n-1)
        P2 = (2*n)**2-(2*n)+1
        ok &= (H == d == T == C == so) and (P2 == 2*d+1)
        if n in (4, 7, 8, 14):
            print(f"  {n:>3} {d:>8} {H:>6} {T:>9} {C:>8} {so:>11} {P2:>14}")
    print(f"\n  all-equal & doubling-bridge hold for n=2..15: {ok}")

    S = sum(1.0/(n*(2*n-1)) for n in range(1, 4_000_000))
    print(f"\n  Sum_n 1/(n(2n-1)) = {S:.7f}  vs  2 ln 2 = {2*math.log(2):.7f} = ln 4 = {math.log(4):.7f}")

    print("\n  Beta-moment: margin = 2*B(2n-1,2) = 2*int_0^1 x^(2n-2)(1-x) dx:")
    for n in (7, 14):
        beta = 2*math.gamma(2*n-1)*math.gamma(2)/math.gamma(2*n+1)
        print(f"    n={n}: 1/(n(2n-1))={1/(n*(2*n-1)):.8f}  2B(2n-1,2)={beta:.8f}  match={abs(beta-1/(n*(2*n-1)))<1e-12}")

    print("\n  LRC14 hinge: T_13 = 91 = margin-denom(7) = 7*13 = (Phi6(14)-1)/2;  14 = 2*7 (apex doubled).")
    print(f"    7*13={7*13}, (Phi6(14)-1)/2={(14*14-14+1-1)//2}, T_13={13*14//2}")


if __name__ == "__main__":
    main()
