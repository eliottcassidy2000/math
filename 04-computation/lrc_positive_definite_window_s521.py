#!/usr/bin/env python3
"""
lrc_positive_definite_window_s521.py   claudebox-2026-06-01-S521

"Global spread => local emptiness": the positive-definite-window obstruction
(reflection: 07-reflections/lrc-global-spread-local-emptiness-s521.md).

If a window w>=0 supported in the safe band B=[1/n,1-1/n] had ALL Fourier coeffs
hhat(c)>=0 (positive-definite = 'globally spread'), then
  int prod_i w(v_i t) dt = sum_{<c,v>=0} prod hhat(c_i) >= hhat(0)^m > 0  => LRC.
IMPOSSIBLE: B excludes the observer 0, so w(0)=0=sum_c hhat(c); hhat>=0 forces w==0.
Demonstrated: the band window is sign-indefinite (many negative coeffs). The signs
are why mu can cancel to 0 (tight/resonant cases). Spread that DOES work =
relation-sparsity (proved sufficient condition).
"""
from math import sin, pi

def hhat_band(c, n):
    if c == 0: return 1 - 2/n
    return ((-1)**c) * sin(pi*c*(1-2/n)) / (pi*c)

def main():
    print("Band-window Fourier coefficients hhat(c): positive-definite (all >=0)?  [the spread test]")
    for n in (5, 7):
        negs = [c for c in range(1, 20) if hhat_band(c, n) < -1e-12]
        print(f"  n={n}: negative coeffs among c=1..19: {negs}")
        print(f"        sample: hhat(1)={hhat_band(1,n):+.4f}, hhat(2)={hhat_band(2,n):+.4f}, "
              f"hhat(0)={hhat_band(0,n):.4f}")
    print("\n  => SIGN-INDEFINITE. The band excludes the observer 0, so any w>=0 on B has")
    print("     w(0)=sum_c hhat(c)=0; positive-definiteness (hhat>=0) then forces w==0.")
    print("     No positive-definite window on the safe band exists -> no pure-spread proof of LRC.")
    print("     The permitted sign-cancellation IS the resonance obstruction (mu can be 0 on tight sets).")
    print("\n  Working spread condition (PROVED elsewhere): relation-sparsity")
    print("     sum_{c!=0,<c,v>=0} |prod hhat(c_i)| < (1-2/n)^m  =>  mu>0  =>  lonely.")

if __name__ == "__main__":
    main()
