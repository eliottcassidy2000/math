#!/usr/bin/env python3
"""Triangular numbers as the add/mult bridge; +2 ladder = odd gnomons; x2 doubling;
the worry-modulus identity 8*C(n,2)+1=(2n-1)^2. opus-2026-06-03-S586."""
def T(k): return k*(k+1)//2
def C2(n): return n*(n-1)//2
if __name__=='__main__':
    for n in range(3,16):
        print(f"n={n}: C(n,2)=T_{{n-1}}={C2(n)}; 8*C(n,2)+1={8*C2(n)+1}=(2n-1)^2={(2*n-1)**2}; "
              f"+2 gnomon(2n+1)={2*n+1 if n<14 else '-'}")
