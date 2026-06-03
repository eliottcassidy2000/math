#!/usr/bin/env python3
"""claudebox-S611: LRC(14) reduction. (1) unit-clock sieve: no multiple of 14 => gap>=1/14 (lonely).
(2) division sieve m in 2..14: no multiple of m => gap>=1/m>=1/14. Residual = covers every m in
2..14 (a covering system). See HYP-2170."""
import random
def caught(V):
    for m in range(2,15):
        if not any(v%m==0 for v in V): return m
    return None
if __name__=="__main__":
    random.seed(14); N=200000; esc=0
    for _ in range(N):
        V=random.sample(range(1,80),13)
        if caught(V) is None: esc+=1
    print(f'division sieve m in 2..14: caught {100*(N-esc)/N:.2f}%, residual (covering) {100*esc/N:.2f}%')
