#!/usr/bin/env python3
"""Benchmark n=8 OCF verification speed."""
import time
from itertools import combinations

N = 8
FULL = (1 << N) - 1
ARC_PAIRS = []
for a in range(N):
    for b in range(a+1, N):
        if not (a == 0 and b == 1):
            ARC_PAIRS.append((a, b))

def ham_count(T):
    dp = [[0]*N for _ in range(1<<N)]
    for v in range(N): dp[1<<v][v] = 1
    for mask in range(1, 1<<N):
        for last in range(N):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(N):
                if mask & (1<<nxt): continue
                if T[last*N+nxt]: dp[mask|(1<<nxt)][nxt] += c
    return sum(dp[FULL])

def count_dc(T, combo):
    m = len(combo); start = combo[0]; others = list(combo[1:]); mo = len(others)
    if mo == 0: return 0
    ofull = (1<<mo)-1
    dp = [[0]*mo for _ in range(1<<mo)]
    for i,o in enumerate(others):
        if T[start*N+o]: dp[1<<i][i] = 1
    for om in range(1, ofull+1):
        for li in range(mo):
            c = dp[om][li]
            if c==0: continue
            for ni in range(mo):
                if om&(1<<ni): continue
                if T[others[li]*N+others[ni]]: dp[om|(1<<ni)][ni] += c
    return sum(dp[ofull][li] for li in range(mo) if T[others[li]*N+start])

def ocf(T):
    ci = []
    for length in range(3, N+1, 2):
        for combo in combinations(range(N), length):
            vm = 0
            for v in combo: vm |= (1<<v)
            cnt = count_dc(T, combo)
            if cnt > 0: ci.append((vm, cnt))
    tc = sum(c for _,c in ci)
    vd = 0
    nc = len(ci)
    for i in range(nc):
        for j in range(i+1, nc):
            if ci[i][0] & ci[j][0] == 0:
                vd += ci[i][1] * ci[j][1]
    return 1 + 2*tc + 4*vd

t0 = time.time()
fails = 0
count = 4096
for mask in range(count):
    T = [0]*(N*N)
    T[0*N+1] = 1
    for bit,(a,b) in enumerate(ARC_PAIRS):
        if mask&(1<<bit): T[a*N+b]=1
        else: T[b*N+a]=1
    if ham_count(T) != ocf(T):
        fails += 1

elapsed = time.time() - t0
rate = count / elapsed
total = 1 << 27
print(f"Checked {count} in {elapsed:.1f}s ({rate:.0f}/s), {fails} fails")
print(f"Full estimate (1 proc): {total/rate:.0f}s = {total/rate/3600:.1f}hr")
print(f"Full estimate (4 proc): {total/rate/4:.0f}s = {total/rate/4/3600:.1f}hr")
