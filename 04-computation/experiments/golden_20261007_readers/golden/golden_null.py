#!/usr/bin/env python3
"""Null model for Fibonacci/Lucas hits on the 27-trajectory (the trunk containing 332, 233, 425, 223's tail).
Fibonacci and Lucas are both round(a*phi^k) sequences; compare with random a in [1, phi)."""
import random
phi = (1 + 5 ** 0.5) / 2
def std_orbit(n):
    o = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        o.append(n)
    return o
traj = set(std_orbit(27))
LO, HI = 29, 10 ** 4
def geo(a):
    s = set(); k = 0
    while True:
        v = round(a * phi ** k); k += 1
        if v > HI: break
        if v >= LO: s.add(v)
    return s
fibset = geo(1 / 5 ** 0.5); lucset = geo(1.0)
obs = len(traj & (fibset | lucset))
print('Fibonacci values in [29,1e4]:', sorted(fibset)); print('Lucas values in [29,1e4]:', sorted(lucset))
print('observed Fib/Lucas hits on 27-trajectory:', sorted(traj & (fibset | lucset)), '=', obs)
random.seed(1)
M = 20000; ge = 0; tot = 0
for _ in range(M):
    a1 = random.uniform(1, phi); a2 = random.uniform(1, phi)
    h = len(traj & (geo(a1) | geo(a2)))
    tot += h; ge += (h >= obs)
print(f'null: two random ratio-phi sequences: mean hits {tot/M:.2f}, P(hits >= {obs}) = {ge/M:.4f}')
# single trajectory vs population: over trajectories of n in [1,1000], how often >= 4 Fib/Lucas hits >= 29
cnt = sum(1 for n in range(1, 1001) if len(set(std_orbit(n)) & (fibset | lucset)) >= obs)
print(f'n <= 1000 with >= {obs} Fib/Lucas values >= 29 on the orbit: {cnt}')
# how many distinct 'trunks' carry them: orbits with >= obs hits that do NOT pass 9232
cnt2 = sum(1 for n in range(1, 1001) if len(set(std_orbit(n)) & (fibset | lucset)) >= obs and 9232 not in set(std_orbit(n)))
print('   of which not through 9232:', cnt2)
