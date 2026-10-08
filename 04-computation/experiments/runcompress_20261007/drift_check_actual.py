#!/usr/bin/env python3
"""Actual-integer checks of two universal-state transitions (THM-4603 (3)):
 (a) a long post-run ones-run of the SOURCE (X_end = -1 mod 2^L): debt drifts by +1/2 per Terras step (child limit 25/27 -> trivial cycle);
 (b) a long post-run ones-run of the CHILD (Y_end = -1 mod 2^L): after 12 steps the chain is anchored at -1 with debt -5."""
import random
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
rnd = random.Random(5)
def make_source(K, J, r, cls_mod, cls_val):
    # residual source with two-run J, type r, and X_end - 1 = 27(Y_end - 1) with Y_end = cls_val mod cls_mod (cls_mod a power of 2)
    j = 2*J + r
    # u (run-end odd part) : Y_end = 1 + 2^r * 3^(J-3) u  => choose u mod cls_mod accordingly
    Mb = cls_mod << (j + 8)
    target_w = ((cls_val - 1) >> r) % cls_mod          # w' = (Y_end - 1)/2^r mod cls_mod
    u = (pow(3, -(J-3), cls_mod) * target_w) % cls_mod
    Mx = 1 << (j + cls_mod.bit_length() + 8)
    x = (1 + (1 << j) * u) % Mx
    t = (pow(3, -(K-1), Mx) * ((x + 1) // 2)) % Mx       # 2*3^(K-1)*t - 1 = x  (mod)
    if t % 2 == 0: t += Mx
    t += Mx * rnd.getrandbits(30)
    return t
ok_a = ok_b = 0
for trial in range(40):
    K = rnd.randint(8, 30); J = rnd.randint(3, 8); L = 40
    # (b) child long ones-run: Y_end = -1 mod 2^(L+2), r = 1 (Y_end = 3 mod 4)
    r = 1; t = make_source(K, J, r, 1 << (L + 2), (1 << (L + 2)) - 1)
    x = 2*3**(K-1)*t - 1
    if v2(x - 1) != 2*J + r: continue
    y = (x + 1)//27 - 1
    u, v, k = x, y, 3
    for s in range(2*J): k += (u & 1) - (v & 1); u, v = T(u), T(v)
    assert k == 3 and u - 1 == 27*(v - 1)
    assert (v + 1) % (1 << (L - 2)) == 0
    for s in range(12): k += (u & 1) - (v & 1); u, v = T(u), T(v)
    assert k == -5 and 243*(u + 1) == v + 1, (k,)
    ok_b += 1
for trial in range(40):
    K = rnd.randint(8, 30); J = rnd.randint(3, 8); L = 40
    # (a) source long ones-run: X_end = -1 mod 2^(L+2) <=> Y_end = 1 + (X_end - 1)/27 = 25/27 (2-adically) mod 2^(L+2)
    r = 1; M = 1 << (L + 2)
    yend = (1 + (-2) * pow(27, -1, M)) % M
    t = make_source(K, J, r, M, yend)
    x = 2*3**(K-1)*t - 1
    if v2(x - 1) != 2*J + r: continue
    y = (x + 1)//27 - 1
    u, v, k = x, y, 3
    for s in range(2*J): k += (u & 1) - (v & 1); u, v = T(u), T(v)
    k0 = k; ks = []
    for s in range(L - 4):
        k += (u & 1) - (v & 1); u, v = T(u), T(v); ks.append(k)
    # after the transient (8 steps) the debt grows by 1 every 2 steps
    assert ks[7] == 6 and all(ks[7 + 2*i] == 6 + i for i in range((L - 4 - 8)//2)), ks[:20]
    ok_a += 1
print(f"(a) source ones-run: debt 6 at step 8 then +1/2 per step, {ok_a} sources; (b) child ones-run: anchored at -1 with debt -5 after 12 steps, {ok_b} sources")
print("ALL CHECKS PASSED")
