#!/usr/bin/env python3
"""Benchmark c_4 Burnside vs closed-form at various n values."""

from a051240_closedform_v3 import compute_ck_burnside, compute_c4_closedform, lcm_list
from time import perf_counter

def gen_parts(rem, maxp, parts):
    if rem == 0:
        yield list(parts)
        return
    for p in range(min(rem, maxp), 0, -1):
        for m in range(1, rem // p + 1):
            parts.append((p, m))
            yield from gen_parts(rem - m * p, p - 1, parts)
            parts.pop()

# Find the partition with the largest L for each n
for n in [15, 20, 25, 30]:
    parts = list(gen_parts(n, n, []))
    max_L_part = max(parts, key=lambda p: lcm_list([r for r, m in p]))
    max_L = lcm_list([r for r, m in max_L_part])
    print(f"\nn={n}: {len(parts)} partitions, max L = {max_L}")
    print(f"  Worst partition: {max_L_part}")

    # Benchmark on just the high-L partitions
    high_L = [p for p in parts if lcm_list([r for r, m in p]) > max_L // 10]
    print(f"  High-L partitions (L > {max_L//10}): {len(high_L)}")

    if n <= 25:
        t0 = perf_counter()
        for part in parts:
            compute_ck_burnside(part, 4)
        dt_burn = perf_counter() - t0
    else:
        # Only benchmark on high-L for large n
        t0 = perf_counter()
        for part in high_L[:20]:
            compute_ck_burnside(part, 4)
        dt_burn = perf_counter() - t0
        print(f"  (Burnside on {min(20,len(high_L))} high-L only)")

    t0 = perf_counter()
    if n <= 25:
        for part in parts:
            compute_c4_closedform(part)
    else:
        for part in high_L[:20]:
            compute_c4_closedform(part)
    dt_cf = perf_counter() - t0

    print(f"  Burnside:    {dt_burn:.3f}s")
    print(f"  Closed-form: {dt_cf:.3f}s")
    if dt_cf > 0:
        print(f"  Speedup: {dt_burn/dt_cf:.1f}x")

    # Single worst-case partition
    t0 = perf_counter()
    ref = compute_ck_burnside(max_L_part, 4)
    dt1 = perf_counter() - t0

    t0 = perf_counter()
    cf = compute_c4_closedform(max_L_part)
    dt2 = perf_counter() - t0

    print(f"  Worst-case single: Burnside {dt1:.4f}s, closed-form {dt2:.4f}s, speedup {dt1/max(dt2,1e-9):.1f}x")
    assert ref == cf, f"MISMATCH: {ref} vs {cf}"
