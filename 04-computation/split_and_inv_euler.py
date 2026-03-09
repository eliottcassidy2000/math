#!/usr/bin/env python3
"""
Split b_matrix_quad.txt into individual b-files and compute inverse Euler transforms.

b_matrix_quad.txt format: n A002724 A006383 A122082 A007139
b000568_v4.txt format: n a(n)

Inverse Euler transform:
  c(n) = n*a(n) - sum_{k=1}^{n-1} c(k)*a(n-k)
  b(n) = (c(n) - sum_{d|n, d<n} d*b(d)) / n
"""

import os

DIR = os.path.dirname(os.path.abspath(__file__))

# --- Part (a): Split b_matrix_quad.txt into 4 individual b-files ---

quad_seqs = ["002724", "006383", "122082", "007139"]
quad_data = {s: [] for s in quad_seqs}

with open(os.path.join(DIR, "b_matrix_quad.txt")) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        n = int(parts[0])
        for i, seq in enumerate(quad_seqs):
            quad_data[seq].append((n, int(parts[i + 1])))

for seq in quad_seqs:
    outfile = os.path.join(DIR, f"b{seq}_quad.txt")
    with open(outfile, "w") as f:
        for n, val in quad_data[seq]:
            f.write(f"{n} {val}\n")
    print(f"Wrote {outfile} ({len(quad_data[seq])} terms)")

# --- Helper: inverse Euler transform ---

def inverse_euler_transform(a_dict):
    """
    Given a_dict mapping n -> a(n) for n >= 1 (with a(0) ignored if present),
    compute b(n) via inverse Euler transform.

    c(n) = n*a(n) - sum_{k=1}^{n-1} c(k)*a(n-k)
    b(n) = (c(n) - sum_{d|n, d<n} d*b(d)) / n
    """
    nmax = max(a_dict.keys())
    a = {n: a_dict[n] for n in a_dict}

    c = {}
    b = {}

    for n in range(1, nmax + 1):
        if n not in a:
            break
        # c(n) = n*a(n) - sum_{k=1}^{n-1} c(k)*a(n-k)
        cn = n * a[n]
        for k in range(1, n):
            if k in c and (n - k) in a:
                cn -= c[k] * a[n - k]
        c[n] = cn

        # b(n) = (c(n) - sum_{d|n, d<n} d*b(d)) / n
        bn = c[n]
        for d in range(1, n):
            if n % d == 0 and d in b:
                bn -= d * b[d]
        assert bn % n == 0, f"b({n}) not integer: remainder {bn % n}"
        b[n] = bn // n

    return b

# --- Part (b): Inverse Euler of A000568 -> A000171 (connected tournaments) ---

a568 = {}
with open(os.path.join(DIR, "b000568_v4.txt")) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) == 2:
            n, val = int(parts[0]), int(parts[1])
            if n >= 1:
                a568[n] = val

b171 = inverse_euler_transform(a568)
outfile = os.path.join(DIR, "b000171_inv.txt")
with open(outfile, "w") as f:
    for n in sorted(b171.keys()):
        f.write(f"{n} {b171[n]}\n")
print(f"Wrote {outfile} ({len(b171)} terms, A000171 = connected tournaments)")

# --- Part (c): Inverse Euler of A002724 -> A001242 (connected binary matrices) ---

a2724 = {}
for n, val in quad_data["002724"]:
    if n >= 1:
        a2724[n] = val

b1242 = inverse_euler_transform(a2724)
outfile = os.path.join(DIR, "b001242_inv.txt")
with open(outfile, "w") as f:
    for n in sorted(b1242.keys()):
        f.write(f"{n} {b1242[n]}\n")
print(f"Wrote {outfile} ({len(b1242)} terms, A001242 = connected binary matrices)")

# Print first few terms for verification
print("\nA000171 (connected tournaments), first 15 terms:")
for n in range(1, min(16, max(b171.keys()) + 1)):
    if n in b171:
        print(f"  b({n}) = {b171[n]}")

print("\nA001242 (connected binary matrices), first 15 terms:")
for n in range(1, min(16, max(b1242.keys()) + 1)):
    if n in b1242:
        print(f"  b({n}) = {b1242[n]}")
