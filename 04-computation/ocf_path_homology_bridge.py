#!/usr/bin/env python3
"""
ocf_path_homology_bridge.py — Connection between OCF and path homology

OCF: H(T) = I(Omega(T), 2) = sum_{k>=0} i_k * 2^k
Path homology: beta_p(T) = dim H_p(path complex of T)

Both involve the cycle/independence structure of Omega(T).
Is there a formula relating beta_p to I(Omega(T), x)?

Possible connections:
1. I(Omega,x) at x=-1 gives Euler characteristic of Omega_3 independence complex
2. beta_p might be related to log-concavity properties of I(Omega,x)
3. The independence polynomial coefficients encode topological information

Let's test: does beta_1 or beta_3 determine (or constrain) I(Omega,x)?

At n=6 exhaustive:
  C-phase (beta_1=1): what are the I(Omega,x) values?
  S-phase (beta_3=1): what are the I(Omega,x) values?
  Are they the same or different?

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, time
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def count_3cycles(A, n):
    """Return list of vertex sets that form directed 3-cycles."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append(frozenset([i,j,k]))
                if A[i][k] and A[k][j] and A[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def independence_poly(cycles_as_sets):
    """Compute I(Omega, x) where vertices are 3-cycle vertex sets.
    Two cycles conflict iff they share a vertex."""
    if not cycles_as_sets:
        return [1]

    # Unique vertex sets
    unique = list(set(cycles_as_sets))
    m = len(unique)

    # Build conflict graph
    conflicts = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if unique[i] & unique[j]:  # Share a vertex
                conflicts[i][j] = True
                conflicts[j][i] = True

    # Compute independence polynomial by enumeration
    # I(x) = sum_{S independent} x^|S|
    coeffs = [0] * (m + 1)
    coeffs[0] = 1

    for mask in range(1, 1 << m):
        # Check if this subset is independent
        independent = True
        bits = []
        for i in range(m):
            if mask & (1 << i):
                bits.append(i)
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if conflicts[bits[i]][bits[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(bits)] += 1

    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()

    return coeffs

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

# ===== n=6 EXHAUSTIVE: correlate path homology with independence polynomial =====
print("=" * 70)
print("n=6: Path homology vs Independence polynomial")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
t0 = time.time()

# Collect data by phase
phase_data = {'P': [], 'C': [], 'S': []}

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0
    phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')

    # Get 3-cycle structure and independence polynomial
    cycles = count_3cycles(A, n)
    unique_cycle_sets = list(set(cycles))
    ip = independence_poly(unique_cycle_sets)
    H = sum(c * 2**k for k, c in enumerate(ip))  # I(Omega, 2)

    phase_data[phase].append({
        'ip': tuple(ip),
        'H': H,
        'c3': len(unique_cycle_sets),
        'beta': [int(beta[k]) for k in range(min(5, len(beta)))]
    })

    if (bits + 1) % 10000 == 0:
        print(f"  {bits+1}/{total} ({time.time()-t0:.1f}s)")

print(f"\nDone in {time.time()-t0:.1f}s")

# Analyze by phase
for phase in ['P', 'C', 'S']:
    data = phase_data[phase]
    ip_counts = Counter(d['ip'] for d in data)
    H_counts = Counter(d['H'] for d in data)
    c3_counts = Counter(d['c3'] for d in data)

    print(f"\n{phase} ({len(data)} tournaments):")
    print(f"  IP distribution (top 10):")
    for ip, cnt in ip_counts.most_common(10):
        H_val = sum(c * 2**k for k, c in enumerate(ip))
        print(f"    I={list(ip)}, H={H_val}: {cnt}")
    print(f"  H values: {dict(sorted(H_counts.items()))}")
    print(f"  c3 values: {dict(sorted(c3_counts.items()))}")

# Cross-tabulate: which IPs appear in which phases?
print("\n" + "=" * 70)
print("IP -> PHASE MAPPING")
print("=" * 70)

all_ips = {}
for phase in ['P', 'C', 'S']:
    for d in phase_data[phase]:
        ip = d['ip']
        if ip not in all_ips:
            all_ips[ip] = Counter()
        all_ips[ip][phase] += 1

for ip in sorted(all_ips.keys(), key=lambda x: sum(c * 2**k for k, c in enumerate(x)), reverse=True)[:20]:
    phases = all_ips[ip]
    H_val = sum(c * 2**k for k, c in enumerate(ip))
    total_count = sum(phases.values())
    print(f"  I={list(ip)}, H={H_val}: {dict(phases)} ({total_count} total)")

# Key question: does IP determine phase?
print("\n" + "=" * 70)
print("DOES IP DETERMINE PHASE?")
print("=" * 70)

ip_determines = 0
ip_mixed = 0
for ip, phases in all_ips.items():
    if len(phases) == 1:
        ip_determines += 1
    else:
        ip_mixed += 1

print(f"IP with unique phase: {ip_determines}")
print(f"IP with mixed phases: {ip_mixed}")

# Show the mixed ones
print("\nMixed-phase IPs:")
for ip, phases in all_ips.items():
    if len(phases) > 1:
        H_val = sum(c * 2**k for k, c in enumerate(ip))
        print(f"  I={list(ip)}, H={H_val}: {dict(phases)}")
