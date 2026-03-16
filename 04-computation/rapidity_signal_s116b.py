#!/usr/bin/env python3
"""rapidity_signal_s116b.py — The H-spectrum as a rapidity signal

The achievable H values form a DISCRETE SIGNAL in rapidity space.
What does this signal look like? What are its gaps? Its density?
"""
from math import sqrt, log, pi, comb, factorial
from itertools import permutations
from collections import Counter

def count_ham_paths(adj, n):
    """Count Hamiltonian paths in tournament given by adjacency matrix."""
    count = 0
    for perm in permutations(range(n)):
        valid = True
        for i in range(n-1):
            if not adj[perm[i]][perm[i+1]]:
                valid = False
                break
        if valid:
            count += 1
    return count

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(arcs)):
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(arcs):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

print("THE H-SPECTRUM AS A RAPIDITY SIGNAL")
print("="*60)
print()

# Compute all achievable H values for small n
h_sets = {}
for n in range(3, 8):
    h_counter = Counter()
    for adj in all_tournaments(n):
        h = count_ham_paths(adj, n)
        h_counter[h] += 1
    h_sets[n] = h_counter
    achievable = sorted(h_counter.keys())
    print(f"  n={n}: {len(achievable)} distinct H values")
    print(f"    range: {achievable[0]} to {achievable[-1]}")
    print(f"    values: {achievable}")
    # Rapidity range
    raps = [log(h)/2 for h in achievable if h > 0]
    print(f"    rapidity range: [{min(raps):.4f}, {max(raps):.4f}]")
    print()

print()

# ============================================================
print("THE GAPS IN RAPIDITY SPACE")
print("="*60)
print()

for n in [5, 6, 7]:
    achievable = sorted(h_sets[n].keys())
    print(f"  n={n}: H-values and rapidity gaps:")
    print(f"  {'H':>6s}  {'rapidity':>10s}  {'gap':>10s}  {'gap/octave':>12s}")
    print("  " + "-"*45)
    octave = log(2)/2
    prev_rap = None
    for h in achievable:
        rap = log(h)/2
        if prev_rap is not None:
            gap = rap - prev_rap
            print(f"  {h:6d}  {rap:10.6f}  {gap:10.6f}  {gap/octave:12.4f}")
        else:
            print(f"  {h:6d}  {rap:10.6f}")
        prev_rap = rap
    print()

# ============================================================
print("FORBIDDEN GAPS: WHERE H CANNOT BE")
print("="*60)
print()

for n in [5, 6, 7]:
    achievable = sorted(h_sets[n].keys())
    # Find all integers NOT achievable
    forbidden = [h for h in range(1, achievable[-1]+1) if h not in achievable]
    if forbidden:
        print(f"  n={n}: Forbidden H values in range [1, {achievable[-1]}]:")
        print(f"    {forbidden}")
        # Rapidity positions of forbidden values
        forb_raps = [(h, log(h)/2) for h in forbidden if h > 0]
        print(f"    Forbidden rapidities:")
        for h, r in forb_raps:
            print(f"      H={h}: rapidity = {r:.6f}")
        print()

# ============================================================
print("THE DENSITY OF H IN RAPIDITY SPACE")
print("="*60)
print()

for n in [5, 6, 7]:
    achievable = sorted(h_sets[n].keys())
    rap_min = log(achievable[0])/2
    rap_max = log(achievable[-1])/2
    rap_range = rap_max - rap_min
    density = len(achievable) / rap_range
    print(f"  n={n}: {len(achievable)} H-values over rapidity range {rap_range:.4f}")
    print(f"    Density = {density:.2f} values per unit rapidity")
    print(f"    Average gap = {rap_range/(len(achievable)-1):.6f}")
    octave = log(2)/2
    print(f"    Values per octave = {density * octave:.2f}")
    print()

# ============================================================
print("MULTIPLICATIVE STRUCTURE OF H")
print("="*60)
print()

print("  H values are always ODD (Redei's theorem).")
print("  In rapidity, all H-values avoid rapidity(2^k) for k >= 1.")
print("  The rapidity line has 'excluded bands' around even numbers.")
print()

for n in [5, 6, 7]:
    achievable = sorted(h_sets[n].keys())
    print(f"  n={n}: H values modulo small primes:")
    for p in [3, 5, 7]:
        residues = Counter(h % p for h in achievable)
        print(f"    mod {p}: {dict(sorted(residues.items()))}")
    print()

# ============================================================
print("H-VALUES AS RAPIDITY VECTORS (prime decomposition)")
print("="*60)
print()

for n in [6, 7]:
    achievable = sorted(h_sets[n].keys())
    print(f"  n={n}: H-values and their prime factorization:")
    basis = [3, 5, 7, 11, 13, 17, 19, 23]
    for h in achievable:
        temp = h
        factors = []
        for p in basis:
            while temp % p == 0:
                factors.append(p)
                temp //= p
        if temp > 1:
            factors.append(temp)
        factor_str = "*".join(str(f) for f in factors) if factors else "1"
        rap = log(h)/2
        print(f"    H={h:5d} = {factor_str:20s}  rapidity = {rap:.5f}")
    print()

# ============================================================
print("THE RAPIDITY FOURIER TRANSFORM OF THE H-SPECTRUM")
print("="*60)
print()
print("  The H-spectrum at each n is a discrete measure on the rapidity line.")
print("  Its Fourier transform reveals periodicities.")
print()

for n in [5, 7]:
    achievable = sorted(h_sets[n].keys())
    total = sum(h_sets[n].values())
    # Compute the "Fourier transform" at a few frequencies
    print(f"  n={n}: Fourier analysis of H-distribution in rapidity")
    print(f"    Total tournaments: {total}")
    print()
    print(f"    freq k   Re(F_k)      Im(F_k)      |F_k|")
    print("    " + "-"*50)

    for k_freq in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12]:
        re_sum = 0
        im_sum = 0
        for h, count in h_sets[n].items():
            rap = log(h)/2
            weight = count / total
            re_sum += weight * cos_approx(2*pi*k_freq*rap)
            im_sum += weight * sin_approx(2*pi*k_freq*rap)
        mag = sqrt(re_sum**2 + im_sum**2)
        print(f"    {k_freq:5d}   {re_sum:+.6f}   {im_sum:+.6f}   {mag:.6f}")
    print()

def cos_approx(x):
    """Cosine for our Fourier computation."""
    from math import cos
    return cos(x)

def sin_approx(x):
    """Sine for our Fourier computation."""
    from math import sin
    return sin(x)

# Redo with proper functions defined
from math import cos as cos_approx, sin as sin_approx

for n in [5, 7]:
    achievable = sorted(h_sets[n].keys())
    total = sum(h_sets[n].values())
    print(f"  n={n}: Fourier analysis of H-distribution in rapidity")
    print(f"    freq k   Re(F_k)      Im(F_k)      |F_k|")
    print("    " + "-"*50)

    for k_freq in [1, 2, 3, 4, 5, 6, 7, 8, 10, 12]:
        re_sum = 0
        im_sum = 0
        for h, count in h_sets[n].items():
            rap = log(h)/2
            weight = count / total
            re_sum += weight * cos_approx(2*pi*k_freq*rap)
            im_sum += weight * sin_approx(2*pi*k_freq*rap)
        mag = sqrt(re_sum**2 + im_sum**2)
        print(f"    {k_freq:5d}   {re_sum:+.6f}   {im_sum:+.6f}   {mag:.6f}")
    print()

# ============================================================
print("THE RAPIDITY GENERATING FUNCTION")
print("="*60)
print()
print("  Define Z_n(s) = sum over all tournaments T on n vertices of H(T)^{-s}.")
print("  = sum over all tournaments T of exp(-2s * rapidity(H(T))).")
print("  This is a PARTITION FUNCTION with rapidity as energy!")
print()

for n in [3, 4, 5, 6, 7]:
    total = sum(h_sets[n].values())
    print(f"  n={n} ({total} tournaments):")
    for s in [0.5, 1.0, 2.0]:
        Z = sum(count * h**(-s) for h, count in h_sets[n].items())
        Z_norm = Z / total
        print(f"    Z_n(s={s:.1f}) = {Z:.4f} (normalized: {Z_norm:.6f})")
    # Compute the "average rapidity" = -d/ds ln Z at s=0
    # = sum (count/total) * ln(H)
    avg_rap = sum(count/total * log(h)/2 for h, count in h_sets[n].items())
    var_rap = sum(count/total * (log(h)/2)**2 for h, count in h_sets[n].items()) - avg_rap**2
    print(f"    Average rapidity: {avg_rap:.6f}")
    print(f"    Rapidity variance: {var_rap:.6f}")
    print(f"    Rapidity std dev: {sqrt(var_rap):.6f}")
    # Expected H and its rapidity
    avg_h = sum(count/total * h for h, count in h_sets[n].items())
    print(f"    E[H] = {avg_h:.4f}, rapidity(E[H]) = {log(avg_h)/2:.6f}")
    print(f"    exp(2*avg_rapidity) = {exp(2*avg_rap):.4f} (geometric mean of H)")
    print()

print()
print("  The geometric mean of H = exp(2 * average rapidity).")
print("  This is DIFFERENT from E[H] because exp is convex (Jensen).")
print("  The gap: E[H] / geo_mean(H) = exp(variance_rapidity + ...) >= 1.")
print()

# ============================================================
print("RAPIDITY OF n! (Stirling)")
print("="*60)
print()
print("  n! grows superexponentially. rapidity(n!) = ln(n!)/2.")
print("  Stirling: ln(n!) ~ n*ln(n) - n + ln(2*pi*n)/2.")
print("  rapidity(n!) ~ n*ln(n)/2 - n/2 + ln(2*pi*n)/4.")
print("  = n*(rapidity(n) - 1/2) + ln(2*pi*n)/4.")
print()

for n in range(1, 13):
    fn = factorial(n)
    rap = log(fn)/2
    stirling = n*log(n)/2 - n/2 + log(2*pi*n)/4
    print(f"  {n:2d}!  = {fn:10d}   rapidity = {rap:8.4f}   Stirling = {stirling:8.4f}")

print()
print("  rapidity(n!) = n * rapidity(n) - n/2 + O(ln(n))")
print("  The dominant term is n * rapidity(n) = n * ln(n) / 2.")
print("  Factorial rapidity is ALMOST quadratic in rapidity(n).")
print()
print("  For H(T) ~ n!/2^{n-1} (Paley maximizer):")
print("  rapidity(H_max) ~ rapidity(n!) - (n-1)*rapidity(2)")
print("  = n*ln(n)/2 - n/2 - (n-1)*ln(2)/2")
print("  = n*(ln(n) - 1 - ln(2))/2 + ln(2)/2")
print("  = n*(ln(n/2) - 1)/2 + ln(2)/2")
print()
for n in [3, 7, 11]:
    fn = factorial(n)
    h_max_est = fn / 2**(n-1)
    rap_est = n*(log(n/2) - 1)/2 + log(2)/2
    print(f"  n={n:2d}: est rapidity(H_max) = {rap_est:.4f}")
