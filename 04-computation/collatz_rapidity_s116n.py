#!/usr/bin/env python3
"""collatz_rapidity_s116n.py — The Collatz conjecture in rapidity/Cayley space.

The Collatz map T(n) = n/2 (even) or 3n+1 (odd) is a dynamical system
on the positive integers. We study it through the Cayley transform lens.

Key insight: rapidity rho(n) = arctanh((n-1)/(n+1)) = (1/2)*ln(n).
In rapidity space, Collatz is:
  - Even step: rho -> rho - ln(2)/2  (subtract half a log)
  - Odd step: rho -> (1/2)*ln(3n+1) ~ rho + ln(3)/2  (add half a log)

The average drift per odd-even cycle: ln(3/4)/2 = -0.1438... < 0.
This negative drift is WHY Collatz converges (heuristically).

The higher-dimensional space: the 2-adic integers Z_2.
Each n's orbit encodes a path through the binary tree of mod-2^k residues.

Session: kind-pasteur-2026-03-16-S116n32 (continued)
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from math import log, sqrt, atanh, gcd, isqrt
from fractions import Fraction

print()
print("  THE COLLATZ CONJECTURE IN RAPIDITY SPACE")
print()
print("=" * 70)
print()

# ============================================================
print("  I. RAPIDITY IS HALF THE LOG")
print("  " + "-" * 50)
print()
print("  Cayley address: x_n = (n-1)/(n+1) = Q^{-1}(n)")
print("  Rapidity: rho(n) = arctanh(x_n) = (1/2)*ln(n)")
print()
print("  Proof: arctanh((n-1)/(n+1)) = (1/2)*ln((1+(n-1)/(n+1))/(1-(n-1)/(n+1)))")
print("       = (1/2)*ln((2n/(n+1)) / (2/(n+1)))")
print("       = (1/2)*ln(n)")
print()

# Verify
for n in [2, 3, 7, 42, 100]:
    x_n = (n-1)/(n+1)
    rho = atanh(x_n)
    half_ln = log(n)/2
    print(f"  n={n:4d}: x_n = {x_n:.6f}, rho = {rho:.6f}, ln(n)/2 = {half_ln:.6f}, "
          f"match: {abs(rho - half_ln) < 1e-10}")
print()

# ============================================================
print("  II. COLLATZ IN RAPIDITY SPACE")
print("  " + "-" * 50)
print()

def collatz_step(n):
    if n % 2 == 0:
        return n // 2
    else:
        return 3*n + 1

def collatz_orbit(n, max_steps=1000):
    orbit = [n]
    while n != 1 and len(orbit) < max_steps:
        n = collatz_step(n)
        orbit.append(n)
    return orbit

print("  The rapidity changes:")
print(f"  Even step n -> n/2: rho -> rho - ln(2)/2 = rho - {log(2)/2:.6f}")
print(f"  Odd step n -> 3n+1: rho -> (1/2)*ln(3n+1)")
print(f"    For large n: ~ rho + ln(3)/2 = rho + {log(3)/2:.6f}")
print(f"    Correction: (1/2)*ln(3n+1) = (1/2)*ln(3) + (1/2)*ln(n+1/3)")
print(f"              = (1/2)*ln(3) + (1/2)*ln(n) + (1/2)*ln(1+1/(3n))")
print(f"              ~ rho + {log(3)/2:.6f} + O(1/n)")
print()
print(f"  After an odd step, the result 3n+1 is always EVEN.")
print(f"  So an odd step is ALWAYS followed by at least one even step.")
print(f"  Combined 'odd + k halvings': rho -> rho + ln(3)/2 - k*ln(2)/2")
print(f"  = rho + (ln(3) - k*ln(2))/2")
print()

# The key: what is the average k?
# After 3n+1, v_2(3n+1) determines k.
# For random odd n, the probability that 3n+1 is divisible by 2^k is 1/2^k
# (since 3n+1 mod 2^k is roughly uniform for random odd n).
# Actually: 3n+1 for odd n. n odd means n = 2m+1, so 3n+1 = 6m+4 = 2(3m+2).
# So 3n+1 is always divisible by 2 but 3m+2 is even iff m is even iff n = 4j+1.
# So v_2(3n+1) = 1 if n = 3 mod 4, v_2(3n+1) >= 2 if n = 1 mod 4.
# For n = 1 mod 4: 3n+1 = 3(4j+1)+1 = 12j+4 = 4(3j+1). v_2 >= 2.
# 3j+1 even iff j odd iff n = 5 mod 8. Then v_2 >= 3.
# This is the 2-adic tree structure!

print("  THE 2-ADIC TREE OF COLLATZ:")
print()
print("  v_2(3n+1) depends on n mod 2^k:")
for k in range(1, 8):
    # For odd n in range, compute v_2(3n+1)
    v2_counts = {}
    for n in range(1, 2**(k+2), 2):  # odd n up to 2^{k+2}
        val = 3*n + 1
        v = 0
        while val % 2 == 0:
            v += 1
            val //= 2
        res = n % (2**k)
        if res not in v2_counts:
            v2_counts[res] = set()
        v2_counts[res].add(v)

    # Show residue classes with constant v_2
    for res in sorted(v2_counts.keys()):
        if res % 2 == 0: continue  # only odd residues
        vs = v2_counts[res]
        if len(vs) == 1:
            v = list(vs)[0]
            if v >= k-1:
                print(f"    n = {res} mod {2**k}: v_2(3n+1) = {v}")
print()

# ============================================================
print("  III. THE RAPIDITY ORBIT")
print("  " + "-" * 50)
print()

# Show the rapidity orbit for a few starting values
for n0 in [27, 42, 97, 871]:
    orbit = collatz_orbit(n0)
    print(f"  n = {n0}: orbit length = {len(orbit)}, max = {max(orbit)}")

    # Rapidity values
    rhos = [log(n)/2 if n > 0 else 0 for n in orbit]

    # Show first 20 and last 5
    rapidity_str = ', '.join(f'{r:.3f}' for r in rhos[:15])
    print(f"    Rapidity: {rapidity_str}, ...")
    print(f"    Start rho = {rhos[0]:.4f}, peak rho = {max(rhos):.4f}, "
          f"end rho = {rhos[-1]:.4f} (=0 at n=1)")

    # Count odd steps and even steps
    odd_steps = sum(1 for i in range(len(orbit)-1) if orbit[i] % 2 == 1)
    even_steps = len(orbit) - 1 - odd_steps
    print(f"    Odd steps: {odd_steps}, Even steps: {even_steps}, "
          f"ratio E/O = {even_steps/odd_steps:.4f}")
    print(f"    Total rapidity change: {rhos[-1] - rhos[0]:.4f} = "
          f"{odd_steps}*ln(3)/2 - {even_steps}*ln(2)/2 + corrections")
    predicted = odd_steps * log(3)/2 - even_steps * log(2)/2
    print(f"    Predicted (no corrections): {predicted:.4f}, "
          f"actual: {rhos[-1] - rhos[0]:.4f}")
    print()

# ============================================================
print("  IV. THE EVEN/ODD RATIO AND ln(3)/ln(2)")
print("  " + "-" * 50)
print()
print(f"  For rapidity to decrease to 0 (reaching n=1):")
print(f"  Need: odd_steps * ln(3)/2 - even_steps * ln(2)/2 < -rho_0")
print(f"  => even/odd > ln(3)/ln(2) = {log(3)/log(2):.6f}")
print(f"  This is the CRITICAL RATIO.")
print()
print(f"  ln(3)/ln(2) = {log(3)/log(2):.10f}")
print(f"  Continued fraction: [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, ...]")
print(f"  Convergents: 1, 2, 3/2, 8/5, 19/12, 65/41, 84/53, ...")
print()

# Verify the ratio converges to ln(3)/ln(2)
print("  Even/Odd ratio for first N steps (starting from n=27):")
orbit = collatz_orbit(27)
running_odd = 0
running_even = 0
for i in range(len(orbit)-1):
    if orbit[i] % 2 == 1:
        running_odd += 1
    else:
        running_even += 1
    if (i+1) % 20 == 0 and running_odd > 0:
        r = running_even / running_odd
        print(f"    After {i+1:4d} steps: E/O = {r:.6f} "
              f"(target: {log(3)/log(2):.6f}, diff: {r - log(3)/log(2):+.6f})")
print()

# Large-scale statistics
print("  Large-scale E/O ratios:")
total_odd = 0
total_even = 0
for n0 in range(3, 10001, 2):  # odd starting values
    orbit = collatz_orbit(n0)
    for i in range(len(orbit)-1):
        if orbit[i] % 2 == 1:
            total_odd += 1
        else:
            total_even += 1

ratio = total_even / total_odd
print(f"  Over all odd n in [3, 10000]: E/O = {ratio:.8f}")
print(f"  ln(3)/ln(2) = {log(3)/log(2):.8f}")
print(f"  Difference: {ratio - log(3)/log(2):+.8f}")
print()

# ============================================================
print("  V. THE 2-ADIC ENCODING: COLLATZ AS A PATH IN Z_2")
print("  " + "-" * 50)
print()

# For each odd n, the "parity sequence" of the orbit encodes a 2-adic integer.
# Start with odd n. Apply T repeatedly. Record O for odd, E for even.
# The sequence of O/E bits (mod 2 of each orbit value) IS the 2-adic expansion.

def parity_sequence(n, length=32):
    """Record the parities of the Collatz orbit."""
    bits = []
    for _ in range(length):
        bits.append(n % 2)
        if n == 1: break
        n = collatz_step(n)
    return bits

print("  Parity sequences (0=even, 1=odd) for small odd n:")
for n in [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 27]:
    bits = parity_sequence(n, 24)
    bit_str = ''.join(str(b) for b in bits)
    # Count transitions
    transitions = sum(1 for i in range(len(bits)-1) if bits[i] != bits[i+1])
    print(f"  n={n:3d}: {bit_str}  (transitions: {transitions})")
print()

# The "compressed" Collatz: skip even steps, only record odd values
def odd_orbit(n, max_steps=100):
    """Return sequence of odd values in Collatz orbit."""
    odds = [n] if n % 2 == 1 else []
    for _ in range(max_steps):
        if n == 1 and len(odds) > 1: break
        n = collatz_step(n)
        if n % 2 == 1:
            odds.append(n)
    return odds

print("  Compressed (odd-only) orbits:")
for n in [3, 7, 27, 41, 97]:
    odds = odd_orbit(n)
    print(f"  n={n:3d}: {' -> '.join(str(x) for x in odds[:12])}{'...' if len(odds)>12 else ''}")
print()

# The "shortcut" map S(n) = (3n+1)/2^{v_2(3n+1)} for odd n
def shortcut(n):
    """S(n) = (3n+1) / 2^v_2(3n+1) for odd n."""
    val = 3*n + 1
    while val % 2 == 0:
        val //= 2
    return val

print("  Shortcut map S(n) = (3n+1)/2^v (odd -> odd):")
print("  n mod 4 = 1: 3n+1 = 0 mod 4, so v >= 2.")
print("  n mod 4 = 3: 3n+1 = 2 mod 4, so v = 1 exactly.")
print()
for n in range(1, 40, 2):
    s = shortcut(n)
    v = 0
    temp = 3*n+1
    while temp % 2 == 0:
        v += 1
        temp //= 2
    direction = "UP" if s > n else "DOWN"
    print(f"  S({n:3d}) = {s:5d}  (v_2 = {v}, {direction:4s}, "
          f"n mod 4 = {n%4}, n mod 8 = {n%8})")
print()

# ============================================================
print("  VI. THE HIGHER-DIMENSIONAL ENCODING")
print("  " + "-" * 50)
print()
print("  The Collatz orbit of n encodes a point in the infinite-dimensional space")
print("  Z_2 x Z_3 x Z_5 x ... (the profinite completion of Z).")
print()
print("  But the DYNAMICS depend only on the 2-adic and 3-adic components:")
print("  - The 2-adic part determines the halving depth (v_2 structure)")
print("  - The 3-adic part tracks how 3n+1 interacts with powers of 3")
print()

# Track v_2 and v_3 through orbits
def v_p(n, p):
    if n == 0: return float('inf')
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v

print("  v_2 and v_3 traces for n=27 orbit:")
orbit = collatz_orbit(27)
print(f"  {'step':>4s}  {'n':>6s}  v_2  v_3  {'rho':>8s}  {'n mod 12':>8s}")
for i, n in enumerate(orbit[:30]):
    v2 = v_p(n, 2)
    v3 = v_p(n, 3)
    rho = log(n)/2 if n > 0 else 0
    print(f"  {i:4d}  {n:6d}  {v2:3d}  {v3:3d}  {rho:8.4f}  {n%12:8d}")
print(f"  ... ({len(orbit)} total steps)")
print()

# ============================================================
print("  VII. THE THREE PRIMES {2, 3, 5} AND COLLATZ")
print("  " + "-" * 50)
print()
print("  The Collatz map involves exactly two primes: 2 (halving) and 3 (tripling).")
print("  These are the first two Hurwitz primes!")
print()
print("  In our framework:")
print("  - 2 = the doubler (formal group evaluation point)")
print("  - 3 = the first Cayley pair of 7 (the forbidden prime)")
print("  - 3n+1 = 3(n + 1/3) = 3(n + Q^{-1}(2))")
print("  The odd step ADDS the Cayley address of 2, then MULTIPLIES by 3.")
print()
print("  The VS chain diverges when 5 enters (via 4+1=5).")
print("  The Collatz map never involves 5 directly.")
print("  But the ORBIT passes through multiples of 5:")

# How often do Collatz orbits hit multiples of 5?
hit_5 = 0
total_steps = 0
for n0 in range(1, 10001):
    orbit = collatz_orbit(n0)
    for n in orbit:
        total_steps += 1
        if n % 5 == 0:
            hit_5 += 1

print(f"  Fraction of orbit values divisible by 5: {hit_5}/{total_steps} = {hit_5/total_steps:.6f}")
print(f"  Expected (random): 1/5 = 0.200000")
print(f"  Ratio: {(hit_5/total_steps)*5:.6f}")
print()

# ============================================================
print("  VIII. COLLATZ AND THE FORMAL GROUP")
print("  " + "-" * 50)
print()
print("  In the Cayley formal group F(x,y) = (x+y)/(1+xy):")
print("  [n](x) = Q^{-1}(Q(x)^n)")
print()
print("  The Collatz even step n -> n/2:")
print("  In formal group: 'divide by 2' = apply the INVERSE of [2].")
print("  [2](x) = 2x/(1+x^2). Inverse: solve 2y/(1+y^2) = x.")
print("  y^2*x - 2y + x = 0. y = (2 +/- sqrt(4-4x^2))/(2x) = (1 +/- sqrt(1-x^2))/x")
print("  The formal inverse exists only as a MULTIVALUED function.")
print("  This is the formal group saying: HALVING IS NOT A HOMOMORPHISM.")
print()
print("  The Collatz odd step n -> 3n+1:")
print("  In Cayley addresses: x_{3n+1} = (3n)/(3n+2) = 3n/(3n+2)")
print("  While x_n = (n-1)/(n+1).")
print("  Relationship: x_{3n+1} = (3*((n+1)*x_n/(1+x_n^{-1}))... complex.")
print("  The odd step is NOT a formal group operation.")
print()
print("  KEY INSIGHT: Collatz is not a formal group map.")
print("  It mixes MULTIPLICATION (3n) with ADDITION (+1),")
print("  but the formal group only sees one operation at a time.")
print("  The +1 is the PERTURBATION that makes Collatz non-algebraic.")
print()

# ============================================================
print("  IX. THE DRIFT IN RAPIDITY SPACE")
print("  " + "-" * 50)
print()

# For each "odd step + subsequent halvings" cycle:
# rho_new = (1/2)*ln((3n+1)/2^k) = (1/2)*ln(3) + (1/2)*ln(n+1/3) - k*ln(2)/2
# ~ rho + (ln(3) - k*ln(2))/2
# where k = v_2(3n+1).

# Distribution of k = v_2(3n+1) for odd n:
v2_dist = {}
for n in range(1, 100000, 2):
    val = 3*n + 1
    v = 0
    while val % 2 == 0:
        v += 1
        val //= 2
    v2_dist[v] = v2_dist.get(v, 0) + 1

total = sum(v2_dist.values())
print("  Distribution of v_2(3n+1) for odd n in [1, 100000):")
expected_drift = 0
for v in sorted(v2_dist.keys()):
    count = v2_dist[v]
    freq = count / total
    drift = (log(3) - v*log(2)) / 2
    expected_drift += freq * drift
    print(f"    v_2 = {v}: {count:6d}/{total} = {freq:.6f}  "
          f"(expected: {1/2**v:.6f})  "
          f"drift = {drift:+.6f}")

print()
print(f"  Expected drift per odd step: {expected_drift:+.8f}")
print("  Theoretical: sum_{k=1}^inf (1/2^k) * (ln(3) - k*ln(2))/2")
print(f"  = (ln(3)/2) * sum(1/2^k) - (ln(2)/2) * sum(k/2^k)")
print(f"  = (ln(3)/2) * 1 - (ln(2)/2) * 2 = ln(3)/2 - ln(2)")
print(f"  = {log(3)/2 - log(2):+.8f}")
print()
print(f"  Average drift per odd step = ln(3)/2 - ln(2) = (ln(3) - 2*ln(2))/2")
print(f"  = ln(3/4)/2 = {log(3/4)/2:+.8f}")
print(f"  This is NEGATIVE: orbits drift DOWNWARD in rapidity on average.")
print()

# ============================================================
print("  X. THE STOPPING TIME AND RAPIDITY RETURN")
print("  " + "-" * 50)
print()

# The stopping time sigma(n) = first time orbit drops below n.
# In rapidity: first time rho(orbit) < rho(n).
# Average stopping time should be proportional to rho(n) / |drift|.

print("  Stopping times vs rapidity for odd n up to 10000:")
print(f"  {'n':>6s}  {'rho(n)':>8s}  {'stop_time':>10s}  {'total_time':>11s}  {'rho/|drift|':>12s}")
for n0 in [3, 7, 27, 97, 255, 871, 2463, 4591, 6171, 8191, 9663]:
    orbit = collatz_orbit(n0)
    rho = log(n0)/2
    # Stopping time: first i where orbit[i] < n0
    stop = len(orbit)
    for i in range(1, len(orbit)):
        if orbit[i] < n0:
            stop = i
            break
    drift = abs(log(3/4)/2)
    print(f"  {n0:6d}  {rho:8.4f}  {stop:10d}  {len(orbit)-1:11d}  {rho/drift:12.2f}")
print()

# ============================================================
print("  XI. THE 2-ADIC TREE: COLLATZ ENCODES Z_2")
print("  " + "-" * 50)
print()
print("  Each odd n defines a path in the binary tree by its mod-2^k behavior.")
print("  The Collatz dynamics depend on n mod 2^k for increasing k.")
print()

# Show how increasing k-refinement determines the trajectory
print("  n = 27 = 11011 in binary. Residues:")
n = 27
for k in range(1, 8):
    res = n % (2**k)
    # What does this tell us about the first step?
    if k == 1:
        print(f"    mod 2^{k} = {res}: odd -> apply 3n+1")
    elif k == 2:
        print(f"    mod 2^{k} = {res}: 3*{res}+1 = {3*res+1}, "
              f"v_2 = {v_p(3*res+1, 2)}")
    else:
        step_val = 3*res + 1
        v = v_p(step_val, 2)
        next_val = step_val // (2**v)
        print(f"    mod 2^{k} = {res}: 3*{res}+1 = {step_val}, "
              f"v_2 = {v}, next odd = {next_val} mod 2^{k}")
print()

# The "tree address" of n: sequence of choices at each level
print("  The TREE ADDRESS of each n (v_2 sequence of shortcut map):")
for n0 in [1, 3, 5, 7, 9, 11, 13, 15, 27, 41]:
    n = n0
    addr = []
    for _ in range(8):
        if n == 1:
            break
        if n % 2 == 0:
            # shouldn't happen in shortcut
            break
        val = 3*n + 1
        v = v_p(val, 2)
        addr.append(v)
        n = val // (2**v)
    print(f"  n={n0:3d}: v_2 sequence = {addr}")
print()

print("  These v_2 sequences ARE coordinates in the 'higher-dimensional space'")
print("  that the Collatz conjecture encodes. Each sequence is a path through")
print("  the tree of possible halvings, indexed by the 2-adic expansion of n.")
print()

# ============================================================
print("  XII. THE CONNECTION: VS CHAIN AND COLLATZ")
print("  " + "-" * 50)
print()
print("  VS chain: n -> prod{p : (p-1)|n}")
print("  Collatz:  n -> n/2 or 3n+1")
print()
print("  Both are iterated arithmetic maps with attractors:")
print("  - VS has attractor 1806 (density 0.64)")
print("  - Collatz has attractor {1,2,4} (conjectured density 1.0)")
print()
print("  The VS basin: ALL odd n converge (vs_prod(odd) = 2).")
print("  The Collatz basin: ALL n converge (conjectured).")
print("  Both are controlled by 2-adic structure:")
print("  - VS diverges when 4|n (prime 5 enters)")
print("  - Collatz behavior depends on v_2(3n+1)")
print()
print("  The formal group connection:")
print("  - VS uses primes p with (p-1)|n: divisibility of TOTIENT")
print("  - Collatz uses v_2(3n+1): powers of 2 in an AFFINE image")
print("  - Both probe the multiplicative structure of n through the prime 2")
print()
print("  42 appears in both:")
print("  - VS: 42 is the bridge to the fixed point 1806")
print("  - Collatz: 42 -> 21 -> 64 -> 32 -> 16 -> 8 -> 4 -> 2 -> 1")

orbit_42 = collatz_orbit(42)
print(f"    Orbit: {' -> '.join(str(x) for x in orbit_42)}")
print(f"    Steps: {len(orbit_42)-1}. Via 21 -> 64 = 2^6 (fast descent).")
print()

# The rapidity journey of 42 in Collatz
print("  42's rapidity journey:")
for i, n in enumerate(orbit_42):
    if n > 0:
        rho = log(n)/2
        print(f"    step {i}: n={n:4d}, rho={rho:+.4f}", end="")
        if i > 0:
            delta = log(n)/2 - log(orbit_42[i-1])/2
            print(f"  (delta = {delta:+.4f})", end="")
        print()
print()

# ============================================================
print("  XIII. SYNTHESIS: THE THREE LAYERS")
print("  " + "-" * 50)
print()
print("  Layer 1: RAPIDITY (the 'height' coordinate)")
print("    rho(n) = (1/2)*ln(n)")
print("    Collatz has drift ln(3/4)/2 = -0.144 per odd step")
print("    Reaching n=1 means reaching rho=0")
print()
print("  Layer 2: THE 2-ADIC TREE (the 'horizontal' coordinates)")
print("    n mod 2^k for k=1,2,3,... determines the halving depth")
print("    v_2 sequence = coordinates in the infinite-dimensional cube {1,2,3,...}^N")
print("    The 'higher-dimensional space' is Z_2 (2-adic integers)")
print()
print("  Layer 3: THE FORMAL GROUP (the algebraic structure)")
print("    F(x,y) = (x+y)/(1+xy), log_F = arctanh = rapidity")
print("    Height infinity at p=2: the 2-adic structure is maximally degenerate")
print("    Collatz is NOT a formal group map (the +1 breaks homomorphism)")
print("    But rapidity = log_F is the right coordinate for tracking orbits")
print()
print("  The Collatz conjecture lives in the tension between:")
print("  - The formal group (algebraic, multiplicative, structured)")
print("  - The +1 perturbation (arithmetic, additive, chaotic)")
print("  The conjecture says: the algebraic drift always wins.")
print()
