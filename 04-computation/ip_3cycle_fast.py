"""
opus-2026-05-22-S4: Fast 3-cycle IP for all-0 staircase T_k.

All 3-cycles in T_k come from pairs a<b: A(a,b)={2a,2a+1,2b}, B(a,b)={2a+1,2b,2b+1}.
Conflict rules:
  A(a,b) ~ A(c,d): share pair iff {a,b}∩{c,d}!=empty
  B(a,b) ~ B(c,d): share pair iff {a,b}∩{c,d}!=empty
  A(a,b) ~ B(c,d): conflict iff a=c or a=d or b=d  (NOT if b=c)
Uses bitmask adjacency for fast backtracking.
"""
import sys
import time
from math import comb

def build_conflict_graph(k):
    """Return list of cycles and adjacency bitmasks."""
    # Cycle indices: A(a,b) = 2*pair_idx, B(a,b) = 2*pair_idx+1
    # pair_idx for (a,b) with a<b: use sorted list
    pairs = [(a, b) for a in range(k) for b in range(a+1, k)]
    P = len(pairs)  # C(k,2)
    pair_idx = {(a,b): i for i, (a,b) in enumerate(pairs)}

    # Total cycles: 2P
    # cycle 2i = A(pairs[i]), cycle 2i+1 = B(pairs[i])
    N = 2 * P

    # Build adjacency bitmasks
    adj = [0] * N
    for i, (a, b) in enumerate(pairs):
        ai, bi = 2*i, 2*i+1  # A(a,b) and B(a,b)

        # A(a,b) always conflicts with B(a,b): share vertices {2a+1, 2b}
        adj[ai] |= (1 << bi)
        adj[bi] |= (1 << ai)

        for j, (c, d) in enumerate(pairs):
            if j == i: continue
            aj, bj = 2*j, 2*j+1

            # A(a,b) vs A(c,d): conflict iff {a,b}∩{c,d}!=empty
            if len({a,b} & {c,d}) > 0:
                adj[ai] |= (1 << aj)

            # B(a,b) vs B(c,d): same
            if len({a,b} & {c,d}) > 0:
                adj[bi] |= (1 << bj)

            # A(a,b) vs B(c,d): conflict iff a=c or a=d or b=d
            if a == c or a == d or b == d:
                adj[ai] |= (1 << bj)

            # B(a,b) vs A(c,d): conflict iff (A(c,d) vs B(a,b)) c=a or c=b or d=b
            if c == a or c == b or d == b:
                adj[bi] |= (1 << aj)

    return N, adj

def compute_ip_bitmask(k, max_degree=None):
    """Compute independence polynomial coefficients via backtracking with bitmask."""
    t0 = time.time()
    N, adj = build_conflict_graph(k)

    if max_degree is None:
        max_degree = N

    ip = [0] * (max_degree + 1)
    ip[0] = 1

    all_bits = (1 << N) - 1

    def backtrack(v, size, forbidden_mask):
        """v = next vertex to consider, forbidden_mask = cannot include."""
        for i in range(v, N):
            if (forbidden_mask >> i) & 1:
                continue
            if size + 1 <= max_degree:
                ip[size + 1] += 1
                backtrack(i + 1, size + 1, forbidden_mask | adj[i])

    backtrack(0, 0, 0)

    # Trim trailing zeros
    while len(ip) > 1 and ip[-1] == 0:
        ip.pop()

    elapsed = time.time() - t0
    return ip, elapsed

def verify_formula(k, ip, formulas):
    """Verify the alpha_m formula against computed IP."""
    errors = []
    for m in range(1, len(ip)):
        if m not in formulas:
            continue
        pred = sum(c * comb(k, 2*m - j) for j, c in enumerate(formulas[m]))
        if pred != ip[m]:
            errors.append(f"alpha_{m}: formula={pred}, computed={ip[m]}")
    return errors

# Known formulas: alpha_m = sum_j coeffs[j] * C(k, 2m-j)
FORMULAS = {
    1: [2],                                         # 2*C(k,2)
    2: [12, 1],                                     # 12*C(k,4) + 1*C(k,3)
    3: [120, 20],                                   # 120*C(k,6) + 20*C(k,5)
    4: [1680, 420, 10],                             # etc.
    5: [30240, 10080, 560],
    6: [665280, 277200, 25200, 280],
    7: [17297280, 8648640, 1108800, 30800],
    # alpha_8: 5 terms, c_{8,0}=518918400 and c_{8,1}=302702400 from pattern
    # c_{8,2..4} from k=12..14 data
    8: [518918400, 302702400, 50450400, 2402400, 15400],
}

if __name__ == '__main__':
    import sys
    max_k = int(sys.argv[1]) if len(sys.argv) > 1 else 14

    print(f"Computing 3-cycle IP for k=2..{max_k}")
    print("=" * 60)

    all_ip = {}
    for k in range(2, max_k + 1):
        ip, t = compute_ip_bitmask(k)
        d = len(ip) - 1
        all_ip[k] = ip
        print(f"k={k}: d={d}, I_3={ip} [{t:.1f}s]")
        errors = verify_formula(k, ip, FORMULAS)
        if errors:
            for e in errors:
                print(f"  ERROR: {e}")

    print("\n=== EXTRACTING ALPHA_8 COEFFICIENTS ===")
    # Use k where d>=8, i.e., floor(2k/3)>=8, i.e., k>=12
    from fractions import Fraction

    m = 8
    data_pts = [(k, all_ip[k][m]) for k in all_ip if len(all_ip[k]) > m]
    print(f"Data for alpha_8: {data_pts}")

    # Determine unknown coefficients c_{8,2..4} from data at k=12,13,14
    # c_{8,0} = 518918400 (pattern), c_{8,1} = 302702400 (pattern)
    # Solve for c_{8,2..4}
    A0, A1 = 518918400, 302702400
    rows = []
    for k_val, alpha_val in data_pts:
        contrib0 = A0 * comb(k_val, 16) + A1 * comb(k_val, 15)
        rhs = alpha_val - contrib0
        row = [Fraction(comb(k_val, 14)), Fraction(comb(k_val, 13)), Fraction(comb(k_val, 12))]
        rows.append((row, Fraction(rhs)))
        print(f"  k={k_val}: alpha_8={alpha_val}, contrib0={contrib0}, rhs={rhs}")

    print("\n=== PREDICTING ALPHA_8 FOR LARGER K ===")
    for k_val in range(max(data_pts, key=lambda x:x[0])[0]+1, max_k+1):
        if k_val in all_ip and len(all_ip[k_val]) > 8:
            pred = sum(c * comb(k_val, 16 - j) for j, c in enumerate(FORMULAS[8]))
            actual = all_ip[k_val][8]
            status = "✓" if pred == actual else f"ERROR: pred={pred}"
            print(f"  k={k_val}: pred={pred}, actual={actual} {status}")
