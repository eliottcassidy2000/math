"""
zeckendorf_pair_dp.py — Corrected pair automaton DP for Zeckendorf tournament tilings.

The pair automaton processes binary strings of length m (no consecutive 1s)
in pairs. States:
  F (free): last emitted bit was 0, or start of string
  C (constrained): last emitted bit was 1

From F: can emit 00→F, 01→C, 10→F  (three choices)
From C: can emit 00→F, 01→C         (two choices; 10 forbidden)
11 always forbidden.

Key formulas:
  cnt_F[k+1] = 2*cnt_F[k] + cnt_C[k]
  cnt_C[k+1] = cnt_F[k] + cnt_C[k]
  with cnt_F[0]=1, cnt_C[0]=0.
  Total after p pair steps = F_{2p+2}.

  sum_H[F; k+1] = sum_H[F; k] * (1+h_left[k]) + sum_H[C; k]
  sum_H[C; k+1] = (sum_H[F; k] + sum_H[C; k]) * h_right[k]
  with sum_H[F;0]=1, sum_H[C;0]=0.
  where h_tile = 1 + 2^{range-1} (product formula approximation).

oracle-2026-05-01
"""

from math import comb, factorial
from itertools import product as iproduct

# ── Fibonacci utilities ─────────────────────────────────────────────────────

def fib(n):
    """Standard Fibonacci: F_1=1, F_2=1, F_3=2, F_4=3, F_5=5, F_6=8, F_7=13, F_8=21, ..."""
    if n <= 0: return 0
    if n == 1: return 1
    a, b = 1, 1
    for _ in range(n - 2): a, b = b, a + b
    return b


# ── Tile ordering for n-vertex tournaments ──────────────────────────────────

def tiles_for_n(n):
    """
    Return list of (left_vertex, right_vertex) tiles for n-vertex tournament.
    Tiles = backward arcs in the staircase. Vertices 0..n-1.
    Base path: 0→1→2→...→(n-1). Tile (b,a) with a-b≥2.
    Ordered matching the Fibonacci bijection: diagonals bottom-up, left to right.

    The pairing: tiles 2k and 2k+1 form pair k.
    tile 2k = left tile of pair k  (contributes when bit is 10)
    tile 2k+1 = right tile of pair k (contributes when bit is 01)
    """
    # Ordering: for range r from 2 to n-1, for b from 0 to n-r-1:
    # tile (b, b+r). This is diagonal ordering.
    tiles = []
    for r in range(2, n):
        for b in range(n - r):
            tiles.append((b, b + r))
    return tiles


def tile_range(tile):
    b, a = tile
    return a - b


def h_val(r):
    """H-value for a single backward tile of range r: 1 + 2^{r-1}."""
    return 1 + (1 << (r - 1))


# ── Zeckendorf string enumeration ───────────────────────────────────────────

def zeckendorf_strings(m):
    """Yield all binary strings of length m with no consecutive 1s."""
    if m == 0:
        yield ()
        return
    def gen(pos, last):
        if pos == m:
            yield ()
            return
        for bit in (0, 1):
            if bit == 1 and last == 1:
                continue
            for rest in gen(pos + 1, bit):
                yield (bit,) + rest
    yield from gen(0, 0)


# ── Brute-force sum of H_approx ─────────────────────────────────────────────

def brute_sum_H_approx(n):
    """Compute sum of H_approx = prod h_k over all Zeckendorf tilings at n."""
    tiles = tiles_for_n(n)
    m = len(tiles)
    h = [h_val(tile_range(t)) for t in tiles]
    total_sum = 0
    total_cnt = 0
    for s in zeckendorf_strings(m):
        H_approx = 1
        for k, bit in enumerate(s):
            if bit:
                H_approx *= h[k]
        total_sum += H_approx
        total_cnt += 1
    return total_sum, total_cnt


# ── Correct pair automaton DP ───────────────────────────────────────────────

def pair_DP(n):
    """
    O(m) pair automaton DP.
    Returns (sum_H_approx, count) over all Zeckendorf tilings.

    Transition (per pair k, using left tile 2k and right tile 2k+1):
      dp_F_new = dp_F * (1 + h_left) + dp_C     [emit 00 or 10 from F; 00 from C]
      dp_C_new = (dp_F + dp_C) * h_right          [emit 01 from F or C]

    Initial state: dp_F=1 (empty tiling, H_approx=1), dp_C=0.
    """
    tiles = tiles_for_n(n)
    m = len(tiles)
    h = [h_val(tile_range(t)) for t in tiles]

    num_pairs = m // 2
    leftover = m % 2  # handle odd m separately if needed

    # For H_approx sums
    dp_F, dp_C = 1, 0
    # For count verification
    cnt_F, cnt_C = 1, 0

    for k in range(num_pairs):
        h_left  = h[2 * k]
        h_right = h[2 * k + 1]

        # H_approx DP
        new_F = dp_F * (1 + h_left) + dp_C
        new_C = (dp_F + dp_C) * h_right
        dp_F, dp_C = new_F, new_C

        # Count DP
        new_cnt_F = 2 * cnt_F + cnt_C
        new_cnt_C = cnt_F + cnt_C
        cnt_F, cnt_C = new_cnt_F, new_cnt_C

    # Handle leftover tile if m is odd
    if leftover:
        h_left = h[m - 1]
        # Only 00 and 10 are available for a single-tile pair (no right tile)
        new_F = dp_F * (1 + h_left) + dp_C
        new_C = 0  # can't emit 01 without a right tile
        dp_F, dp_C = new_F, new_C

        new_cnt_F = 2 * cnt_F + cnt_C
        new_cnt_C = 0
        cnt_F, cnt_C = new_cnt_F, new_cnt_C

    return dp_F + dp_C, cnt_F + cnt_C


# ── H_approx closed form ────────────────────────────────────────────────────

def H_approx_closed_form(n):
    """
    The sum ∑_{T∈Z_n} H_approx(T) has a closed form via the pair DP.

    The DP shows that for m=2p pairs with h-values h_0,...,h_{m-1}:

    ∑ H_approx = dp_F + dp_C after p pair steps,
    where dp_F = dp_F(1+h_{2k}) + dp_C and dp_C = (dp_F+dp_C)*h_{2k+1}.

    The generating function is:
    ∑_{T} ∏_{k active} h_k = ∏_{pairs k} [(1+h_{left,k})(1+h_{right,k}) + h_{right,k}]
    ... (not quite a simple product due to the automaton memory)
    """
    pass  # complex, use DP directly


# ── Verification ─────────────────────────────────────────────────────────────

def verify_pair_DP(n_max=7):
    print("=" * 60)
    print("PAIR AUTOMATON DP VERIFICATION")
    print("=" * 60)

    for n in range(2, n_max + 1):
        m = comb(n - 1, 2)
        tiles = tiles_for_n(n)
        h = [h_val(tile_range(t)) for t in tiles]

        # Brute force
        brute_sum, brute_cnt = brute_sum_H_approx(n)

        # Pair DP (for even m)
        if m % 2 == 0:
            dp_sum, dp_cnt = pair_DP(n)
            cnt_ok = dp_cnt == brute_cnt
            sum_ok = dp_sum == brute_sum
            fib_ok = brute_cnt == fib(m + 2)
            print(f"\nn={n}, m={m} tiles, {m//2} pairs:")
            print(f"  Tile ranges: {[tile_range(t) for t in tiles]}")
            print(f"  h-values:   {h}")
            print(f"  Brute count = {brute_cnt}, F_{{m+2}} = F_{{{m+2}}} = {fib(m+2)} {'✓' if fib_ok else '✗'}")
            print(f"  DP count    = {dp_cnt} {'✓' if cnt_ok else '✗'}")
            print(f"  Brute sum(H_approx) = {brute_sum}")
            print(f"  DP sum(H_approx)    = {dp_sum} {'✓' if sum_ok else '✗'}")
        else:
            print(f"\nn={n}, m={m} (odd, skipping pair DP)")


# ── Full bijection table ──────────────────────────────────────────────────────

def bijection_table(n):
    """Print the N <-> T_N bijection table with H_approx for each tiling."""
    tiles = tiles_for_n(n)
    m = len(tiles)
    h = [h_val(tile_range(t)) for t in tiles]

    print(f"\nBIJECTION TABLE n={n}, m={m} tiles, F_{{m+2}}=F_{{{m+2}}}={fib(m+2)} tilings")
    print(f"Tile ranges: {[tile_range(t) for t in tiles]}")
    print(f"h-values: {h}")
    print(f"  {'N':>3} | {'bits':>20} | {'H_approx':>10}")
    N = 0
    for s in zeckendorf_strings(m):
        H_approx = 1
        for k, bit in enumerate(s):
            if bit:
                H_approx *= h[k]
        print(f"  {N:3d} | {str(list(s)):>20} | {H_approx:10d}")
        N += 1


# ── DP speedup analysis ───────────────────────────────────────────────────────

def speedup_table():
    import time

    print("\n" + "=" * 60)
    print("PAIR DP SPEEDUP: O(m) vs O(F_{m+2}) naive")
    print("=" * 60)
    print(f"  {'m':>5} | {'F_{m+2}':>14} | {'DP time (μs)':>14} | {'Speedup':>14}")
    print("-" * 60)

    for m in [6, 10, 15, 20, 25, 30, 35, 40, 45]:
        fm2 = fib(m + 2)

        # Simulate pair DP for m even with dummy h=2 for all tiles
        t0 = time.perf_counter()
        iters = 1000
        for _ in range(iters):
            dp_F, dp_C = 1, 0
            for k in range(m // 2):
                h_l, h_r = 2, 2  # dummy
                new_F = dp_F * (1 + h_l) + dp_C
                new_C = (dp_F + dp_C) * h_r
                dp_F, dp_C = new_F, new_C
        dt = (time.perf_counter() - t0) / iters * 1e6  # microseconds

        speedup = fm2 / (m // 2)  # rough: O(m/2) vs O(F_{m+2})
        print(f"  {m:5d} | {fm2:>14,} | {dt:14.2f} | {speedup:>14,.0f}×")


# ── H-ordering via Zeckendorf ─────────────────────────────────────────────────

def H_ordering(n):
    """Show how N ordering relates to H_approx ordering."""
    tiles = tiles_for_n(n)
    m = len(tiles)
    h = [h_val(tile_range(t)) for t in tiles]

    pairs = []
    N = 0
    for s in zeckendorf_strings(m):
        H_approx = 1
        for k, bit in enumerate(s):
            if bit:
                H_approx *= h[k]
        pairs.append((N, H_approx, list(s)))
        N += 1

    print(f"\nH_approx ordering at n={n}:")
    sorted_by_H = sorted(pairs, key=lambda x: x[1])
    for N, H_a, s in sorted_by_H[:10]:
        print(f"  N={N:3d}, H_approx={H_a:5d}, bits={s}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # 1. Verify pair DP correctness
    verify_pair_DP(n_max=7)

    # 2. Full bijection table at n=5
    bijection_table(5)

    # 3. H_approx ordering
    H_ordering(5)

    # 4. Speedup table
    speedup_table()

    # 5. Key theorem statement
    print("\n" + "=" * 60)
    print("KEY THEOREM: Pair automaton aggregate computation")
    print("=" * 60)
    print("""
THEOREM (Pair DP for Zeckendorf Aggregate):
  Let Z_n be the set of Zeckendorf tilings of the n-vertex staircase.
  Let m = C(n-1,2) be the number of tiles, ordered as pairs (t_0,t_1), (t_2,t_3), ...
  Let h_k = 1 + 2^{r_k - 1} where r_k is the range of tile k.

  Define dp_F[0] = 1, dp_C[0] = 0 and for k = 0, 1, ..., m/2 - 1:
    dp_F[k+1] = dp_F[k] * (1 + h_{2k}) + dp_C[k]
    dp_C[k+1] = (dp_F[k] + dp_C[k]) * h_{2k+1}

  Then:
    ∑_{T ∈ Z_n} H_approx(T) = dp_F[m/2] + dp_C[m/2]

  where H_approx(T) = ∏_{k active in T} h_k.

  Setting h_k = 1 for all k gives:
    |Z_n| = dp_F[m/2] + dp_C[m/2] = F_{m+2}

  Time complexity: O(m) vs O(F_{m+2}) = O(φ^m) brute force.
""")
