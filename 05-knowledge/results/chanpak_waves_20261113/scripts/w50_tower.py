"""Wave 50: the owner's mod-3 tower, checked, and the place of 7 and 21.
row(n) = 2*floor(n/3) + [3 does not divide n]: multiples of 3 get short rows, 3k+1 | 3k+2 share a long row (left | right)."""
def row(n): return 2 * (n // 3) + (1 if n % 3 else 0)
N = 3000
# (1) odd n: the accelerated image (3n+1)/2 is 2 mod 3 (right side of a long row) and sits in row n
ok1 = all(((3 * n + 1) // 2) % 3 == 2 and row((3 * n + 1) // 2) == n for n in range(1, N, 2))
# (2) the three odd families: n = 6j+1, 6j+3, 6j+5 start at rows 4j+1, 4j+2, 4j+3 and span 2j, 2j+1, 2j+2 rows
ok2 = all(row(6*j+r) == 4*j + {1: 1, 3: 2, 5: 3}[r] and row((3*(6*j+r)+1)//2) - row(6*j+r) == 2*j + {1: 0, 3: 1, 5: 2}[r]
          for j in range(400) for r in (1, 3, 5))
# (3) doubling chains of multiples of 3: row(3o*2^k) = 2o*2^k ; doubling links n <- 2n span 2k (n = 3k+1) or 2k+2 (n = 3k+2)
ok3 = all(row(3 * o * 2 ** k) == 2 * o * 2 ** k for o in range(1, 50, 2) for k in range(8))
ok4 = all(row(2 * n) - row(n) == (2 * (n // 3) if n % 3 == 1 else 2 * (n // 3) + 2) for n in range(1, N) if n % 3)
print("odd images land right of row n:", ok1, "| family starts/spans:", ok2, "| 3-chains:", ok3, "| doubling spans:", ok4)
# (4) half-counts c = (H-1)/2: inverse accelerated tree on counts; five states
def state(H):
    c = (H - 1) // 2
    if c % 3 == 0: return 'sterile'
    if c % 3 == 1: return 'live unary'
    return {2: 'fork, odd child live unary', 5: 'fork, odd child sterile', 8: 'fork, odd child fork'}[c % 9]
print("states by H mod 18:", {r: state(r) for r in range(1, 19, 2)})
first_sterile_fork = min(H for H in range(3, 200, 2) if state(H) == 'fork, odd child sterile')
print("first sterile fork:", first_sterile_fork, "children:", 2 * first_sterile_fork - 1, (2 * first_sterile_fork - 1) // 3)
print("even children of the holes:", 2 * 7 - 1, 2 * 21 - 1)
