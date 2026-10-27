"""Even balanced analogue of start_set_conditions.py (HYP-9426 (2)).  n = 2m, scores m-1 and m, v of score m (the other
role by the converse).  T - v has 2m-1 vertices and a Hamilton cycle c_0 ... c_{2m-2} (Session75 section 4: T - v is
strong for m >= 3); v has m-1 in-neighbours, U = {t : c_{t-1} in N^-(v)} has |U| = m-1 (indices mod M = 2m-1).
A rotation (c_t, ..., c_{t-1}, v) is black unless s(c_t) = m-1 (endpoint scores of a blue path sum to n-1).
If every rotation is blue, (a) t + 2 not in U and (b) exactly one of 2t - t' + 1, 2t' - t + 1 in U, as in the odd case.
Prints the solutions U (with 0 in U) of (a) and (b)."""
import sys
def search(m):
    M = 2 * m - 1; size = m - 1; sols = []
    def full(U):
        S = set(U)
        if any((t + 2) % M in S for t in U): return False
        return all(((2*t - u + 1) % M in S) + ((2*u - t + 1) % M in S) == 1 for i, t in enumerate(U) for u in U[i+1:])
    def rec(U, nxt):
        if len(U) == size:
            if full(U): sols.append(tuple(U))
            return
        for x in range(nxt, M):
            if M - x < size - len(U): break
            S = set(U)
            if (x + 2) % M in S or (x - 2) % M in S: continue
            ok = True
            for t in U:   # partial (b): both witnesses already decided and in/out wrongly
                a, b = (2*t - x + 1) % M, (2*x - t + 1) % M
                if a <= x and b <= x and ((a in S or a == x) + (b in S or b == x)) != 1: ok = False; break
            if ok: rec(U + [x], x + 1)
    rec([0], 1)
    return sols
for m in range(3, int(sys.argv[1]) + 1):
    s = search(m)
    print(f"m={m} n={2*m}: start sets satisfying (a),(b) with 0 in U: {len(s)}", s[:6])
