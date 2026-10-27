"""Necessary conditions on the start set U (THM-4819, prescribed black endpoint lemma).
If every Hamilton path of a regular T (n = 2m+1) ending at v is blue, fix a Hamilton cycle c_0 ... c_{2m-1} of T - v and
let U = {pos(u) + 1 : u in N^-(v)} (|U| = m).  For t in U the reversal sigma_t of (c_t, ..., c_{t-1}, v) is an
anti-automorphism with sigma_t(v) = c_t, so N^+(c_t) = sigma_t(N^-(v)): c_t -> c_x iff 2t - x + 1 in U (x != t).
Hence (a) t + 2 not in U (the cycle arc c_{t-1} -> c_t), and (b) for distinct t, t' in U exactly one of 2t - t' + 1,
2t' - t + 1 lies in U.  Prints the m-subsets of Z_{2m} satisfying (a) and (b), up to rotation."""
import sys
def search(m):
    M = 2 * m; sols = []
    def ok_with(U, x):
        if (x + 2) % M in U or (x - 2) % M in U: return False
        return True
    def full_check(U):
        S = set(U)
        for t in U:
            if (t + 2) % M in S: return False
        for i, t in enumerate(U):
            for t2 in U[i+1:]:
                if ((2*t - t2 + 1) % M in S) + ((2*t2 - t + 1) % M in S) != 1: return False
        return True
    def partial_ok(U):
        # (b) can only be refuted once both witnesses are decided; check pairs whose witnesses are < current max
        S = set(U); mx = U[-1]
        for i, t in enumerate(U):
            for t2 in U[i+1:]:
                a, b = (2*t - t2 + 1) % M, (2*t2 - t + 1) % M
                ina, inb = a in S, b in S
                da, db = a <= mx, b <= mx
                if da and db and ina + inb != 1: return False
                if ina and inb: return False
        return True
    def rec(U, nxt):
        if len(U) == m:
            if full_check(U): sols.append(tuple(U))
            return
        for x in range(nxt, M):
            if M - x < m - len(U): break
            if not ok_with(U, x): continue
            U2 = U + [x]
            if partial_ok(U2): rec(U2, x + 1)
    rec([0], 1)   # rotate so that 0 is in U
    return sols
for m in range(2, int(sys.argv[1]) + 1):
    s = search(m)
    print(f"m={m} n={2*m+1}: start sets satisfying (a),(b) with 0 in U: {len(s)}", s[:4])
