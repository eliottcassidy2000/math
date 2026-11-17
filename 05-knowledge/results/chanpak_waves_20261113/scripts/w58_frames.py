"""Wave 58: two-hub counts (transfer DP) and their frame decomposition (wave 6), checked on the unit-switch witnesses."""
def count_order(A, B, ab):
    dp = {(0, 0, 0): 1}
    for i in range(len(A)):
        nd = {}
        for (s1, s2, s3), c in dp.items():
            t = ((1 if A[i] == 0 else 2), s2, s3); nd[t] = nd.get(t, 0) + c
            if s2 != 0 or A[i] == 1:
                t = (s1, 1 if B[i] == 0 else 2, s3); nd[t] = nd.get(t, 0) + c
            if s3 != 0 or B[i] == 1:
                t = (s1, s2, 1); nd[t] = nd.get(t, 0) + c
        dp = nd
    return sum(c for (s1, s2, s3), c in dp.items() if 2 not in (s1, s2, s3) and (s2 != 0 or ab))
def H(X, Y, e):
    """e = 1 iff u1 -> u2"""
    return count_order(X, Y, e) + count_order(Y, X, 1 - e)
"""Frame decomposition of two-hub counts: H = sum over hub orders and boundary data of prod_i (#segments available to i)."""
import itertools
def frames(X, Y, e):
    m = len(X); out = []
    for order in (0, 1):
        A, B, ab = (X, Y, e) if order == 0 else (Y, X, 1 - e)   # first hub word A, second B, ab = first->second
        # S0: empty or last r0 with A[r0]==0 ; S1: empty (needs ab) or [a1..b1] with A[a1]==1, B[b1]==0 ; S2: empty or first c2 with B[c2]==1
        R0 = [None] + [i for i in range(m) if A[i] == 0]
        S1 = ([None] if ab else []) + [(a, b) for a in range(m) for b in range(a, m) if A[a] == 1 and B[b] == 0]
        C2 = [None] + [i for i in range(m) if B[i] == 1]
        for r0 in R0:
            for s1 in S1:
                for c2 in C2:
                    w = 1; ok = True; alpha = beta = 0
                    for i in range(m):
                        av = []
                        in0 = r0 is not None and i <= r0; in1 = s1 is not None and s1[0] <= i <= s1[1]; in2 = c2 is not None and i >= c2
                        fixed = [r0 == i, s1 is not None and i in (s1[0], s1[1]), c2 == i]
                        if sum(fixed) > 1: ok = False; break
                        if any(fixed): continue
                        k = in0 + in1 + in2
                        if k == 0: ok = False; break
                        if k == 2: alpha += 1
                        if k == 3: beta += 1
                    if not ok: continue
                    # segments that are declared nonempty must contain their boundary vertices: r0 in S0, a1,b1 in S1, c2 in S2 (fixed); fine
                    out.append((order, r0, s1, c2, alpha, beta))
    return out
def check(X, Y, e):
    fr = frames(X, Y, e); tot = sum(2**a * 3**b for *_, a, b in fr)
    return tot, H(X, Y, e), sum(1 for f in fr if f[-1] > 0)
for (x, y) in (("101000", "011010"), ("1011101100", "0001111011"), ("1100100010", "0010000111"), ("100010", "001011")):
    X = [int(c) for c in x]; Y = [int(c) for c in y]
    for e in (0, 1):
        tot, h, n3 = check(X, Y, e)
        print(f"x={x} y={y} e={e}: frame sum {tot}, H {h}, frames with a triply covered vertex: {n3}")
