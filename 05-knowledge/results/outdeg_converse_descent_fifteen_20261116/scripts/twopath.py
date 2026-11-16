#!/usr/bin/env python3
"""twopath.py: the two-path lemma needed for the out-degree-1 law.  For every strong tournament C of order 3..8 and every
ordered pair of distinct vertices (x, x'), test
  (L1) one-sided: C has vertex-disjoint directed paths P_a starting at x and P_b ending at x' that together cover C;
  (L2) two-sided: (L1) holds for (x, x') or for (x', x).
Also the simpler sufficient conditions: x' in the terminal component of C - x, or x in the initial component of C - x'.
Reports the number of failures of each."""
import subprocess, itertools
exec(open('nearjoin_check.py').read().split('# (1)')[0])   # parse, classes, ham_all, starts, strong, sub, scc_of
def ends(adj, t):
    # paths covering a subset S ending at t: from ham_all's dp
    n = len(adj); dp, H = ham_all(adj); N = 1 << n
    return lambda S: dp[S][t]
def two_path_one_sided(C, x, xp):
    n = len(C); N = 1 << n; dp, H = ham_all(C)       # dp[S][last]: paths covering S ending at last
    st, dps = starts(C, x)                             # st[S]: paths covering S starting at x
    full = N - 1
    for W in range(1, N):
        if not (W >> x) & 1 or (W >> xp) & 1: continue
        if st[W] == 0: continue
        Wc = full & ~W
        if dp[Wc][xp] > 0: return True
    return False
def term_comp(C, v):   # terminal component of C (all vertices reach it) membership of v: v reaches nothing outside its SCC... simpler: v in terminal SCC iff reach(v) == scc(v)
    n = len(C)
    def reach(s):
        r = {s}; f = [s]
        while f:
            a = f.pop()
            for b in range(n):
                if C[a][b] and b not in r: r.add(b); f.append(b)
        return r
    return reach(v) == scc_of(C, v)
def init_comp(C, v):
    n = len(C); Ct = [[C[j][i] for j in range(n)] for i in range(n)]
    return term_comp(Ct, v)
tot = 0; fail1 = 0; fail2 = 0; fail_simple = 0; examples = []
for n in range(3, 9):
    for C in classes(n):
        if not strong(C): continue
        for x in range(n):
            for xp in range(n):
                if x == xp: continue
                tot += 1
                a = two_path_one_sided(C, x, xp)
                if not a:
                    fail1 += 1
                    b = two_path_one_sided(C, xp, x)
                    if not b:
                        fail2 += 1
                        if len(examples) < 5: examples.append((n, x, xp))
                # simple sufficient conditions
                Cmx = sub(C, [v for v in range(n) if v != x]); idx = [v for v in range(n) if v != x].index(xp)
                Cmxp = sub(C, [v for v in range(n) if v != xp]); idy = [v for v in range(n) if v != xp].index(x)
                if not (term_comp(Cmx, idx) or init_comp(Cmxp, idy)): fail_simple += 1
print(f"strong tournaments orders 3-8, ordered pairs (x, x'): {tot}; one-sided lemma (L1) fails {fail1}; two-sided lemma (L2) fails {fail2}; the two simple sufficient conditions both fail {fail_simple} times")
print("examples of (L2) failures:", examples)
