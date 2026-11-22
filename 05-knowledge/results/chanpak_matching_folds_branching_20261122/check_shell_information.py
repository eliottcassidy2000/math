#!/usr/bin/env python3
"""Exact ordinary counts and Z-polynomials for two strong self-converse probes."""

from itertools import combinations, permutations, product

U = (58, 4, 1, 6, 14, 30)
V = (38, 8, 2, 5, 15, 30)
ALPHA_U = (2, 5, 0, 4, 3, 1)
ALPHA_V = (3, 5, 4, 0, 2, 1)

def matrix(rows):
    return [[bool(row >> j & 1) for j in range(len(rows))] for row in rows]

def check_strong_anti(a, sigma):
    n = len(a)
    for i in range(n):
        if sigma[sigma[i]] != i or not all(
            (a[i][j] != a[j][i]) and (a[i][j] == a[sigma[j]][sigma[i]])
            for j in range(n) if i != j
        ):
            raise ArithmeticError("not a tournament with the stated anti-involution")
        discovered = {i}
        while True:
            nxt = discovered | {k for j in discovered for k in range(n) if a[j][k]}
            if nxt == discovered:
                break
            discovered = nxt
        if len(discovered) != n:
            raise ArithmeticError("not strong")

def all_h(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    out = [1] + [0]*((1 << n) - 1)
    for s in range(1, 1 << n):
        for u, val in enumerate(dp[s]):
            if not val:
                continue
            out[s] += val
            for v in range(n):
                if not (s >> v & 1) and a[u][v]:
                    dp[s | (1 << v)][v] += val
    return out

def specializations(h, r_max):
    old = [1] + [0]*(len(h)-1)
    for _ in range(r_max):
        new = [0]*len(h)
        for s in range(len(h)):
            part = s
            while True:
                new[s] += h[part] * old[s ^ part]
                if not part:
                    break
                part = (part-1) & s
        yield new[-1]
        old = new

def inflate_pair(a, sigma):
    n = 2 * len(a)
    return ([[a[i//2][j//2] if i//2 != j//2 else i < j
              for j in range(n)] for i in range(n)],
            [2*sigma[i//2] + 1-i%2 for i in range(n)])

def mark(a, sigma):
    n = len(a)
    choices = [i for i, row in enumerate(a) if sum(row) == 2]
    v, = choices
    w = sigma[v]
    return ([row+[i==v, i!=w] for i, row in enumerate(a)] +
            [[j!=v for j in range(n)]+[0,0], [j==w for j in range(n)]+[1,0]],
            list(sigma)+[n+1,n])

def reverse_shell(a, sigma):
    n = len(a)
    return ([row+[False, True] for row in a] +
            [[True]*n+[False, False], [False]*n+[True, False]],
            list(sigma)+[n+1,n])

def center_counts(a, sigma):
    n = len(a)
    pairs = [(i, sigma[i]) for i in range(n) if i < sigma[i]]
    vals = []
    for picks in product(*pairs):
        b = ([row+[i not in picks] for i,row in enumerate(a)] +
             [[i in picks for i in range(n)]+[False]])
        check_strong_anti(b, list(sigma)+[n])
        vals.append(all_h(b)[-1])
    return vals

def cycle_sets(a, length):
    found = []
    for s in combinations(range(len(a)),length):
        for perm in permutations(s[1:]):
            cyc = s[:1]+perm
            if all(a[i][j] for i,j in zip(cyc, cyc[1:]+cyc[:1])):
                found.append(frozenset(cyc))
    return found

def main():
    for rows, sigma, expected, marked, centers, cycles in (
        (U,ALPHA_U,[5080,14288,6944,1136,64,1],333,
         [33,49,51,67,93,105,111,123], (4,4,0)),
        (V,ALPHA_V,[4024,11528,6608,1136,64,1],325,
         [27,41,65,69,75,87,115,129], (4,2,1)),
    ):
        a = matrix(rows)
        check_strong_anti(a,sigma)
        base = list(specializations(all_h(a),3))
        if base != [r**6+8*r**4+8*r**2 for r in range(1,4)]:
            raise ArithmeticError(("base polynomial",base))
        tri = cycle_sets(a,3)
        if (len(tri),len(cycle_sets(a,5)),
                sum(not (x & y) for x,y in combinations(tri,2))) != cycles:
            raise ArithmeticError("packing witness")
        shell, sig = a, sigma
        shell_counts = []
        for i in range(4):
            check_strong_anti(shell,sig)
            shell_counts.append(all_h(shell)[-1])
            if i < 3:
                shell,sig = reverse_shell(shell,sig)
        if shell_counts != [17,241,4001,80929]:
            raise ArithmeticError(("shell tower", shell_counts))
        if sorted(center_counts(a,sigma)) != centers:
            raise ArithmeticError("center spectra")
        print(f"shell {shell_counts}; centered spectrum {centers}")
        b,sigma_b = inflate_pair(a,sigma)
        check_strong_anti(b,sigma_b)
        nums = list(specializations(all_h(b),6))
        want = [sum(c*r**(2*i) for i,c in enumerate(expected,1))
                for r in range(1,7)]
        if nums != want:
            raise ArithmeticError(("inflated polynomial",nums))
        c,sigma_c = mark(a,sigma)
        check_strong_anti(c,sigma_c)
        if all_h(c)[-1] != marked:
            raise ArithmeticError("marked probes")
        print(f"inflated cover evaluations {nums}; marked count {marked}")

if __name__ == "__main__":
    main()
