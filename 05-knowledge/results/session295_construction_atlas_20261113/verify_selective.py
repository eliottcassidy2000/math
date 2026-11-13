#!/usr/bin/env python3
"""Independent literal-path controls for the selective-port factorial filtration."""
from itertools import permutations
from math import factorial
from random import Random

def checked(condition):
    if not condition:
        raise ArithmeticError("selective-port check failed")

def hamilton_paths(rows, core_n):
    n = len(rows)
    ans = [0] * (core_n + 1)
    for word in permutations(range(n)):
        if not all((rows[a] >> b) & 1 for a, b in zip(word, word[1:])):
            continue
        runs = (word[0] < core_n) + sum(
            a >= core_n and b < core_n for a, b in zip(word, word[1:])
        )
        ans[runs] += 1
    return ans

def exceptional_count(rows, core_n):
    mask = (1 << core_n) - 1
    return sum((row & mask) not in (0, mask) for row in rows[core_n:])

def check(rows, core_n):
    q = exceptional_count(rows, core_n)
    ans = hamilton_paths(rows, core_n)
    for k in range(1, core_n + 1):
        d = factorial(max(k - 1 - q, 0)) * factorial(max(k - q, 0))
        checked(ans[k] % d == 0)
    return ans

def random_rows(n, c, q, rng):
    rows = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if i < c <= j and j - c >= q:
                continue
            u, v = (i, j) if rng.getrandbits(1) else (j, i)
            rows[u] |= 1 << v
    for j in range(c, n):
        if j - c >= q:
            if rng.getrandbits(1):
                for i in range(c):
                    rows[i] |= 1 << j
            else:
                rows[j] |= (1 << c) - 1
    return rows

def main():
    a = [56,117,129,22,100,76,141,59]
    for zmask, expect in [(0, 513), (32, 1147), (255, 513)]:
        rows = [row | ((1 << 8) if not zmask >> i & 1 else 0)
                for i, row in enumerate(a)] + [zmask]
        got = check(rows, 8)
        checked(sum(got) == expect)
        print(f"witness {zmask}: sectors {got}; H={sum(got)}")
    rng = Random(295)
    count = 0
    for n in range(3, 9):
        for c in range(1, n):
            for q in range(min(2, n-c) + 1):
                for _ in range(2):
                    check(random_rows(n, c, q, rng), c)
                    count += 1
    print(f"literal enumeration checked {count} reproducibly generated cases")

if __name__ == "__main__":
    main()

def run_dp(rows, c):
    """Subset recurrence independent of the alternating-word proof."""
    n=len(rows)
    dp=[{} for _ in range(1<<n)]
    for x in range(n):
        dp[1<<x][(x,int(x<c))]=1
    for S in range(1,1<<n):
        for (x,k),a in dp[S].items():
            avail=rows[x] & ~S
            while avail:
                b=avail & -avail
                avail-=b
                y=b.bit_length()-1
                key=(y,k+int(x>=c>y))
                nxt=dp[S|b]
                nxt[key]=nxt.get(key,0)+a
    ans=[0]*(c+1)
    for (x,k),a in dp[-1].items():
        ans[k]+=a
    return ans

def sharp():
    c=6
    n=15
    rows=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if i<c and j>=c:
                forward = i!=0 if j==n-1 else (j-c)%2==0
            else:
                forward=True
            x,y=(i,j) if forward else (j,i)
            rows[x]|=1<<y
    a=run_dp(rows,c)
    checked(a[c]==5*factorial(5)*factorial(4))
    checked(a[c]%9==0 and a[c]%27!=0)
    print(f"sharp 6-run sector: {a[c]}; v3=2")

if __name__ == "__main__":
    sharp()
