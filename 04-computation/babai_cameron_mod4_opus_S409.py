"""opus-2026-07-20-S409 -- BABAI-CAMERON REMARK 7.4 AT n = 1 MOD 4, AND THE 3-vs-1 SPLIT.

Babai-Cameron (EJC 7 (2000) #R38, already cited by THM-474) Remark 7.4 asks to enumerate
switching classes of tournaments in which some automorphism fixes no member, and says
"We cannot do this."

THE OWNER'S OBSERVATION, to be verified: at n = 1 (mod 4) the answer is ZERO, because every
automorphism must fix the UNIQUE EVEN MEMBER of the switching class.

WHY THERE SHOULD BE AN EVEN MEMBER AT ALL, and why the residue matters:
    sum of scores = C(n,2).  All scores even  =>  C(n,2) must be EVEN.
    n = 1 (mod 4): n odd and (n-1)/2 even  =>  C(n,2) EVEN   -> all-even tournaments CAN exist
    n = 3 (mod 4): n odd and (n-1)/2 odd   =>  C(n,2) ODD    -> all-even tournaments CANNOT exist
So the residue of n mod 4 decides whether the canonical "even representative" exists at all.
That is exactly the 3-and-7 versus 1-and-5-and-9 split: primes 3, 7, 11 (= 3 mod 4) live in
the ODD/skew world with NO even representative; 5, 9, 13 (= 1 mod 4) have one.
It is also the Paley split: Paley TOURNAMENTS need p = 3 (mod 4), Paley GRAPHS need
p = 1 (mod 4) -- the skew/symmetric dichotomy of THM-1445 IS the mod-4 dichotomy.

TESTS.
  (1) at n = 5, 9: does every switching class of tournaments contain EXACTLY ONE all-even
      tournament?  (n=5 exhaustive; n=9 by sampling switching classes)
  (2) at n = 3, 7: verify NO all-even tournament exists (C(n,2) odd).
  (3) if (1) holds, Babai-Cameron's count is 0 at n = 1 (mod 4), by the uniqueness argument.
"""
import itertools, random
import numpy as np
from collections import defaultdict

def scores(S, n):
    return [int(sum(1 for j in range(n) if j != i and S[i, j] == 1)) for i in range(n)]

def switch(S, n, W):
    T = S.copy()
    for i in range(n):
        for j in range(n):
            if i != j and (((W >> i) & 1) ^ ((W >> j) & 1)):
                T[i, j] = -S[i, j]
    return T

def rand_tournament(n, rng):
    S = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i+1, n):
            v = 1 if rng.random() < 0.5 else -1
            S[i, j] = v; S[j, i] = -v
    return S

print("="*78)
print("(0) THE PARITY OBSTRUCTION:  sum of scores = C(n,2)")
print("="*78)
print("   n   n mod 4   C(n,2)   C(n,2) parity   all-even tournament possible?")
for n in range(3, 14):
    c = n*(n-1)//2
    print(f"  {n:2d}      {n%4}      {c:4d}       {'even' if c%2==0 else 'ODD '}"
          f"            {'YES' if c%2==0 else 'NO -- impossible'}")

print()
print("="*78)
print("(1) EXHAUSTIVE n = 5:  does every switching class contain exactly ONE all-even?")
print("="*78)
n = 5
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
cls = defaultdict(list)
for m in range(1 << len(pairs)):
    S = np.zeros((n, n), dtype=np.int64)
    for b, (i, j) in enumerate(pairs):
        v = 1 if (m >> b) & 1 else -1
        S[i, j] = v; S[j, i] = -v
    # canonical switching-class key: switch so vertex 0 beats nobody it can avoid
    best = None
    for W in range(1 << n):
        T = switch(S, n, W)
        key = tuple(T[i, j] for (i, j) in pairs)
        if best is None or key < best: best = key
    cls[best].append(m)
print(f"   labelled tournaments: {1 << len(pairs)}")
print(f"   switching classes:    {len(cls)}   (expected 2^C(n-1,2) = {2**((n-1)*(n-2)//2)})")
cnt = defaultdict(int)
for key, members in cls.items():
    ev = 0
    for m in members:
        S = np.zeros((n, n), dtype=np.int64)
        for b, (i, j) in enumerate(pairs):
            v = 1 if (m >> b) & 1 else -1
            S[i, j] = v; S[j, i] = -v
        if all(s % 2 == 0 for s in scores(S, n)): ev += 1
    cnt[ev] += 1
print(f"   all-even members per switching class: {dict(sorted(cnt.items()))}")
print(f"   => {'UNIQUE even member in every class -- CONFIRMED' if list(cnt)==[1] else 'NOT unique'}")

print()
print("="*78)
print("(2) n = 9 (also 1 mod 4), by sampling switching classes")
print("="*78)
rng = random.Random(11)
n = 9
dist = defaultdict(int)
for trial in range(60):
    S = rand_tournament(n, rng)
    ev = 0
    for W in range(1 << n):
        T = switch(S, n, W)
        if all(s % 2 == 0 for s in scores(T, n)): ev += 1
    dist[ev // 2] += 1      # W and its complement give the same tournament
print(f"   all-even members per switching class over 60 random classes: {dict(sorted(dist.items()))}")
print(f"   => {'UNIQUE -- CONFIRMED at n=9' if list(dist)==[1] else 'NOT unique'}")

print()
print("="*78)
print("(3) n = 3, 7 (= 3 mod 4): verify NO all-even tournament exists")
print("="*78)
for n in (3, 7):
    pr = [(i, j) for i in range(n) for j in range(i+1, n)]
    tot = 0
    lim = 1 << len(pr)
    for m in range(min(lim, 1 << 21)):
        S = np.zeros((n, n), dtype=np.int64)
        for b, (i, j) in enumerate(pr):
            v = 1 if (m >> b) & 1 else -1
            S[i, j] = v; S[j, i] = -v
        if all(s % 2 == 0 for s in scores(S, n)): tot += 1
    print(f"   n={n}: all-even tournaments found = {tot}  (C(n,2)={n*(n-1)//2} is odd, so 0 is forced)")

print()
print("="*78)
print("(4) CONCLUSION FOR BABAI-CAMERON REMARK 7.4")
print("="*78)
print("   If each switching class at n = 1 (mod 4) has a UNIQUE all-even member E, then any")
print("   automorphism g of the switching class permutes the class and preserves the")
print("   all-even property (relabelling does not change the multiset of scores), so g(E) is")
print("   an all-even member, hence g(E) = E by uniqueness.  So g FIXES a member.")
print("   => the number of switching classes with an automorphism fixing NO member is ZERO")
print("      at n = 1 (mod 4).  Babai-Cameron's 'we cannot do this' has a clean answer on")
print("      that residue class.")
print("   At n = 3 (mod 4) there is NO even member to anchor the argument -- which is exactly")
print("      where the enumeration stays hard, and exactly where Paley TOURNAMENTS live.")
