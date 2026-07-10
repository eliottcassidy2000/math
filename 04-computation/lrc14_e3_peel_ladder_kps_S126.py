# kind-pasteur-2026-07-09-S126
# Numerical certificate for the E3-side Freiman-stability ladder (LRCSchurPeel.lean).
# Confirms every Lean theorem's claim + the empirical capstone dist <= deficit.
from itertools import combinations
from math import comb

def E3(S):
    Sset = set(S)
    return sum(1 for a in S for b in S if a + b in Sset)

def deficit(S):
    return comb(len(S), 2) - E3(S)

def repCount(S, m):
    return sum(1 for a in S for b in S if a + b == m)

def totalPeelCost(S):
    S = sorted(S); tot = 0
    while S:
        m = S[-1]; k = len(S)
        tot += (k - 1) - repCount(S, m)   # peel cost of the top rung
        S = S[:-1]
    return tot

def reflection_symmetric(S, m):   # every a in S, a != m, has m-a in S
    Sset = set(S)
    return all((m - a) in Sset for a in S if a != m)

def all_peels_full(S):
    S = sorted(S)
    while len(S) >= 2:
        if repCount(S, S[-1]) != len(S) - 1:
            return False
        S = S[:-1]
    return True

def is_dilated(S):
    S = sorted(S); d = S[0]
    return all(S[i] == (i + 1) * d for i in range(len(S)))

def dist_to_dilated(S):
    k = len(S); Sset = set(S); best = k
    for d in range(1, max(S) + 1):
        best = min(best, k - len(Sset & set(d * (i + 1) for i in range(k))))
    return best

print("=== E3-side peel ladder: Lean-theorem claims + empirical capstone ===")
print("(all sets of positive integers, 0 not in S)\n")
for k in [3, 4, 5, 6, 7]:
    N = 2 * k + 3
    n = 0
    rungA = rungB = peelrec = mono = peelcost = accum = fullpeel = base = capstone = 0
    for S in combinations(range(1, N + 1), k):
        S = list(S); n += 1
        m = max(S); Se = [x for x in S if x != m]
        # Rung A: E3 S = E3(S\{max}) + repCount(max)
        if E3(S) == E3(Se) + repCount(S, m): rungA += 1
        # Rung B: repCount(max) <= k-1
        if repCount(S, m) <= k - 1: rungB += 1
        # deficit peel: deficit S + repCount = deficit(erase) + (k-1)
        if deficit(S) + repCount(S, m) == deficit(Se) + (k - 1): peelrec += 1
        # monotonicity: deficit(erase) <= deficit S
        if deficit(Se) <= deficit(S): mono += 1
        # peelCost_le_deficit: (k-1)-repCount <= deficit S
        if (k - 1) - repCount(S, m) <= deficit(S): peelcost += 1
        # capstone accumulation: deficit = totalPeelCost
        if deficit(S) == totalPeelCost(S): accum += 1
        # full-peel characterization: repCount(max)=k-1 <=> reflection symmetric
        if (repCount(S, m) == k - 1) == reflection_symmetric(S, m): fullpeel += 1
        # base: deficit=0 <=> dilated
        if (deficit(S) == 0) == is_dilated(S): base += 1
        # EMPIRICAL capstone: dist_to_dilated <= deficit
        if dist_to_dilated(S) <= deficit(S): capstone += 1
    def tag(c): return "OK " if c == n else "FAIL"
    print(f"k={k} ({n:5d} sets):")
    print(f"  [PROVEN] Rung A  E3 peel            {tag(rungA)} {rungA}/{n}")
    print(f"  [PROVEN] Rung B  repCount<=k-1       {tag(rungB)} {rungB}/{n}")
    print(f"  [PROVEN] deficit peel recursion      {tag(peelrec)} {peelrec}/{n}")
    print(f"  [PROVEN] deficit monotone (peel)     {tag(mono)} {mono}/{n}")
    print(f"  [PROVEN] peelCost <= deficit         {tag(peelcost)} {peelcost}/{n}")
    print(f"  [PROVEN] CAPSTONE deficit=totalPeel  {tag(accum)} {accum}/{n}")
    print(f"  [PROVEN] full-peel <=> reflection    {tag(fullpeel)} {fullpeel}/{n}")
    print(f"  [PROVEN] base deficit=0 <=> dilated  {tag(base)} {base}/{n}")
    print(f"  [EMPIR ] capstone dist <= deficit    {tag(capstone)} {capstone}/{n}")

# Bonus: (reflection-symmetric at EVERY peel) <=> dilated -- corollary of accumulation + base
print("\n=== corollary: (all peels full) <=> dilated ===")
for k in [3, 4, 5, 6]:
    n = ok = 0
    for S in combinations(range(1, 3 * k), k):
        S = list(S); n += 1
        if all_peels_full(S) == is_dilated(S): ok += 1
    print(f"  k={k}: {'OK ' if ok==n else 'FAIL'} {ok}/{n}")
