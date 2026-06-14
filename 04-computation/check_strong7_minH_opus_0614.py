#!/usr/bin/env python3
import subprocess, itertools
GEN="/opt/homebrew/bin/gentourng"

def gen(n, strong=False):
    args=[GEN]+(["-c"] if strong else [])+[str(n)]
    out=subprocess.run(args,capture_output=True,text=True)
    pairs=list(itertools.combinations(range(n),2))
    for line in out.stdout.split():
        bits=line.strip()
        if len(bits)!=len(pairs): continue
        beats=[0]*n
        for (i,j),b in zip(pairs,bits):
            if b=='1': beats[i]|=1<<j
            else: beats[j]|=1<<i
        yield beats

def Hdp(beats):
    n=len(beats); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1: continue
            cv=dp[mask][v]
            if not cv: continue
            succ=beats[v]&~mask&full
            w=0; s=succ
            while s:
                if s&1: dp[mask|(1<<w)][w]+=cv
                s>>=1; w+=1
    return sum(dp[full])

def Hbrute(beats):
    n=len(beats); c=0
    for p in itertools.permutations(range(n)):
        if all((beats[p[a]]>>p[a+1])&1 for a in range(n-1)): c+=1
    return c

# m=7 strong: full enumeration, find min, check 23 and 25, brute-check the minimizers
vals={}
minH=10**9; minbeats=None
for beats in gen(7,strong=True):
    h=Hdp(beats)
    vals[h]=vals.get(h,0)+1
    if h<minH: minH=h; minbeats=beats
print("m=7 strong: distinct H values (sorted):", sorted(vals))
print("min strong-H at m=7:", minH)
print("count of strong tournaments with each H (first few):", dict(sorted(vals.items())[:6]))
print("is 23 a strong-H value at m=7?", 23 in vals)
print("is 25 a strong-H value at m=7?", 25 in vals)
# brute-force confirm the minimizer
print("brute-force H of a minimizer:", Hbrute(minbeats), "(dp said", Hdp(minbeats),")")

# Also: does 23 appear as strong-H at ANY m<=8? It IS a strong-H at m=6 per earlier run.
for m in [6,7,8]:
    s=set()
    for beats in gen(m,strong=True):
        s.add(Hdp(beats))
    print(f"m={m}: 23 in strong-H? {23 in s}; min strong-H = {min(s)}; formula m^2-5m+9 = {m*m-5*m+9}")
