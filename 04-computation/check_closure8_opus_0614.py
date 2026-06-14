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
# strong atoms m=1,3..8
atoms={1:{1}}
for m in range(3,9):
    atoms[m]=set(Hdp(b) for b in gen(m,strong=True))
def closure(N):
    reach={0:{1}}
    for v in range(1,N+1):
        s=set()
        for c in [1]+list(range(3,v+1)):
            if c not in atoms: continue
            if (v-c) in reach:
                for p in reach[v-c]:
                    for a in atoms[c]: s.add(p*a)
        if s: reach[v]=s
    return reach.get(N,set())
R8=set(Hdp(b) for b in gen(8))
cl8_full=closure(8)
print("closure(8) WITH m=8 atoms == R_8 ?", cl8_full==R8)
print("min strong-H m=8:", min(atoms[8]), " (8 strong components singleton give factor presence)")
print("|R_8|=",len(R8),"|closure8_full|=",len(cl8_full))
print("symmetric diff:", sorted(R8 ^ cl8_full))
