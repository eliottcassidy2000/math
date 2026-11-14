# Needs the compiled census: cc -O2 -o runwhite runwhite.c (same directory).
# Greedy elimination of maximal blocks from a prefix-closed Z-free block set (cap 4), keeping the census
# "apparently cofinite": no misses in [2^19, 2^21]. Each candidate removal drops one maximal word (a word that is
# not a proper prefix of another kept word). Reports the round-by-round choice.
import subprocess, itertools
X=2097152; K1=19
def run(words):
    out=subprocess.run(['nice','-n','5','./runwhite',str(X),'4',str(K1)]+words,capture_output=True,text=True).stdout
    # "... | misses M, largest L, misses >= 2^K1: Z"
    tail=out.split('|')[1]
    miss=int(tail.split('misses ')[1].split(',')[0]); largest=int(tail.split('largest ')[1].split(',')[0]); late=int(tail.split(': ')[1])
    return miss,largest,late
def closure(ws):
    s=set()
    for w in ws:
        for i in range(1,len(w)+1): s.add(w[:i])
    return s
def maximal(s): return [w for w in s if not any(v!=w and v.startswith(w) for v in s)]
cur=closure(['OZZZ','OZZO','OZOZ','OOZZ','OOOZ','OOOO'])   # the 14-block set (with empty block and prefixes)
print("start: %d nonempty blocks + empty"%len(cur), sorted(cur), run(maximal(cur)), flush=True)
while True:
    best=None
    for w in sorted(maximal(cur)):
        if w=='O': continue
        cand=cur-{w}
        r=run(maximal(cand))
        print("   try remove %s -> misses %d largest %d late %d"%(w,*r), flush=True)
        if r[2]==0 and (best is None or r[1]<best[1][1]): best=(w,r)
    if best is None: print("no removal keeps [2^19,2^21] miss-free; stop", flush=True); break
    cur=cur-{best[0]}
    print("REMOVE %s: now %d nonempty blocks + empty: %s  (misses %d, largest %d)"%(best[0],len(cur),sorted(cur),best[1][0],best[1][1]), flush=True)
