"""The 30 necklaces of the -17 cycle's shape (11,7) for 3x+1 (clock 2^11 - 3^7 = -139):
run-length orderings, carries c_w, residues mod 139, and the residue histogram over all 330 words.
Rotation acts on residues by c -> (3^{w_0} c)/2 mod 139, so residue classes are unions of rotation orbits."""
from itertools import combinations
from math import gcd
K,X,q=11,7,3; D=2**K-q**X
def carry(w):
    c=0
    for j,b in enumerate(w):
        if b: c=q*c+(1<<j)
    return c
def runs(w):
    # cyclic run-length pairs (ones, zeros) starting at a one-run boundary
    ww=list(w); 
    while not (ww[0]==1 and ww[-1]==0): ww=ww[1:]+ww[:1]
    out=[]; i=0
    while i<K:
        j=i
        while j<K and ww[j]==1: j+=1
        k=j
        while k<K and ww[k]==0: k+=1
        out.append((j-i,k-j)); i=k
    return out
neck={}
hist={}
for ones in combinations(range(K),X):
    w=[0]*K
    for p in ones: w[p]=1
    c=carry(w); hist[c%139]=hist.get(c%139,0)+1
    m=min(tuple(w[i:]+w[:i]) for i in range(K))
    if m not in neck: neck[m]=[]
    neck[m].append((c,tuple(w)))
print(f"shape ({K},{X}), clock {D}; {len(neck)} necklaces, {sum(len(v) for v in neck.values())} words")
print("necklace  runs(ones,zeros)  #runs  min carry over rotations  x_min = -c/139  integral?")
rows=[]
for m,lst in sorted(neck.items()):
    cs=sorted(c for c,_ in lst)
    integral=all(c%139==0 for c in cs)
    rows.append((len(runs(m)),''.join(map(str,m)),runs(m),cs[0],cs[0]/(-139),integral))
for r in sorted(rows):
    print(f"  {r[1]}  {r[2]}  m={r[0]}  c_min={r[3]}  x_min={r[4]:.4f}  {'INTEGRAL (the -17 cycle)' if r[5] else ''}")
print("residue histogram mod 139 over 330 words: classes hit =",len(hist),"; multiplicities:",sorted(set(hist.values())))
print("class 0 multiplicity =",hist.get(0,0),"; max nonzero-class multiplicity =",max(v for k,v in hist.items() if k!=0))
print("number of runs m over the 30 necklaces:",sorted(r[0] for r in rows))
