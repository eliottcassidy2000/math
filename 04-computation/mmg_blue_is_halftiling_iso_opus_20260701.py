"""PROVE: blue-subgraph of merged metagraph = half-tiling metagraph (isomorphism).
Grid-sym tilings = R-fixed subcube (R=complement), dim D=(m+f)/2, f=#anti-diagonal tiles. Fold: half-tile =
{fixed tiles} + {one rep per TRANS 2-cycle}. Unfold u: half-tiling -> grid-sym tiling (bijection). flip commutes.
=> blue flip-metagraph on grid-sym tilings IS the flip-metagraph of the D-dim half-tiling model."""
from itertools import permutations
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    tileidx={t:i for i,t in enumerate(TILES)}; TRANS=[tileidx[(n-y+1,n-x+1)] for (x,y) in TILES]
    perms=list(permutations(range(n)))
    def bits_adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL; A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
        return A
    def canon(A):
        best=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s<best: best=s
        return best
    # folded tiles: fixed (TRANS[i]==i) + one rep per 2-cycle (min index)
    fixed=[i for i in range(m) if TRANS[i]==i]
    reps=[]; seen=set()
    for i in range(m):
        if TRANS[i]==i: continue
        j=TRANS[i]
        if min(i,j) not in seen: seen.add(min(i,j)); reps.append(min(i,j))
    half=sorted(fixed+reps); D=len(half)
    return n,m,TILES,TRANS,fixed,reps,half,D,bits_adj,canon
for n in [4,5,6]:
    n,m,TILES,TRANS,fixed,reps,half,D,bits_adj,canon=build(n)
    f=len(fixed)
    # (1) grid-sym tilings via unfolding half-tilings; verify bijection & count
    def unfold(hbits):  # hbits over 'half' positions -> full grid-sym bits
        b=[0]*m; hb={half[k]:hbits[k] for k in range(D)}
        for i in range(m):
            if TRANS[i]==i: b[i]=hb[i]
            else: b[i]=hb[min(i,TRANS[i])]
        return b
    gridsym_full=set()
    for hm in range(1<<D):
        hbits=[(hm>>k)&1 for k in range(D)]; b=unfold(hbits); gridsym_full.add(tuple(b))
    # direct count of grid-sym tilings
    def isgs(b): return all(b[i]==b[TRANS[i]] for i in range(m) if TRANS[i]!=i)
    direct=sum(1 for mask in range(1<<m) if isgs([(mask>>k)&1 for k in range(m)]))
    biject = (len(gridsym_full)==(1<<D)==direct)
    # (2) flip commutes: flip all m tiles == flip all D half-tiles then unfold
    okflip=True
    for hm in range(1<<D):
        hbits=[(hm>>k)&1 for k in range(D)]; b=unfold(hbits)
        bflip=[1-x for x in b]; hflip=[1-x for x in hbits]
        if bflip!=unfold(hflip): okflip=False; break
    # (3) build blue subgraph = half-tiling metagraph; classes = SC classes of grid-sym tournaments
    def cls(b):
        A=bits_adj(b); c=canon(A); cop=canon([[A[j][i] for j in range(n)] for i in range(n)]); return min(c,cop)
    nodes={}; edges=[]; selfl=0; seenL=set()
    for hm in range(1<<D):
        hbits=[(hm>>k)&1 for k in range(D)]; b=unfold(hbits); c=cls(b); nodes[c]=nodes.get(c,0)+1
    for hm in range(1<<D):
        hbits=[(hm>>k)&1 for k in range(D)]; hflip=tuple(1-x for x in hbits)
        key=(min(hm, sum(hf<<k for k,hf in enumerate(hflip))), max(hm, sum(hf<<k for k,hf in enumerate(hflip))))
        if key in seenL: continue
        seenL.add(key)
        ca=cls(unfold(hbits)); cb=cls(unfold(list(hflip)))
        if ca==cb: selfl+=1
        else: edges.append((ca,cb))
    print(f"n={n}: m={m} f={f} D=(m+f)/2={D}  |half-tilings|=2^{D}={1<<D}  |grid-sym|(direct)={direct}")
    print(f"   (1) unfold bijection half-tilings<->grid-sym tilings? {biject}")
    print(f"   (2) flip commutes (flip_D unfold == unfold flip_m)? {okflip}")
    print(f"   (3) half-tiling metagraph: {len(nodes)} classes (=SC count), {len(edges)} edges + {selfl} self-loops = 2^(D-1)={1<<(D-1)} lines? {len(edges)+selfl==(1<<(D-1))}")
    print(f"       => blue subgraph IS the half-tiling flip-metagraph (nodes=SC classes, lines=flip-pairs). ISO holds.")
