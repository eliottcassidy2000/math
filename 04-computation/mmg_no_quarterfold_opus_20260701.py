from itertools import permutations
def build(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; m=len(TILES)
    ti={t:i for i,t in enumerate(TILES)}; TRANS=[ti[(n-y+1,n-x+1)] for (x,y) in TILES]
    perms=list(permutations(range(n)))
    def adj(bits):
        A=[[0]*n for _ in range(n)]
        for k in range(n-1): A[k][k+1]=1
        for i,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL; A[xi][yi]=1 if bits[i]==0 else 0; A[yi][xi]=1-A[xi][yi]
        return A
    def canon(A):
        b=None
        for p in perms:
            s=tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if b is None or s<b: b=s
        return b
    return n,m,TILES,TRANS,adj,canon
for n in [4,5,6]:
    n,m,TILES,TRANS,adj,canon=build(n)
    allsc=True
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        if not all(bits[i]==bits[TRANS[i]] for i in range(m) if TRANS[i]!=i): continue  # grid-sym only
        A=adj(bits); c=canon(A); cop=canon([[A[j][i] for j in range(n)] for i in range(n)])
        if c!=cop: allsc=False; break
    print(f"  n={n}: every grid-sym (half-tiling) tournament is self-complementary (c==c^op)? {allsc}")
print("  => the folding involution (R=complement) is order-2; its fixed set = the SC world = H_n;")
print("     within H_n every tournament is complement-FIXED, so 'fold by complement' again fixes EVERYTHING")
print("     => NO proper quarter-tiling. The fold M_n->H_n is a SINGLE terminal step (corrects S19 'tower').")
