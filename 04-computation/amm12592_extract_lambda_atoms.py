import sys
from math import comb
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation')
from amm12592_block_profile_ilp_thm3006 import feasible_ilp
PROF={8:[3,4,4,4,3,2,1,0], 16:[8,8,9,9,10,10,9,8,7,6,5,4,3,2,1,0]}
CAT=[1,1,2,5,14,42,132,429]; BALLOT=[comb(2*n,n-1) for n in range(1,8)]
print("Catalan", CAT); print("ballot C(2n,n-1)", BALLOT); print()
for m,a in PROF.items():
    ok,w=feasible_ilp(m,a); assert ok
    print(f"=== m={m} ===")
    for k in range(m):
        e=[2*w[(k,i)]-comb(a[k],i) for i in range(a[k]+1)]
        # Lam_k(w) = sum_i e_i w^i (1-w)^{a_k-i}
        L=[0]*(a[k]+1)
        for i,c in enumerate(e):
            if not c: continue
            d=a[k]-i
            for j in range(d+1):
                L[i+j]+=c*comb(d,j)*(-1)**j
        while len(L)>1 and L[-1]==0: L.pop()
        pos=sum(x for x in L if x>0); neg=-sum(x for x in L if x<0)
        print(f"  k={k:2d} a={a[k]:2d}  Lam_k = {L}   (pos={pos}, neg={neg})")
