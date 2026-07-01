"""n=7 category counts (B,M,K) via vectorized canon, to complete the recursion table + girth of black even-graph."""
import numpy as np
from itertools import combinations, permutations
def cats(n):
    TILES=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]; mt=len(TILES)
    tileidx={t:i for i,t in enumerate(TILES)}; TRANS=[tileidx[(n-y+1,n-x+1)] for (x,y) in TILES]
    arcs=list(combinations(range(n),2)); ai={a:i for i,a in enumerate(arcs)}; M=len(arcs)
    pw=(1<<np.arange(M)).astype(np.int64); perms=list(permutations(range(n)))
    srcpos=np.empty((len(perms),M),np.int64); oflip=np.empty((len(perms),M),np.int64)
    for pi,p in enumerate(perms):
        dest=np.empty(M,np.int64); flip=np.empty(M,np.int64)
        for e,(i,j) in enumerate(arcs):
            a,b=p[i],p[j]
            if a<b: dest[e]=ai[(a,b)]; flip[e]=0
            else: dest[e]=ai[(b,a)]; flip[e]=1
        sp=np.empty(M,np.int64); sp[dest]=np.arange(M); srcpos[pi]=sp; oflip[pi]=flip[sp]
    def canon_batch(X):
        best=np.full(X.shape[0],1<<62,np.int64)
        for pi in range(len(perms)): np.minimum(best,(X[:,srcpos[pi]]^oflip[pi])@pw,out=best)
        return best
    # build arc-bitvector per tiling: base path (k,k+1) forward; tile (xL,yL): vertices xi=n-xL,yi=n-yL
    NT=1<<mt; X=np.zeros((NT,M),np.int8)
    for e in range(n-1):  # base path arcs (e,e+1): forward => arc (e,e+1) bit 0 (i<j i->j). indices e<e+1 so ai[(e,e+1)]
        pass
    def arcbits(bits):
        v=np.zeros(M,np.int8)
        # base path k->k+1
        for k in range(n-1):
            i,j=k,k+1  # i<j, i->j => bit 0
            v[ai[(i,j)]]=0
        for t,(xL,yL) in enumerate(TILES):
            xi=n-xL; yi=n-yL  # xi<yi? xL>yL => n-xL<n-yL => xi<yi
            i,j=min(xi,yi),max(xi,yi)
            # bits[t]==0 => xi->yi ; orientation of arc (i,j): i->j is bit0
            if bits[t]==0:  # xi->yi
                v[ai[(i,j)]]= 0 if xi<yi else 1
            else:           # yi->xi
                v[ai[(i,j)]]= 1 if xi<yi else 0
        return v
    for mask in range(NT):
        bits=[(mask>>k)&1 for k in range(mt)]
        X[mask]=arcbits(bits)
    ci=canon_batch(X)
    # complement: reverse all arcs => flip all M bits
    cop=canon_batch(X^1)
    # merged node
    uniq=sorted(set(ci.tolist())); idx={c:k for k,c in enumerate(uniq)}
    tpt={}
    for r in range(NT): tpt[int(ci[r])]=int(cop[r])
    def mn(c): return min(c,tpt[c])
    gscnt={}; tau={}
    for mask in range(NT):
        bits=[(mask>>k)&1 for k in range(mt)]
        gs=all(bits[i]==bits[TRANS[i]] for i in range(mt) if TRANS[i]!=i)
        v=mn(int(ci[mask])); tau[v]=tau.get(v,0)+1; gscnt[v]=gscnt.get(v,0)+(1 if gs else 0)
    B=M_=K=0
    for v in tau:
        if gscnt[v]==0: K+=1
        elif gscnt[v]==tau[v]: B+=1
        else: M_+=1
    return B,M_,K,len(tau)
B,M,K,nodes=cats(7)
print(f"n=7: pureBLUE B={B}, mixed M={M}, pureBLACK K={K}, total nodes={nodes} (expect K=184)")
print(f"RECURSION table (B,M,K): n4=(1,1,1) n5=(3,5,2) n6=(2,10,22) n7=({B},{M},{K})")
print(f"  SC=B+M: {2},{8},{12},{B+M}  ; K(NS-merged): 1,2,22,{K}")
