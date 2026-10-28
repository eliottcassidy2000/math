"""Independent finite-domain checks of the signed affine carrier theorem."""
from itertools import product, permutations
from pathlib import Path
import json

def solve(n,edges,L,U):
    coeff=[None]*n; coeff[0]=(1,0)
    for _ in range(n):
        for u,v,s,w in edges:
            if coeff[u] is not None and coeff[v] is None:
                e,t=coeff[u]; coeff[v]=(s*e,s*t+w)
            if coeff[v] is not None and coeff[u] is None:
                e,t=coeff[v]; coeff[u]=(s*e,s*(t-w))
    if any(c is None for c in coeff): raise ValueError('connected graph required')
    lo=max(L-t if e==1 else t-U for e,t in coeff)
    hi=min(U-t if e==1 else t-L for e,t in coeff)
    forbidden=set()
    for u,v,s,w in edges:
        eu,tu=coeff[u]; ev,tv=coeff[v]
        a,b=ev-s*eu,s*tu+w-tv
        if not a:
            if b: return []
        else:
            if b%a: return []
            z=b//a; lo=max(lo,z); hi=min(hi,z)
    for i in range(n):
        ei,ti=coeff[i]
        for j in range(i):
            ej,tj=coeff[j]
            if ei==ej:
                if ti==tj: return []
            elif (tj-ti)%(ei-ej)==0: forbidden.add((tj-ti)//(ei-ej))
    return [tuple(e*z+t for e,t in coeff) for z in range(lo,hi+1) if z not in forbidden]

def brute(n,edges,L,U):
    return sorted(x for x in permutations(range(L,U+1),n) if all(x[v]==s*x[u]+w for u,v,s,w in edges))

def main():
    checks=0
    settings=list(product((-1,1),range(-3,4)))
    for skeleton in [[(0,1),(1,2)],[(0,1),(1,2),(2,0)],[(0,1),(1,2),(2,3)],[(0,1),(0,2),(0,3)]]:
        n=1+max(max(e) for e in skeleton)
        for params in product(settings,repeat=len(skeleton)):
            edges=[(u,v,s,w) for (u,v),(s,w) in zip(skeleton,params)]
            a=sorted(solve(n,edges,0,4)); b=brute(n,edges,0,4)
            if a!=b: raise RuntimeError((edges,a,b))
            checks+=1
    square_edges=[(0,1,-1,4),(1,2,-1,9),(2,3,-1,16),(3,4,-1,25),(4,0,-1,16)]
    square=solve(5,square_edges,1,15)
    if square!=[(1,3,6,10,15)]: raise RuntimeError('square cycle')
    odd_hostile=solve(3,[(0,1,-1,4),(1,2,-1,9),(2,0,-1,16)],0,20)
    if odd_hostile: raise RuntimeError('odd cycle integrality')
    # Complete sign and edge-length choices on the five-vertex path.
    graceful=[]
    for lengths in permutations(range(1,5)):
        for signs in product((-1,1),repeat=4):
            edges=[(i,i+1,1,s*w) for i,(s,w) in enumerate(zip(signs,lengths))]
            graceful+=solve(5,edges,0,4)
    direct=[x for x in permutations(range(5)) if {abs(x[i+1]-x[i]) for i in range(4)}=={1,2,3,4}]
    if sorted(graceful)!=sorted(direct): raise RuntimeError('graceful completeness')
    out={'status':'FINITE-EXACT','full_signed_systems_checked':checks,'domain':[0,4],'skeletons':'three-vertex path and triangle, four-vertex path and star','edge_signs':[-1,1],'edge_constants':[-3,-2,-1,0,1,2,3],'square_cycle_solution':square,'odd_square_cycle_no_integral_root':not odd_hostile,'graceful_five_path_labelings':len(graceful)}
    Path(__file__).with_name('signed.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))

if __name__=='__main__': main()
