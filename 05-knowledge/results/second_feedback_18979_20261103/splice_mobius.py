"""All-different connector recurrence, with exact square and parity controls."""
from itertools import permutations,product
from functools import lru_cache
from math import factorial
import sys,json
from support import HERE,PREV,need,save
sys.path.insert(0,str(PREV))
from square_splices import square,components,bridges,BASE

def injective(menus):
    menus=tuple(map(frozenset,menus));n=len(menus)
    @lru_cache(None)
    def meet(mask):
        b=mask&-mask;i=b.bit_length()-1
        return menus[i] if mask==b else menus[i]&meet(mask^b)
    @lru_cache(None)
    def f(mask):
        if not mask:return 1
        b=mask&-mask;rest=mask^b;sub=rest;out=0
        while True:
            B=sub|b;k=B.bit_count()
            out+=(-1)**(k-1)*factorial(k-1)*len(meet(B))*f(mask^B)
            if not sub:break
            sub=(sub-1)&rest
        return out
    return f((1<<n)-1)

def matching_count(menus):
    # Independent subset matching DP over connector names.
    universe=sorted(set().union(*map(set,menus)));D={0:1}
    for row in menus:
        N={}
        for mask,c in D.items():
            for j,v in enumerate(universe):
                if v in row and not mask>>j&1:N[mask|1<<j]=N.get(mask|1<<j,0)+c
        D=N
    return sum(D.values())

def square_count(V):
    C=components(V);blocks=[tuple(p for p in permutations(c) if all(square(a+b) for a,b in zip(p,p[1:]))) for c in C]
    total=0;guards=0;nonzero=0
    for order in permutations(range(len(C))):
        for chosen in product(*(blocks[i] for i in order)):
            menus=[bridges(A[-1],B[0],set(V)) for A,B in zip(chosen,chosen[1:])]
            a=injective(menus);need(a==matching_count(menus),'independent shared-connector matching DP')
            total+=a;guards+=1;nonzero+=a>0
    need(total%2==0,'global path reversal pairs')
    return {'anchors':V,'component_count':len(C),'ordered_endpoint_guards':guards,'nonzero_guards':nonzero,'oriented_paths':total,'paths_modulo_reversal':total//2}

def main():
    cases=[square_count(BASE),square_count(BASE+(8,))]
    need([x['paths_modulo_reversal'] for x in cases]==[9,17],'complete inherited square catalogues')
    controls=[]
    for r in range(2,6):
        U=[tuple(x for x in range(1<<r) if x.bit_count()%2==parity) for parity in (0,1)]
        M=[[{x for x in u if x>>i&1} for i in range(r)] for u in U]
        for mask in range(1,(1<<r)-1):
            cols=[i for i in range(r) if mask>>i&1]
            sizes=[len(set.intersection(*(z[i] for i in cols))) for z in M]
            need(sizes[0]==sizes[1],'all proper joint intersections equal')
        counts=[injective(z) for z in M];ind=[matching_count(z) for z in M]
        need(counts==ind and abs(counts[0]-counts[1])==factorial(r-1),'sharp top-intersection witness')
        controls.append({'slots':r,'universe_sizes':list(map(len,U)),'proper_intersections_equal':True,'injective_counts':counts,'difference_absolute':factorial(r-1)})
    save('splice_mobius.json',{'status':'FINITE_EXACT_WITH_ALL_SLOT_RECURRENCE','squares':cases,'parity_hostiles':controls,
        'scope':'minimum square connectors on the specified anchors; arbitrary guarded slot menus in the theorem; not a classification of all square-sum graphs or tournaments'})
    print('PASS square splice recurrence and sharp joint-depth parity hostiles',flush=True)
if __name__=='__main__':main()
