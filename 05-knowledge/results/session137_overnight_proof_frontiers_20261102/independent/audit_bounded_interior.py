"""Independent SCC compression and fixed-(h,L) affine-dyadic source audit."""
from tournament_kernel import *
from functools import lru_cache
import json

def native_data(a):
    n=len(a);N=1<<n;full=N-1
    @lru_cache(None)
    def starts(s,i):
        rest=s^(1<<i)
        return 1 if not rest else sum(starts(rest,j) for j in range(n) if rest>>j&1 and a[i][j])
    H=[1]+[sum(starts(s,i) for i in range(n) if s>>i&1) for s in range(1,N)]
    qi=[sum(starts(s,i)*H[full^s] for s in range(1,N) if s>>i&1) for i in range(n)]
    # Separate endpoint recursion for last-vertex counts.
    @lru_cache(None)
    def ends(s,i):
        rest=s^(1<<i)
        return 1 if not rest else sum(ends(rest,j) for j in range(n) if rest>>j&1 and a[j][i])
    gi=[ends(full,i) for i in range(n)]
    def state(cut):
        h=H[full];q=sum(qi[i] for i in range(n) if cut>>i&1);g=sum(gi[i] for i in range(n) if cut>>i&1)
        F=hp(append(a,cut));require((q-F+h-g)%2==0,'L integrality');L=(q-F+h-g)//2
        require(L>=0,'L nonnegative');return h,F,q,g,L
    return state

def sccs(a):
    n=len(a);reach=[]
    for i in range(n):
        seen={i};front=[i]
        while front:
            v=front.pop()
            for j in range(n):
                if a[v][j] and j not in seen:seen.add(j);front.append(j)
        reach.append(seen)
    remaining=set(range(n));groups=[]
    while remaining:
        i=min(remaining);group={j for j in remaining if j in reach[i] and i in reach[j]}
        groups.append(sorted(group));remaining-=group
    return sorted(groups,key=lambda c:-sum(a[c[0]][d[0]] for d in groups if d!=c))

def run():
    cases=positive_L=tail_cases=prefix_cases=0;bases11=set();largest_tail=0
    for n,reps in class_bank(6):
        if not n:continue
        for code in reps:
            a=matrix(n,code);old=native_data(a);components=sccs(a)
            for cut in range(1<<n):
                h,F,q,g,L=old(cut);cases+=1
                if not cut:require(F==h and L==0,'empty-cut exception');continue
                selected=[j for j,c in enumerate(components) if any(cut>>i&1 for i in c)]
                first,last=min(selected),max(selected)
                prefix={c[0] for j,c in enumerate(components) if len(c)==1 and j<first}
                suffix={c[0] for j,c in enumerate(components) if len(c)==1 and j>last}
                keep=[i for i in range(n) if i not in prefix|suffix]
                base=induced(a,keep);basecut=sum(((cut>>v)&1)<<j for j,v in enumerate(keep))
                hb,Fb,qb,gb,Lb=native_data(base)(basecut);k=len(suffix)
                require(hb==h and Lb==L,'fixed h,L compression')
                if L:
                    positive_L+=1;r=sum(len(c)==1 for c in components[first+1:last])
                    if r:require(2**(r-1)<=L,'intermediate singleton pair budget')
                    require(len(base)<=h+2+L.bit_length(),'positive-L base bound')
                else:require(len(base)<=h+1,'L0 unified base bound')
                if k:
                    tail_cases+=1;largest_tail=max(largest_tail,k)
                    require((h,F,q,g)==(h,2**k*qb+h-2*L,2**k*qb,0),'affine-dyadic output ray')
                else:require((h,F,q,g)==(hb,Fb,qb,gb),'prefix-neutral reduction')
                prefix_cases+=bool(prefix)
                if h==L==1:bases11.add((hb,Fb,qb,gb))
    require(bases11=={(1,1,3,1),(1,3,5,1)},'h1 L1 complete base bank')
    print(json.dumps({'scope':'FINITE-EXACT every tournament isomorphism class through6 and every labelled cut; independent endpoint-cover counts and SCC reconstruction',
                      'marked_sources':cases,'positive_L_sources':positive_L,'nonempty_deleted_suffixes':tail_cases,
                      'nonempty_deleted_prefixes':prefix_cases,'max_deleted_suffix_singletons':largest_tail,
                      'h1_L1_base_bank':sorted(bases11),'verdict':'PASS exact native-state reencoding, all base-order bounds, and affine-dyadic output law'},sort_keys=True,indent=2))
if __name__=='__main__':run()
