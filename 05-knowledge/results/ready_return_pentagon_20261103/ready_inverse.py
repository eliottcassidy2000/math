"""Complete finite inversion of the classified ready-value images."""
from explore_returns import *
from functools import lru_cache
from verify import ratio_coeff


@lru_cache(None)
def two_start(n):
    result=set()
    for m in range(3,n):
        t=n-m-1
        words=[()] if t==0 else [w+(1,) for w in product((0,1),repeat=t-1)]
        cores=[(T,(0,)) for T in two_start(m)]
        for T in unique_cycle_classes(m):
            starts=[r for r in range(m) if unique_start_path(T,r) is not None]
            cores.extend((T,J) for J in combinations(starts,2))
        for T,J in cores:
            for word in words:
                U,r=funnel(T,set(J),word)
                paths=start_paths_at_most_two(U,r)
                require(paths is not None and len(paths)==2 and strong(U),'constructive two-start grammar')
                result.add(min(induced(U,p) for p in paths))
    return tuple(sorted(result))


def busch(n):
    q,r=divmod(n,3)
    return 3*5**(q-1) if r==0 else 5**q if r==1 else 9*5**(q-1)


def ready_values(limit):
    # Both lanes use the larger necessary parent bound; final tests are exact.
    n=3
    while busch(n)<=(limit+2)//2:n+=1
    bound=n-1;values={1:{},2:{}};counts=[]
    if limit>=1:values[1][1]=dict(T=(2,0),root=1,construction='weak transitive boundary')
    for n in range(3,bound+1):
        banks={1:[unrank(n,j) for j in range(class_counts(n)[n])],2:two_start(n)}
        counts.append(dict(n=n,one=len(banks[1]),two=len(banks[2])))
        for b,bank in banks.items():
            for T in bank:
                h,s,q,F=rooted_data(T,0)
                require(s==b,'constructed parent start count')
                child=h+F-b
                if child<=limit and child not in values[b]:
                    U,r=birth(T,0)
                    require(rooted_data(U,r)[0]==child and rooted_data(U,r)[3]-child==b,'inverse witness')
                    values[b][child]=dict(T=U,root=r,parent=T,construction='newborn')
                tables.cache_clear()
    tails=((),((1,0),),((0,1),),((1,0),(0,1)),((0,1),(1,0)))
    for tail in tails:
        k=1
        while True:
            T,r=incoming_pair_source(((1,1),)*k+tail);h,_,_,F=rooted_data(T,r)
            require(F-h==2,'exception source')
            tables.cache_clear()
            if h>limit:break
            values[2].setdefault(h,dict(T=T,root=r,construction='exception',word=((1,1),)*k+tail));k+=1
    return bound,counts,values


def main():
    V=ratio_coeff([0,0,0,0,1,-1,-1,1],[1,-6,11,-6,1],8)
    for n in range(3,9):require(len(two_start(n))==V[n],'two-start generating function')
    # An independent complete chart census identifies every generated class <=6.
    for n in range(3,7):
        chart=set()
        for bits in range(1<<((n-1)*(n-2)//2)):
            T=tiling(n,bits)
            if not strong(T):continue
            paths=start_paths_at_most_two(T,0)
            if paths is not None and len(paths)==2:chart.add(min(induced(T,p) for p in paths))
        require(chart==set(two_start(n)),'grammar vs independent chart universe')
    bound,counts,values=ready_values(127)
    require(35 not in values[1] and 35 not in values[2] and 17 not in values[2],'all-size hostile images')
    require(11 in values[1] and 11 in values[2],'both branches at the obstruction source')
    print(json.dumps(dict(status='PASS',limit=127,all_size_parent_bound=bound,parent_class_counts=counts,
        lane_values={b:sorted(d) for b,d in values.items()},witnesses=values,
        neither_ready=[h for h in range(5,128,2) if all(h not in d for d in values.values())]),indent=2,sort_keys=True))


if __name__=='__main__':main()
