"""Finite head of an all-size no-ready-source theorem at H=35."""
from explore_returns import *


def main():
    records=[];chart_classes={1:{},2:{}}
    # Every rooted class with one/two starts has such a fixed-spine chart.
    for n in range(3,7):
        for bits in range(1<<((n-1)*(n-2)//2)):
            T=tiling(n,bits)
            if not strong(T):continue
            paths=start_paths_at_most_two(T,0)
            if paths is None or len(paths) not in (1,2):continue
            b=len(paths);key=min(induced(T,p) for p in paths)
            chart_classes[b][key]=n
    one_expected=sum(class_counts(6)[3:7])
    require(len(chart_classes[1])==one_expected,'complete one-start head')
    require(len(chart_classes[2])==1+5+18,'complete two-start head')
    independent={unrank(n,j) for n in range(3,7) for j in range(class_counts(n)[n])}
    require(independent==set(chart_classes[1]),'independent one-start grammar head')
    for b,classes in chart_classes.items():
        for T,n in sorted(classes.items(),key=lambda item:(item[1],item[0])):
            h,s,q,F=rooted_data(T,0)
            require(s==b and F>=h,'parent and injection bound')
            U,r=birth(T,0);ch,cs,cq,cF=rooted_data(U,r)
            require(ch==h+F-b and cF-ch==b,'actual newborn')
            records.append(dict(parent=T,n=n,lane=b,h=h,F=F,child_H=ch))
            tables.cache_clear()
    eligible=[r for r in records if r['h']<=18]
    require(all(r['child_H']!=35 for r in eligible),'no newborn ready H35')
    require(all(r['child_H']!=17 for r in eligible if r['lane']==2),'no newborn 2-ready H17')
    exceptions=[(3,-1),(6,1),(5,-1),(13,3),(14,3)]
    exceptional_values=sorted({a*2**k+b for a,b in exceptions for k in range(1,7) if a*2**k+b<=35})
    require(35 not in exceptional_values and 17 not in exceptional_values,'exception exclusion')
    # Independent SCC-factorization check for weak roots in all small charts.
    factor_checks=0
    for n in range(1,6):
        for bits in range(1<<(max(0,(n-1)*(n-2)//2))):
            T=tiling(n,bits);blocks=ordered_components(T)
            for r in range(n):
                i=next(i for i,B in enumerate(blocks) if r in B)
                predicted=1
                for j,B in enumerate(blocks):
                    U=induced(T,B)
                    if j<i:predicted*=bin_values(U,0,2)[0][2]
                    elif j==i:predicted*=rooted_data(U,B.index(r))[3]
                    else:predicted*=tables(U)[0][-1]
                require(predicted==rooted_data(T,r)[3],'SCC E(2) product')
                factor_checks+=1
            tables.cache_clear()
    print(json.dumps(dict(status='PASS',scope='all-size reduction, exact finite head',
        head_classes={b:len(c) for b,c in chart_classes.items()},eligible_parents=eligible,
        all_head_parents=records,exceptional_values_through35=exceptional_values,
        weak_factor_controls=factor_checks),indent=2,sort_keys=True))


if __name__=='__main__':main()
