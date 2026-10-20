#!/usr/bin/env python3
"""Independent reduced-word, cycle, and dyadic controls; no producer imports."""
from itertools import product, combinations
from functools import lru_cache
import json, re
from independent_checks import check, graph, literal, inversion_count


def inversion_edges(w, bit):
    return [frozenset((i,j)) for i in range(len(w)) for j in range(i+1,len(w))
            if (w[i] & (1 << bit)) and not (w[j] & (1 << bit))]


def cross(w):
    a,b = inversion_edges(w,0), inversion_edges(w,1)
    return bool(a and b) and all(e & f for e in a for f in b)


def grammar(w):
    s = ''.join(map(str,w))
    return bool(re.fullmatch(r'(1+20|2+10|1+21+|2+12+|1+3+0|2+3+0|30+1+|30+2+|312+|321+|30+|3+0|310|320)',s))


def triangles(adj):
    return [frozenset(t) for t in combinations(range(len(adj)),3)
            if all(sum(adj[i][j] for j in t) == 1 for i in t)]


def count_dp(adj):
    n=len(adj)
    @lru_cache(None)
    def extend(last, unused):
        if not unused:
            return 1
        return sum(extend(v,unused ^ (1 << v)) for v in range(n)
                   if unused & (1 << v) and adj[last][v])
    return sum(extend(v,((1 << n)-1) ^ (1 << v)) for v in range(n))


def odd_cycles(adj):
    # Each directed cycle starts at its least vertex; feedback labels are 0,1.
    totals=[0,0,0]
    n=len(adj)
    for start in range(n):
        def grow(path,used):
            if len(path)>=3 and len(path)%2 and adj[path[-1]][start]:
                x,y=0 in used,1 in used
                check(x or y,'ordinary chain cycle')
                totals[2 if x and y else 0 if x else 1] += 1
            for v in range(start+1,n):
                if v not in used and adj[path[-1]][v]:
                    grow(path+[v],used|{v})
        grow([start],{start})
    return totals


def main():
    record={}; accepted=[]; triangle_tests=0
    for m in range(2,9):
        n=0
        for w in product(range(4),repeat=m):
            if w[0]==0 or w[-1]==3:
                continue
            is_cross=cross(w)
            check(is_cross == grammar(w),('regular expression grammar',w))
            if is_cross:
                n += 1
            if m<=5:
                adj=graph(w,2,1);ts=triangles(adj)
                intercyclic=all(e & f for e,f in combinations(ts,2))
                aa,bb=inversion_edges(w,0),inversion_edges(w,1)
                check(intercyclic == (not aa or not bb or is_cross),('triangles vs inversions',w))
                triangle_tests+=1
            if is_cross and m<=7:
                adj=graph(w,2,1);ts=triangles(adj)
                common=set(range(len(adj))).intersection(*ts)
                h=count_dp(adj)
                if not common:
                    check(h in {3*2**t+1 for t in range(2,12)} | {9*2**t-3 for t in range(1,12)},('nononeapex rays',w,h))
                if m<=5:
                    check(h==literal(adj),('literal reduced grammar H',w))
        accepted.append([m,n])
    record['reduced_grammar_counts']=accepted
    record['actual_triangle_vs_inversion_controls']=triangle_tests
    ray_tests=0
    for a in range(1,5):
        families=[((1,)*a+(2,0),[3*2**(a-1)-1,1,3*2**(a-1)],6*2**a+1),
                  ((2,)*a+(1,0),[1,3*2**(a-1)-1,3*2**a-2],9*2**a-3)]
        for b in range(1,5):
            if a+b>5:
                continue
            families += [((2,)*a+(1,)+(2,)*b,[2**(b-1),2**(a-1),3*2**(a+b-1)-2**(a-1)-2**(b-1)],3*2**(a+b)+1),
                         ((2,)*a+(3,)*b+(0,),[2**(b-1),2**(a+b-1),2**(a+b)-2**(b-1)],3*2**(a+b)+1)]
        for w,c,h in families:
            adj=graph(w,2,1)
            check(odd_cycles(adj)==c,('DFS odd cycles',w))
            check(count_dp(adj)==h==1+2*sum(c),('ray H',w))
            ray_tests+=1
    record['independent_DFS_cycle_and_subset_recursion_controls']=ray_tests
    witnesses=[]
    for rows,expected in [((50,16,3,7,12,30),29),((126,124,248,240,96,192,128,19),87)]:
        adj=[[bool(rows[i] & (1 << j)) for j in range(len(rows))] for i in range(len(rows))]
        check(literal(adj)==count_dp(adj)==expected,('explicit witness',rows))
        ts=triangles(adj)
        disjoint_pairs=[[sorted(e),sorted(f)] for e,f in combinations(ts,2) if e.isdisjoint(f)]
        check(bool(disjoint_pairs)==(expected==29),('whole-count vs divisibility',expected))
        witnesses.append(dict(H=expected,disjoint_triangle_pairs=disjoint_pairs))
    record['witnesses']=witnesses
    # The finite set at chain length six independently excludes 29.
    length_six={inversion_count((1,)+mid+(0,)) for mid in product(range(2),repeat=4)}
    check(29 not in length_six,'middle order excludes29')
    record['strong_oneapex_order_seven_values']=sorted(length_six)
    density=[]
    for L in range(5,15):
        X=1 << L;values={1};orbits=set();excluded=0;fixed=0
        for m in range(2,L+2):
            for mid in product(range(2),repeat=m-2):
                w=(1,)+mid+(0,)
                dual=tuple(1-b for b in reversed(w))
                check(inversion_count(w)==inversion_count(dual),('duality',w))
                fixed += w==dual
                forbidden=m==L+1 and w[:2]==(1,1) and w[-2:]==(0,0)
                excluded += forbidden
                if not forbidden:
                    orbits.add(min(w,dual))
                h=inversion_count(w)
                if h<=X:
                    check(not forbidden,('four outer inversions',w))
                    values.add(h)
        bound=7*X//16+2**((L+1)//2-1)
        ray1={3*2**t+1 for t in range(2,L+1) if 3*2**t+1<=X}
        ray2={9*2**t-3 for t in range(1,L+1) if 9*2**t-3<=X}
        check(excluded==X//8 and fixed==2**((L+1)//2)-1,'bank sizes')
        check(len(values)<=1+len(orbits)<=bound,'orbit upper bound')
        check(len(ray1)==L-3 and len(ray2)==L-4,'ray cutoffs')
        total=len(values|ray1|ray2)
        check(total<=bound+2*L-7,'intercyclic dyadic bound')
        density.append(dict(L=L,oneapex_values=len(values),intercyclic_values=total,bound=bound+2*L-7))
    record['dyadic_controls']=density
    record['status']='PASS'
    print(json.dumps(record,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
