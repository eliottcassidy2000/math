"""Canonical unique-start census and independent permutation/bin controls."""
import json
import importlib.util
from pathlib import Path
from itertools import permutations, product
from explore import (tables, opposite, tiling, strong, require, bin_values, induced,
                     funnel_poly, as_polynomials, evaluate, funnel,
                     tail_coefficients, seam_poly, seam, out_birth)


def weak_compositions(n,m):
    if m==1:
        yield (n,)
    else:
        for k in range(n+1):
            for tail in weak_compositions(n-k,m-1):
                yield (k,)+tail


def brute_bins(T,r,m):
    counts=[0,0,0]
    for word in permutations(range(len(T))):
        for sizes in weak_compositions(len(T),m):
            offset=0;valid=True
            for size in sizes:
                segment=word[offset:offset+size]
                if any(not(T[a] >> b & 1) for a,b in zip(segment,segment[1:])):
                    valid=False;break
                offset+=size
            if valid:
                counts[0]+=1
                if sizes[0]:
                    counts[1]+=word[sizes[0]-1]==r
                    counts[2]+=word[0]==r
    return counts


def main():
    result={'canonical_unique_start':{},'independent_checks':{},'ready31':{}}
    # Fix the unique path to 0->1->...->n-1. This is exactly one chart per
    # rooted unique-start isomorphism class, not division by guessed weights.
    for n in range(1,8):
        hist={}; all_count=0
        for bits in range(1 << (max(0,(n-1)*(n-2)//2))):
            T=tiling(n,bits)
            starts=tables(opposite(T))[1][-1][0]
            if starts==1:
                all_count+=1
                if n>1 and strong(T):
                    h=tables(T)[0][-1]
                    hist[h]=hist.get(h,0)+1
            if n<=6 and starts in (1,2) and strong(T):
                reverse_child,r=out_birth(opposite(T),0)
                child=opposite(reverse_child)
                if tables(child)[0][-1]==31 and str(starts) not in result['ready31']:
                    vals=bin_values(child,r,2)
                    require(vals[1][2]-31==starts,('ready',T,starts))
                    result['ready31'][str(starts)]={'parent':T,'root':0,'child':child,'new_root':r,'child_F':vals[1][2],'reset_H':3*31+2*starts}
            tables.cache_clear()
        expected=[0,1,1,2,5,13,34,89][n]
        expected_strong=[0,0,0,1,3,8,21,55][n]
        require(all_count==expected and sum(hist.values())==expected_strong,('Fibonacci census',n,all_count,hist))
        result['canonical_unique_start'][str(n)]={'all':all_count,'strong':sum(hist.values()),'strong_H_histogram':hist}
    brute_cases=0
    for n in range(1,5):
        for bits in range(1 << (max(0,(n-1)*(n-2)//2))):
            T=tiling(n,bits)
            for r in range(n):
                vals=bin_values(T,r,3)
                for m in range(1,4):
                    require(brute_bins(T,r,m)==[v[m] for v in vals],('permutation bins',T,r,m))
                    brute_cases+=1
    result['independent_checks']['permutation_bins']=brute_cases
    automaton=0
    for t in range(6):
        for word in product((0,1),repeat=t):
            for m in range(1,4):
                observed=0
                for assignment in product(range(m+1),repeat=t):
                    before=[i for i,c in enumerate(assignment) if c==0]
                    after=[i for i,c in enumerate(assignment) if c==1]
                    observed+=int((not before or word[before[-1]]) and (not after or not word[after[0]]))
                require(observed==tail_coefficients(word,m)[2],('automaton',word,m))
                automaton+=1
    result['independent_checks']['word_automaton']=automaton
    X=tiling(3,1);Y=tiling(4,1)
    x=as_polynomials(bin_values(X,0,len(X)))
    y=as_polynomials(bin_values(Y,2,len(Y)))
    z=seam_poly(x,y);Z,c=seam(X,0,Y,2)
    observed=bin_values(Z,c,4)
    require(all(evaluate(z[i],m)==observed[i][m] for i in range(3) for m in range(1,5)), 'polynomial seam n8')
    # Composed grammar, checked beyond the root-count specialization.
    word=(0,1,0,1)
    child,r=funnel(X,{0},word)
    predicted=funnel_poly(x,word)
    observed=bin_values(child,r,4)
    require(all(evaluate(predicted[i],m)==observed[i][m] for i in range(3) for m in range(1,5)), 'polynomial funnel n8')
    result['independent_checks']['order8_polynomial_examples']=2
    T=(26,44,57,48,34,1);r=1
    v=bin_values(T,r,2)
    require((v[0][1],v[2][1],v[1][2])==(31,5,33),'second H31 entry')
    require(strong(T) and sum(row >> r & 1 for row in T)==2,'indegree-two entry')
    require(tables(induced(T,[i for i in range(len(T)) if i!=r]))[0][-1]==5,'H31 deletion')
    result['ready31']['indegree_two']={'rows':T,'root':r,'h':31,'s':5,'q':5,'F':33,'reset_H':97}
    # Actually apply the inherited N then R reset to both entries.
    for row in result['ready31'].values():
        if 'child' in row:
            graph,root=tuple(row['child']),row['new_root']
        else:
            graph,root=tuple(row['rows']),row['root']
        first,new=out_birth(opposite(graph),root);first=opposite(first)
        second,newest=out_birth(opposite(first),new);second=opposite(second)
        require(strong(second) and tables(second)[0][-1]==row['reset_H'],'literal ternary reset')
    # Preserve the inherited all-five-Q4 orientation convention on new inputs.
    path=Path(__file__).resolve().parent.parent/'golden_inverse_defects_20261103/inverse.py'
    spec=importlib.util.spec_from_file_location('golden_reference',path)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    def observer(T,r=0):
        p,e,s=bin_values(T,r,2)
        require((s[2]-e[1])%2==0,'portal parity')
        return p[1],(s[2]-e[1])//2,(s[2]+e[1])//2
    def field(v):
        h,x,y=v
        return ((4*(h+2*x-y))%7,y%7,(h-y)%7)
    q4checks=0
    sources=[(tiling(3,1),0),(tiling(4,1),2),(T,r),(child,0)]
    for graph,root in sources:
        order=[root]+[i for i in range(len(graph)) if i!=root]
        graph=tuple(sum(((graph[u] >> w)&1) << j for j,w in enumerate(order)) for u in order)
        a=observer(graph)
        for name in ('q0','q1','q2','q4','q5'):
            block=old.literal((name,));d=len(block)
            rows=list(block)+[(row << d)|((1 << d)-1) for row in graph]
            rows[d-1]|=1 << d;rows[d]&=~(1 << (d-1))
            b=observer(tuple(rows))
            require(b==old.step(a,old.BANK[name]),('Q4',name,a,b))
            if name=='q2':
                u,v,g=field(a)
                require(field(b)==(2*v%7,(u+4*v)%7,g),'F49 q2 phase')
            q4checks+=1
    result['independent_checks']['all_five_Q4_and_q2_phase']=q4checks
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
