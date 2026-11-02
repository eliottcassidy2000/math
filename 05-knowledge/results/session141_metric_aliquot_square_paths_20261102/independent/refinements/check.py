from itertools import combinations
from collections import Counter
from math import isqrt, comb
import json
from functools import reduce
from math import gcd
from fractions import Fraction

def run(n):
    vertices=tuple(range(1,n+1))
    edges=tuple((i,j) for i in vertices for j in vertices if i<j and isqrt(i+j)**2==i+j)
    adj={i:tuple(j for j in vertices if j!=i and isqrt(i+j)**2==i+j) for i in vertices}
    paths=[]
    calls=0
    def visit(path,used):
        nonlocal calls
        calls+=1
        if len(path)==n: paths.append(tuple(path)); return
        for j in adj[path[-1]]:
            b=1<<(j-1)
            if not used&b: visit(path+[j],used|b)
    for v in vertices: visit([v],1<<(v-1))
    arcs=[frozenset(zip(p,p[1:])) for p in paths]
    compatible=lambda a,b: not any((v,u) in b for u,v in a)
    moments=Counter({0:2**len(edges)})
    profiles={}
    subsets=[]
    def walk(chosen,candidates,union):
        if chosen:
            k=len(chosen)
            profiles.setdefault(k,Counter())[len(union)]+=1
            moments[k]+=2**(len(edges)-len(union))
            subsets.append((chosen,len(union)))
        for pos,j in enumerate(candidates):
            walk(chosen+(j,),[k for k in candidates[pos+1:] if compatible(arcs[j],arcs[k])],union|arcs[j])
    walk((),list(range(len(paths))),frozenset())
    coeff=[sum((-1)**(j-k)*comb(j,k)*a for j,a in moments.items() if j>=k) for k in range(max(moments)+1)]
    out={'n':n,'guard_edges':len(edges),'directed_paths':len(paths),'dfs_calls':calls,'subset_union_profiles':profiles,'binomial_moments':moments,'guard_polynomial':coeff,'free_tournament_edges':comb(n,2)-len(edges),'paths':paths,'compatible_subsets':subsets}
    if sum(coeff)!=2**len(edges) or any(x<0 for x in coeff): raise RuntimeError('polynomial gate')
    return out


def require(c,label):
    if not c:raise RuntimeError(label)

def main():
    reports=[]
    for n in (25,26):
        r=run(n);paths=r['paths'];lookup={p:i for i,p in enumerate(paths)}
        reversal={i:lookup[p[::-1]] for i,p in enumerate(paths)}
        subsets={tuple(q):u for q,u in r['compatible_subsets']}
        for q,u in subsets.items():
            opposite=tuple(sorted(reversal[i] for i in q))
            require(opposite!=q,'nonempty compatible subset fixed by reversal')
            require(subsets[opposite]==u,'reversal changed compatibility or union count')
        N=comb(n,2); umax=max(subsets.values()); guaranteed=N-umax+1
        coeff=[a*2**r['free_tournament_edges'] for a in r['guard_polynomial']]
        moments=[r['binomial_moments'][j]*2**r['free_tournament_edges'] for j in range(max(r['binomial_moments'])+1)]
        content=reduce(gcd,coeff);moment_content=reduce(gcd,moments)
        require(content==moment_content,'content not invariant under integral basis shift')
        require(content==2**(276 if n==25 else 296),'exact content')
        require(all(c%2**guaranteed==0 for c in coeff),'reversal divisibility')
        mean=Fraction(sum(i*c for i,c in enumerate(coeff)),2**N)
        positive=Fraction(sum(coeff[1:]),2**N)
        require(mean==(Fraction(5,2**22) if n==25 else Fraction(3,2**22)),'mean')
        require(positive==(Fraction(19,2**24) if n==25 else Fraction(351,2**29)),'positive probability')
        selected=next(q for q in subsets if len(q)==(2 if n==25 else 3) and all(paths[i][0]==18 for i in q))
        bank=[paths[i] for i in selected]; arcsets=[set(zip(p,p[1:])) for p in bank]
        union=set.union(*arcsets);common=set.intersection(*arcsets)
        b=n-len(common);J=len(union)-len(common)
        outgoing={v:sum(u==v for u,w in union) for v in range(1,n+1)}
        z=sum(k==0 for k in outgoing.values());fixed=sum(max(0,k-1) for k in outgoing.values())
        require(fixed==len(union)-n+z,'fixed direction deletion budget')
        require((b,J,z,fixed)==((3,4,0,1) if n==25 else (5,8,0,3)),'cross-lane examples')
        reversed_union={(v,u) for u,v in union}
        reverse_degrees=[sum(u==v for u,w in reversed_union) for v in range(1,n+1)]
        require(reverse_degrees.count(0)==1 and sum(max(0,d-1) for d in reverse_degrees)==fixed+1,'reversal fixed-arrow loss')
        reports.append({'n':n,'mean':str(mean),'positive_probability':str(positive),'common_block_budget':{'b':b,'J':J,'sinks':z,'cycle_rank':J-b+1,'reorientable_defect':max(0,J-b),'fixed_arrow_defect':fixed,'reversed_fixed_arrow_defect':fixed+1},'max_compatible_union':umax,'reversal_bound_exponent':guaranteed,'exact_content_exponent':content.bit_length()-1,'binomial_moments':moments,'normalized_coefficients':[c//content for c in coeff],'paired_compatible_nonempty_subsets':len(subsets)//2})
    require(reports[1]['exact_content_exponent']-reports[0]['exact_content_exponent']==20,'content exponent increase')
    bank=[(1,2,3,4,5),(3,4,1,2,5)];union=set().union(*(set(zip(p,p[1:])) for p in bank))
    require(not any((v,u) in union for u,v in union),'same-terminal hostile compatibility')
    degrees=[sum(u==v for u,w in union) for v in range(1,6)]
    require((len(union)-5,sum(max(0,d-1) for d in degrees),degrees.count(0))==(1,2,1),'same-terminal budgets must differ')
    hostile=run(1)
    require(hostile['guard_polynomial']==[0,1],'n1 hostile')
    print(json.dumps({'status':'PASS','rows':reports,'n1_reversal_fixed_hostile_polynomial':[0,1],'same_terminal_hostile':{'paths':bank,'reorientable_budget':1,'fixed_arrow_budget':2},'scope':'unpruned complete square-path banks n25 and26; all compatibility cliques and reversal pairs; n1 hostile'},indent=2,sort_keys=True))
if __name__=='__main__':main()
