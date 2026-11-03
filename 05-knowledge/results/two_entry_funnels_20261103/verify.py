"""Exact controls for the two-entry classification; run with Python 3, no packages."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import importlib.util
import json

from research import *


def start_paths_at_most_two(T, r):
    """Return all paths if their number is 0,1,2; None means more than two."""
    if len(T)==1:
        return [(r,)]
    keep=[v for v in range(len(T)) if v!=r]
    U=induced(T,keep); blocks=ordered_components(U); core=blocks[0]
    entries=[j for j,v in enumerate(core) if T[r] >> keep[v] & 1]
    if not entries:
        return []
    if len(entries)>2 or any(len(b)>1 for b in blocks[1:]):
        return None
    C=induced(U,core)
    if len(entries)==1:
        inner=start_paths_at_most_two(C,entries[0])
        if inner is None:return None
    else:
        inner=[unique_start_path(C,j) for j in entries]
        if any(p is None for p in inner):return None
    tail=tuple(keep[b[0]] for b in blocks[1:])
    return [(r,)+tuple(keep[core[j]] for j in p)+tail for p in inner]


def ratio_coeff(numerator, denominator, bound):
    a=[]
    for n in range(bound+1):
        a.append((numerator[n] if n<len(numerator) else 0)
                 -sum(denominator[j]*a[n-j] for j in range(1,min(n,len(denominator)-1)+1)))
    return a


def exceptional_word(word):
    k=0
    while k<len(word) and word[k]==(1,1):k+=1
    return k>=1 and tuple(word[k:]) in ((),((1,0),),((0,1),),((1,0),(0,1)),((0,1),(1,0)))


def joint_pair_word(word):
    """Positive two-bin automaton, first letter 11."""
    require(word[0]==(1,1),'first letter')
    g1=d1=g2=d2=0;ap=bp=1
    for a,b in word[1:]:
        g1,d1,g2,d2=(a*(g1+d1)+g1+(1-b)*ap,
                       (1-a)*(g1+d1)+d1+(1-b)*(1-ap),
                       b*(g2+d2)+g2+(1-a)*bp,
                       (1-b)*(g2+d2)+d2+(1-a)*(1-bp))
        ap,bp=a,b
    return ap+bp+g1+g2


def rooted_data(T,r):
    H,ends=tables(T);_,starts=tables(opposite(T))
    full=(1<<len(T))-1
    F=0;sub=full
    while sub:
        if sub >> r & 1:F+=ends[sub][r]*H[full^sub]
        sub=(sub-1)&full
    return H[-1],starts[-1][r],H[full^(1<<r)],F


def exceptional_signature(T,r):
    J=[v for v in range(len(T)) if T[v] >> r & 1]
    if len(J)!=2:return None
    a,b=J
    if T[b] >> a & 1:a,b=b,a
    keep=[v for v in range(len(T)) if v not in (a,b,r)]
    C=induced(T,keep);blocks=ordered_components(C)
    if any(len(c)!=1 for c in blocks):return None
    order=[keep[c[0]] for c in blocks]
    return tuple(((T[v]>>a)&1,(T[v]>>b)&1) for v in order)


def main():
    D=[1,-3,1];D2=[1,-6,11,-6,1]
    V=ratio_coeff([0,0,0,0,1,-1,-1,1],D2,20)
    P=ratio_coeff([0,0,0,1,0,-1],D,12)
    Q=ratio_coeff([0,0,0,1,-1],D,12)
    # Independent endpoint DP for every fixed-spine chart through six.
    census=[];decisions=0;ready=Counter();reset_checks=0
    for n in range(1,8):
        one=two=strong_two=0;hc_split=Counter()
        for bits in range(1 << (max(0,(n-1)*(n-2)//2))):
            T=tiling(n,bits);paths=start_paths_at_most_two(T,0)
            one+=paths is not None and len(paths)==1
            two+=paths is not None and len(paths)==2
            st=n>=3 and strong(T)
            if paths is not None and len(paths)==2 and st:
                strong_two+=1
                hc_split[sum((T[p[-1]]>>p[0])&1 for p in paths)]+=1
            if n<=6:
                actual=tables(opposite(T))[1][-1][0]
                require((None if actual>2 else actual)==(None if paths is None else len(paths)),('recognizer DP',T,actual,paths))
                if paths is not None:
                    require(len(set(paths))==len(paths) and all(p[0]==0 and sorted(p)==list(range(n)) and all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in paths),('actual paths',T,paths))
                decisions+=1
                if st:
                    for r in range(n):
                        h,s,q,F=rooted_data(T,r);gap=F-h
                        if gap in (1,2):
                            d=sum(T[v]>>r&1 for v in range(n));delta=gap-q+s
                            require(d==1 or (gap==2 and d==2 and delta==2 and q==s),('ready classification',T,r,(h,s,q,F),d))
                            if d==1:
                                keep=[v for v in range(n) if v!=r];parent=induced(T,keep)
                                a=next(v for v in range(n) if T[v]>>r&1)
                                require(strong(parent) and tables(opposite(parent))[1][-1][keep.index(a)]==gap,('ready parent',T,r))
                            else:
                                require(exceptional_word(exceptional_signature(T,r)),('exceptional classification',T,r))
                            ready[n,gap,d]+=1
                    # One source root per chart: complete strong chart universe.
                    if n<=5:
                        R,x=reset(T,0)
                        gaps=[rooted_data(R,r)[3]-rooted_data(R,r)[0] for r in range(n+2)]
                        require(n==3 or not any(g in (1,2) for g in gaps),('reset trap',T,gaps))
                        reset_checks+=1
            tables.cache_clear()
        require(strong_two==2*V[n],('V coefficients',n,strong_two,V[n]))
        require(hc_split.get(1,0)==2*(Q[n-1] if n else 0),('cycle-one stratum',n,hc_split))
        census.append(dict(n=n,one_charts=one,two_charts=two,strong_two_charts=strong_two,cycle_split=dict(hc_split)))
    # Classical unique-cycle grammar independently reached through rooted unrank.
    entry_classes=[]
    for n in range(3,10):
        cores=unique_cycle_classes(n);pairs=adjacent=0;hist=Counter()
        for T in cores:
            J=[v for v in range(n) if unique_start_path(T,v) is not None]
            require(len(J) in (2,3),('entry cardinality',T,J))
            mins=[v for v in range(n) if T[v].bit_count()==1]
            require(n==3 or set(J)==set(mins)|{(v+1)%n for v in mins},('successor entry law',T,J,mins))
            require(n==3 or all(induced(T,list(range(k,n))+list(range(k)))!=T for k in range(1,n)),('rigid core',T))
            hist[len(J)]+=1
            pairs+=len(list(combinations(J,2)))
            adjacent+=sum((a-b)%n in (1,n-1) for a,b in combinations(J,2))
        if n==3:pairs//=3;adjacent//=3
        require((pairs,adjacent)==(P[n],Q[n]),('pair GF',n,pairs,adjacent))
        entry_classes.append(dict(n=n,cores=len(cores),unique_start_hist=dict(hist),pairs=pairs,adjacent_pairs=adjacent))
    # Exhaustive word language and independent subset counting of joint covers.
    word_checks=direct_words=0
    for m in range(1,8):
        for tail in product(((0,0),(0,1),(1,0),(1,1)),repeat=m-1):
            word=((1,1),)+tail;C=joint_pair_word(word)
            require((C==2 and word[-1]!=(0,0))==exceptional_word(word),('five languages',word,C))
            if m<=4:
                T,r=incoming_pair_source(word);U=induced(T,list(range(r)))
                require(joint_values(U,{0,1},1)==C,('independent C2 DP',word,C))
                direct_words+=1
            tables.cache_clear();word_checks+=1
    tails=[(),((1,0),),((0,1),),((1,0),(0,1)),((0,1),(1,0))]
    formulas=[(3,-1),(6,1),(5,-1),(13,3),(14,3)]
    family_checks=[]
    for tail,(a,b) in zip(tails,formulas):
        vals=[]
        for k in range(1,7):
            T,r=incoming_pair_source(((1,1),)*k+tail)
            h,s,q,F=rooted_data(T,r)
            require(strong(T) and h==a*2**k+b and F-h==2 and s==q,('exception family',tail,k,(h,s,q,F)))
            vals.append(h);tables.cache_clear()
        family_checks.append(dict(tail=tail,A=a,B=b,H_k1_to6=vals))
    # The exceptional n=3 boundary to the reset theorem is intentional.
    C3=tiling(3,1);R,x=reset(C3,0)
    boundary=[rooted_data(R,r)[3]-rooted_data(R,r)[0] for r in range(5)]
    require(2 in boundary,('C3 boundary',boundary))
    result=dict(status='PASS',census=census,independent_DP_decisions=decisions,
                entry_classes=entry_classes,V_coefficients=V,word_checks=word_checks,
                independent_joint_words=direct_words,exceptional_families=family_checks,
                ready_chart_occurrences=[dict(n=n,gap=b,indegree=d,count=v) for (n,b,d),v in sorted(ready.items())],
                reset_chart_controls=reset_checks,C3_reset_gaps=boundary)
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
