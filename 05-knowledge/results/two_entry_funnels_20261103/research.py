"""Two-entry funnel construction experiments; all counters are exact."""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import importlib.util
import json
import sys

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'polarized_path_calculus_20261103'
sys.path.insert(0,str(PREV))
from explore import (tables, opposite, induced, tiling, strong, require,
                     bin_values, joint_values, funnel, out_birth,
                     tail_coefficients, transitive)
from recognize import unique_start_path, unrank, class_counts, ordered_components


def birth(T,r):
    U,x=out_birth(opposite(T),r)
    return opposite(U),x


def reset(T,r):
    U,x=birth(T,r)
    V,y=birth(U,x)
    return V,x


def cyclic_code(T):
    n=len(T)
    return min(induced(T,list(range(k,n))+list(range(k))) for k in range(n))


def unique_cycle_classes(n):
    if n==3:
        return [tiling(3,1)]
    return sorted({cyclic_code(unrank(n,j)) for j in range(class_counts(n)[n])})


def hamilton_paths(T,root=None):
    for word in permutations(range(len(T))):
        if root is not None and word[0]!=root:
            continue
        if all(T[u] >> v & 1 for u,v in zip(word,word[1:])):
            yield word


def rooted_code(T,r):
    return min(induced(T,p) for p in hamilton_paths(T,r))


def small_core_bank():
    bank={1:[],2:[]}
    for n in range(3,6):
        for T in unique_cycle_classes(n):
            starts=[r for r in range(n) if unique_start_path(T,r) is not None]
            for r in starts:
                bank[1].append((T,(r,)))
            for ports in combinations(starts,2):
                bank[2].append((T,ports))
        seen=set()
        for mask in range(1 << ((n-1)*(n-2)//2)):
            T=tiling(n,mask)
            if strong(T) and tables(opposite(T))[1][-1][0]==2:
                code=rooted_code(T,0)
                if code not in seen:
                    seen.add(code);bank[2].append((code,(0,)))
    return bank


def core_data(T,J):
    values=[bin_values(T,r,2) for r in J]
    return dict(h=values[0][0][1],p2=values[0][0][2],
                e=sum(v[1][1] for v in values),F=sum(v[1][2] for v in values),
                G=sum(v[2][2] for v in values),b=sum(v[2][1] for v in values),
                C=joint_values(T,set(J),1))


def ready_count(data,word):
    A=tail_coefficients(word,2)[0]
    _,D,B=tail_coefficients(word,1)
    return A*(data['G']+data['p2'])+B*data['h']-D*data['e']-data['C']-data['F']-data['b']


def incoming_pair_source(word):
    """Ports a=0 -> b=1, ordered chain 2.., root last loses only to a,b."""
    m=len(word);n=m+3;r=n-1
    T=[0]*n;T[0]|=1 << 1
    for i,(A,B) in enumerate(word):
        v=i+2
        for port,bit in enumerate((A,B)):
            if bit:T[v]|=1 << port
            else:T[port]|=1 << v
        for j in range(i+1,m):T[v]|=1 << (j+2)
        T[r]|=1 << v
    T[0]|=1 << r;T[1]|=1 << r
    return tuple(T),r


def main():
    bank=small_core_bank()
    print('core bank', {b:len(rows) for b,rows in bank.items()},flush=True)
    words=[w for k in range(4) for w in product((0,1),repeat=k)]
    hits={1:[],2:[]}
    families={1:{},2:{}}
    for b,cores in bank.items():
        for T,J in cores:
            data=core_data(T,J)
            for left in words:
                for right in words:
                    for bit in (0,1):
                        sample=[left+(bit,)*t+right for t in (1,2,3)]
                        if not sample[0][-1]:continue
                        vals=[ready_count(data,w) for w in sample]
                        require((vals[1]-vals[0])%2==0,'integral dyadic coefficient')
                        aa=(vals[1]-vals[0])//2;bb=vals[0]-2*aa
                        require(vals[2]==8*aa+bb,('dyadic template',data,left,bit,right,vals))
                        if aa>0:
                            key=(aa,bb)
                            rec=dict(core=T,ports=J,left=left,bit=bit,right=right,values=vals)
                            if key not in families[b] or len(T)+len(left)+len(right)<len(families[b][key]['core'])+len(families[b][key]['left'])+len(families[b][key]['right']):
                                families[b][key]=rec
                        if vals==[97,181,349]:
                            hits[b].append(dict(core=T,ports=J,left=left,bit=bit,right=right,values=vals))
    print('uniform candidate counts', {b:len(v) for b,v in hits.items()},flush=True)
    print(json.dumps({str(b):v[:4] for b,v in hits.items()},indent=2),flush=True)
    common=set(families[1])&set(families[2])
    from functools import lru_cache
    @lru_cache(None)
    def depth(key):
        if key not in common:return -1
        a,b=key
        return 1+min(depth((3*a,3*b+2)),depth((3*a,3*b+4)))
    best=sorted(common,key=lambda key:(-depth(key),2*key[0]+key[1],key))
    print('common affine banks',len(common),'best',[(key,depth(key)) for key in best[:12]],flush=True)
    result={'bank':{b:len(v) for b,v in bank.items()},'hits':hits,
            'family_counts':{b:len(v) for b,v in families.items()},'common_count':len(common),
            'selected':[dict(A=key[0],B=key[1],depth=depth(key),lane1=families[1][key],lane2=families[2][key]) for key in best[:12]]}
    (HERE/'template_search.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
