from itertools import combinations
from math import comb
import json

def need(x,m):
    if not x:raise RuntimeError(m)
def matchings(labels):
    if not labels:yield ();return
    a=labels[0]
    for b in labels[1:]:
        for tail in matchings(tuple(x for x in labels if x not in (a,b))):
            yield tuple(sorted(((a,b),)+tail))
def audit(V,act,expected):
    index={x:i for i,x in enumerate(V)};s=[index[act(v)] for v in V]
    f=sum(s[i]==i for i in range(15));b=(15-f)//2
    need((f,b)==expected,'actual fixed-point character')
    E=list(combinations(range(15),2));arcs={};seen=set();orbits=[]
    for e in E:
        if e in seen:continue
        a,c=e;mapped=tuple(sorted((s[a],s[c])))
        if mapped==e:
            arcs[e]=True;seen.add(e)
            orbits.append((1,1 if (s[a],s[c])==e else -1))
        else:
            arcs[e]=True;arcs[mapped]=s[a]<s[c]
            seen.update((e,mapped));orbits.append((2,1))
    def edge(a,c):return arcs[(a,c)] if a<c else not arcs[(c,a)]
    defects=[e for e in E if edge(*e)!=edge(s[e[0]],s[e[1]])]
    need(len(defects)==b,'sharp literal orientation certificate')
    counts={k:orbits.count(k) for k in ((1,1),(1,-1),(2,1))}
    m=b*(f+b-1);exponent=comb(f,2)+b+m
    need(counts=={(1,1):comb(f,2),(1,-1):b,(2,1):m},'pair orbit structure')
    need(sum((1<<exponent)*comb(m,j) for j in range(m+1))==1<<105,'complete enumerator total')
    return {'permutation':s,'fixed_vertices':f,'swapped_pairs':b,'pair_orbits':[{'length':k[0],'sign':k[1],'count':v} for k,v in counts.items()],'minimum_defect':b,'polynomial':f'2^{exponent} z^{b} (1+z^2)^{m}','one_minimizer_rows':[sum((1<<j) for j in range(15) if j!=i and edge(i,j)) for i in range(15)],'forced_defect_pairs':defects}
D=list(combinations(range(6),2));M=sorted(matchings(tuple(range(6))))
need(len(D)==len(M)==15 and len(set(M))==15,'exact carriers')
swap=lambda x:1-x if x<2 else x
pairact=lambda v:tuple(sorted(map(swap,v)))
matchact=lambda v:tuple(sorted(pairact(e) for e in v))
a=audit(D,pairact,(7,4));b=audit(M,matchact,(3,6))
triple=lambda x:x^1
c=audit(D,lambda v:tuple(sorted(map(triple,v))),(3,6))
need(a['minimum_defect']!=b['minimum_defect'],'equal size does not preserve action')
print(json.dumps({'status':'PASS','duads_under_transposition':a,'synthemes_under_transposition':b,'duads_under_triple_transposition':c,'scope':'The first two are different S6 carriers; one matched outer-cycle character is not a proof of the full outer isomorphism.'},sort_keys=True,indent=2))
