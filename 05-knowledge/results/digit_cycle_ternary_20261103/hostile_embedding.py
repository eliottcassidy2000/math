"""A literal 54-group on the shared triple structure, with orientation guards."""
import importlib.util
from itertools import combinations
from collections import Counter
from support import ROOT,need,save
def compose(p,q):return tuple(p[q[i]] for i in range(len(p)))
def power(p,n):
    out=tuple(range(len(p)))
    for _ in range(n):out=compose(p,out)
    return out
def order(p):
    q=tuple(range(len(p)))
    for n in range(1,200):
        q=compose(p,q)
        if q==tuple(range(len(p))):return n
    raise RuntimeError('bounded group order')
def group(gens):
    identity=tuple(range(len(gens[0])));seen={identity};todo=[identity]
    while todo:
        p=todo.pop()
        for g in gens:
            q=compose(g,p)
            if q not in seen:seen.add(q);todo.append(q)
    return seen
def tournament(flip):
    rows=[0]*12
    for u,v in combinations(range(12),2):
        b,s=divmod(u,4);c,t=divmod(v,4)
        arc=(c-b)%3==1 if b!=c else s==0 or (s>0 and (t-s)%3==1)
        if flip and b==c==2:arc=not arc
        rows[u if arc else v]|=1<<(v if arc else u)
    return tuple(rows)
def main():
    path=ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('actual_hostile',path);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    A,B=tournament(False),tournament(True);GA=set(h.iso_maps(A,A));GB=set(h.iso_maps(B,B))
    g=next(p for p in sorted(GA) if order(p)==9);v=next(i for i in range(12) if len({power(g,j)[i] for j in range(9)})==9)
    orbit=[power(g,j)[v] for j in range(9)];a=(v//4)*4;outer=[power(g,j)[a] for j in range(3)]
    need(len(set(orbit+outer))==12,'actual nine-plus-three labelled source frame')
    D=[None]*12
    for j,x in enumerate(orbit):D[x]=orbit[2*j%9]
    for j,x in enumerate(outer):D[x]=outer[2*j%3]
    D=tuple(D);G=group([g,D]);P=group([g,power(D,2)])
    triples=lambda T:{t for t in combinations(range(12),3) if all(sum(T[u]>>v&1 for v in t)==1 for u in t)}
    C=triples(A);need(C==triples(B) and len(C)==67,'same actual labelled cyclic triples')
    need(len(G)==54 and len(P)==27 and order(D)==6,'literal full and ternary affine groups')
    need(compose(D,g)==compose(power(g,2),D),'doubling/translation group relation')
    need(all({tuple(sorted(p[v] for v in t)) for t in C}==C for p in G),'every54 permutation preserves the actual common triples')
    need(G&GA==P and len(G&GB)==9,'orientation-preserving intersections retain the hostile distinction')
    hist=lambda X:dict(sorted(Counter(order(p) for p in X).items()))
    need(hist(GB)=={1:1,3:26} and hist(P)=={1:1,3:8,9:18},'two order27 groups have different element-order spectra')
    original=h.hostile_checks()
    save('hostile_embedding.json',{'status':'FINITE_EXACT_ACTUAL_GROUP_EMBEDDING','translation_permutation':g,
      'doubling_permutation':D,'inner_clock_orbit':orbit,'outer_orbit':outer,'common_triples':sorted(C),
      'affine_group_order':len(G),'ternary_subgroup_order':len(P),'intersection_with81':len(G&GA),'intersection_with27':len(G&GB),
      'order_histograms':{'actual81':hist(GA),'actual27':hist(GB),'affine54':hist(G),'affine27':hist(P)},'original_control':original,
      'scope':'54 preserves the cyclic-triple relation; only27 preserve the first tournament, only9 preserve the second'})
    print('PASS actual54 on common triples; oriented intersections27/9; affine27 differs from hostile27',flush=True)
if __name__=='__main__':main()
