#!/usr/bin/env python3
"""Independent generic-isomorphism self-controls and complete strong spectra."""
from audit import *

def graph(n,word,base=2):
    a=[[0]*n for _ in range(n)]
    for u,v in combinations(range(n),2):a[u][v]=a[v][u]=word%base;word//=base
    return a

def graph_key(a):
    n=len(a)
    return min(tuple(a[p[i]][p[j]] for i,j in combinations(range(n),2)) for p in permutations(range(n)))

def strongly_connected(a):
    n=len(a);reach=[(1<<i)|sum(a[i][j]<<j for j in range(n)) for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]&(1<<k):reach[i]|=reach[k]
    return all(x==(1<<n)-1 for x in reach)

def run_additional():
    graph_pairs=0;graph_auto=0
    for n,base in [(i,2) for i in range(5)]+[(3,4)]:
        data=[graph(n,b,base) for b in range(base**comb(n,2))]
        keys=[graph_key(a) for a in data]
        for i,a in enumerate(data):
            brute=sum(transport(a,p)==a for p in permutations(range(n)))
            require(graph_isomorphisms(a,a,True)==brute,'generic graph automorphism brute control');graph_auto+=1
            for j in range(i,len(data)):
                require(bool(graph_isomorphisms(a,data[j]))==(keys[i]==keys[j]),'generic graph class brute control');graph_pairs+=1
    spectra={};count=0
    for n in range(3,7):
        hs=set()
        for word in range(1<<comb(n,2)):
            a=orient(n,word);h=hp(a);require(h not in (7,21),'spectrum hole');count+=1
            if strongly_connected(a):hs.add(h)
        spectra[n]=sorted(hs)
    require(spectra=={3:[3],4:[5],5:[9,11,13,15],6:[15,17,19,23,25,27,29,31,33,37,41,43,45]},'complete strong spectra')
    # n=0/1/2 exact outputs and empty carriers; reject each all-isolate nonempty carrier.
    boundary=[]
    for n in range(3):
        for word in range(1<<comb(n,2)):
            a=orient(n,word)
            for kind in ('joint','gamma'):
                g=make(a,kind);require(strict_decode(g,kind)==a,'small boundary')
                boundary.append({'n':n,'word':word,'kind':kind,'vertices':len(g),'edges_by_color':dict(Counter(g[i][j] for i,j in combinations(range(len(g)),2) if g[i][j])),'automorphisms':graph_isomorphisms(g,g,True)})
    for n in (1,7,9,14,18):
        for kind in ('joint','gamma'):
            try:strict_decode([[0]*n for _ in range(n)],kind)
            except ValueError:pass
            else:raise RuntimeError('nonempty edgeless graph incorrectly recognized')
    print(json.dumps({'generic_iso_brute_controls':{'all_graph_pair_comparisons':graph_pairs,'all_automorphism_counts':graph_auto,'universe':'all simple graphs through n=4; all four-state symmetric relations on n=3'},'independent_all_labelled_spectra':{'tournaments_orders_3_through_6':count,'strong_spectra':spectra},'boundaries':boundary,'all_isolate_rejections':10},indent=2,sort_keys=True))

if __name__=='__main__':run_additional()
