"""Independent finite controls for the proved all-order growth theorem."""
from collections import Counter
from fractions import Fraction
from itertools import combinations,permutations,product
from math import factorial
from pathlib import Path
from functools import lru_cache
import argparse,json,importlib.util
from port_scout import add,converse,encode,rows
from growth_identity import sink,edge_set
HERE=Path(__file__).resolve().parent
checks=0
def need(x,why):
    global checks
    checks+=1
    if not x:raise RuntimeError(why)

@lru_cache(maxsize=512)
def count(a):
    n=len(a);dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in dp[mask].items():
            opts=a[v]&~mask
            while opts:
                w=(opts&-opts).bit_length()-1;opts&=opts-1
                d=dp[mask|(1<<w)];d[w]=d.get(w,0)+c
    return sum(dp[-1].values())

def cylinders(n):
    edges=list(combinations(range(n),2));idx={p:i for i,p in enumerate(edges)};full=(1<<len(edges))-1;out=[0]*(full+1)
    for word in permutations(range(n)):
        on=used=0
        for u,v in zip(word,word[1:]):
            bit=1<<idx[min(u,v),max(u,v)];used|=bit
            if u<v:on|=bit
        free=full^used;s=free
        while True:
            out[on|s]+=1
            if not s:break
            s=(s-1)&free
    return out

def growth_terms(a,x,y):
    b=sink(a)
    return [(1,add(b,y)),(-1,add(b,x)),(-1,add(add(b,x),y)),(1,add(add(b,y),x))]

def growth_jet(a,es):
    terms=[(1,a)]
    for x,y in es:terms=[(sg*t,b) for sg,v in terms for t,b in growth_terms(v,x,y)]
    return sum(sg*count(tuple(a)) for sg,a in terms)

def direct_jet(a,es):
    if not es:return count(tuple(a))
    x,y=es[0]
    return direct_jet(edge_set(a,x,y,1),es[1:])-direct_jet(edge_set(a,x,y,0),es[1:])

def hessian(n):
    edges=list(combinations(range(n),2));idx={p:i for i,p in enumerate(edges)};m=len(edges);raw=[[0]*m for _ in range(m)]
    for word in permutations(range(n)):
        signed=[(idx[min(u,v),max(u,v)],1 if u<v else -1) for u,v in zip(word,word[1:])]
        for (e,sg),(f,tg) in combinations(signed,2):raw[e][f]+=sg*tg;raw[f][e]+=sg*tg
    D=[[int(v==b)-int(v==a) for a,b in edges] for v in range(n)]
    for e in range(m):
        for f in range(m):
            expected=2*factorial(n-2)*(2*int(e==f)-sum(row[e]*row[f] for row in D))
            need(raw[e][f]==expected,'literal permutation Hessian')
    return {'n':n,'c':str(Fraction(2*factorial(n-2),2**(n-3))),'cycle_multiplicity':(n-1)*(n-2)//2,'cut_multiplicity':n-1}

def hostile():
    out=[]
    for flip in (False,True):
        a=[0]*12
        for u,v in combinations(range(12),2):
            b,i=divmod(u,4);c,j=divmod(v,4)
            arc=(c-b)%3==1 if b!=c else i==0 or (i>0 and (j-i)%3==1)
            if flip and b==c==2:arc=not arc
            a[u if arc else v]|=1<<(v if arc else u)
        out.append(a)
    cyclic=lambda a:{t for t in combinations(range(12),3) if all(sum(a[v]>>w&1 for w in t if w!=v)==1 for v in t)}
    need(cyclic(out[0])==cyclic(out[1]),'actual hostile cyclic triples')
    f=HERE.parent/'regular_parent_symmetry_20261029/verify.py'
    spec=importlib.util.spec_from_file_location('retained_hostile',f);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    orders=[len(mod.iso_maps(a,a)) for a in out];need(orders==[81,27],'actual automorphism groups')
    values=[count(tuple(a)) for a in out]
    probes=[count(tuple(add(sink(a),8))) for a in out]
    need(values[0]==values[1],'same raw H hostile')
    need(probes[0]!=probes[1],'named birth separates hostile')
    return {'groups':orders,'common_H':values,'port8_after_sink':probes,'common_cyclic_triples':len(cyclic(out[0]))}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='verification.json');args=parser.parse_args()
    ports=json.loads((HERE/'port_scout.json').read_text());jets=[]
    for n in range(2,7):
        vals=cylinders(n);m=n*(n-1)//2;seen={}
        for mask,H in enumerate(vals):
            sig=(H,)+tuple(vals[mask|(1<<e)]-vals[mask&~(1<<e)] for e in range(m))
            seen.setdefault(sig,[]).append(mask)
            if n in (3,4):
                edges=list(combinations(range(n),2))
                incidence=[[int(v==b)-int(v==a) for a,b in edges] for v in range(n)]
                g=sig[1:]
                dg=[sum(t*x for t,x in zip(row,g)) for row in incidence]
                c=Fraction(2*factorial(n-2),2**(n-3))
                recovered=[Fraction(1,2)+(g[e]-Fraction(sum(incidence[v][e]*dg[v] for v in range(n)),n-2))/(2*c) for e in range(m)]
                need(recovered==[(mask>>e)&1 for e in range(m)],'explicit small-order Boolean decoder')
        for entry in ports[n-2]['entries']:
            a=rows(n,entry['mask']);need(vals[entry['mask']]==entry['H'],'path-cylinder class count')
            need([count(tuple(add(a,x))) for x in range(n)]==entry['children'],'literal point births')
        jets.append({'n':n,'states':len(seen),'orientations':len(vals),'max_fibre':max(map(len,seen.values())),'collision':next((v for v in seen.values() if len(v)>1),None)})
    hs=[hessian(n) for n in range(3,8)]
    growth=[]
    cases=[(3,0,[(0,1)]),(3,2,[(0,1),(1,2)]),(3,5,[(0,2),(0,1)]),(4,4,[(0,1),(2,3)]),(4,17,[(1,3),(0,2)])]
    for n,m,es in cases:
        a=rows(n,m);u=growth_jet(a,es);v=direct_jet(a,es);need(u==v,'iterated growth derivative')
        growth.append({'n':n,'mask':m,'edges':es,'value':u,'max_final_order':n+3*len(es)})
    tree=json.loads((HERE/'tree_scan11.json').read_text(encoding='utf-8-sig'))
    old=json.loads((HERE/'tree_scan.json').read_text(encoding='utf-8-sig'))
    need(tree['layers'][:10]==old['layers'],'two full scans agree through ten')
    need(all(z['presentations']==factorial(z['n']-1) for z in tree['layers']),'entire recursive-tree presentation universe')
    need(all(str(x) in tree['witnesses'] for x in range(1,20407,2) if x not in (7,21)),'literal witness prefix')
    for i,(target,word) in enumerate(tree['witnesses'].items()):
        if i%97==0 or int(target) in (3,5,11,13,17,19,35,49,63,343,20405,25827):
            a=[0]
            for x in word:a=add(a,x)
            need(count(tuple(a))==int(target),'independent recursive-tree witness')
    source=hostile()
    parity_contexts=0
    for k in range(4):
        for word in product(range(4),repeat=k):
            values=[]
            for mask in range(8):
                a=rows(3,mask)
                for op in word:a=sink(a) if op==3 else add(a,op)
                values.append(count(tuple(a)))
            need(sum(v*(-1)**m.bit_count() for m,v in enumerate(values))==0,'parity laws have equal context means')
            parity_contexts+=1
            if word==(0,1):parity_values=values
    moments=[[Fraction(sum(v**r for m,v in enumerate(parity_values) if m.bit_count()%2==p),4) for r in (1,2)] for p in (0,1)]
    need(moments==[[10,113],[10,111]],'sharp second-moment control')
    for mask in range(8):
        a=add(add(rows(3,mask),0),1);k=len(a)
        joined=[r|(((1<<k)-1)<<k) for r in a]+[r<<k for r in a]
        need(count(tuple(joined))==parity_values[mask]**2,'actual shared-copy ordinal product')
    first7=json.loads((HERE/'first_jet7.json').read_text(encoding='utf-8-sig'))
    need(first7['orientations']==1<<21 and first7['first_jet_states']==1<<21 and first7['exact_collision_comparisons'],'complete first jet at seven receipt')
    all_tree=json.loads((HERE/'tree_independent.json').read_text(encoding='utf-8-sig'))
    need(all_tree['verified_witnesses']==len(tree['witnesses'])==11888,'independent all-witness replay')
    result={'checks':checks,'first_jet':jets,'first_jet7':first7,'literal_hessians':hs,'iterated_growth':growth,'hostile':source,'tree_prefix_through':20405,'tree_classes_are_presentations':True,'parity_contexts':parity_contexts,'parity_example':{'values':parity_values,'moments':[[str(v) for v in p] for p in moments],'all_eight_shared_copies_replayed':True}}
    (HERE/args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
