"""Independent local-system audit on actual graphs, without producer imports."""
from collections import Counter
from itertools import permutations,product
from pathlib import Path
import hashlib,json

checks=0
def need(value,why):
    global checks
    checks+=1
    if not value: raise RuntimeError(why)
def parity(x): return x.bit_count()%2
EVEN=(0,3,5,6)
PS=list(permutations(range(3)))
def pmap(p,x): return sum(((x>>j)&1)<<p[j] for j in range(3))
def rank(columns):
    basis={}
    for c in columns:
        while c:
            top=c.bit_length()-1
            if top in basis: c ^= basis[top]
            else:
                basis[top]=c
                break
    return len(basis)
def eencode(x):
    need(x in EVEN,'even module closes')
    return ((x>>1)&1)|(((x>>2)&1)<<1)
def differential_rank(n,edges,linear,even):
    basis=(3,5) if even else (1,2,4)
    width=len(basis)
    columns=[]
    for vertex in range(n):
        for value in basis:
            out=0
            for k,((a,b),p) in enumerate(zip(edges,linear)):
                v=(value if vertex==b else 0) ^ (pmap(p,value) if vertex==a else 0)
                out |= (eencode(v) if even else v)<<(width*k)
            columns.append(out)
    return rank(columns)

# Independently generate the 48 cubical maps by fixing all pair distances.
cube_edges=[(x,x^(1<<j)) for x in range(8) for j in range(3) if x<(x^(1<<j))]
aut=[p for p in permutations(range(8))
     if all((p[a]^p[b]).bit_count()==1 for a,b in cube_edges)]
need(len(aut)==48,'cube group order')
split={}
for p in aut:
    flip=parity(p[0])
    tetra=tuple(EVEN.index(p[e] ^ (7 if parity(p[e]) else 0)) for e in EVEN)
    split[p]=(flip,tetra)
    for x in range(8):
        eta=parity(x); e=x^(7 if eta else 0)
        out_even=EVEN[tetra[EVEN.index(e)]]
        expected=out_even^(7 if eta^flip else 0)
        need(expected==p[x],'all eight states split')
need(len(set(split.values()))==48,'C2 times all S4')
for p in aut:
    for q in aut:
        composition=tuple(p[q[i]] for i in range(8))
        a,A=split[p]; b,B=split[q]
        need(split[composition]==(a^b,tuple(A[B[i]] for i in range(4))),
             'split is a group homomorphism')

graphs={
    'triangle':(3,((0,1),(1,2),(2,0))),
    'three_parallel_edges':(2,((0,1),(0,1),(0,1))),
}
summary={}
for name,(n,edges) in graphs.items():
    e=len(edges); b1=e-n+1
    census=Counter()
    instances=0
    for linear in product(PS,repeat=e):
        rankE=differential_rank(n,edges,linear,True)
        rankA=differential_rank(n,edges,linear,False)
        eps=2*n-rankE
        hidden=2*e-rankE
        need(hidden==2*(b1-1)+eps,'Euler local-system dimension')
        need(3*n-rankA==1+eps,'constant odd splitting of invariants')
        need(3*e-rankA==b1+hidden,'H1 splitting')
        scalar_ok=frame_ok=0
        for gains in product(range(8),repeat=e):
            instances+=1
            if name=='triangle':
                scalar=(parity(gains[0]^gains[1]^gains[2])==0)
            else:
                scalar=(parity(gains[0])==parity(gains[1])==parity(gains[2]))
            scalar_ok+=scalar
            sections=0
            for root in range(8):
                if name=='triangle':
                    x1=pmap(linear[0],root)^gains[0]
                    x2=pmap(linear[1],x1)^gains[1]
                    valid=(pmap(linear[2],x2)^gains[2])==root
                else:
                    x1=pmap(linear[0],root)^gains[0]
                    valid=all((pmap(linear[k],root)^gains[k])==x1 for k in (1,2))
                sections+=valid
            if sections:
                frame_ok+=1
                need(scalar and sections==2**(1+eps),'all sections have invariant-space size')
            census[(eps,scalar,bool(sections))]+=1
        need(frame_ok==2**rankA,'full compatible gain data are im d0')
        need(scalar_ok==2**(3*e-b1),'all scalar-compatible data counted')
        need(scalar_ok==frame_ok*2**hidden,'hidden cohomology fibres')
    summary[name]={'vertices':n,'edges':e,'cycle_rank':b1,'instances':instances,
                   'census':[{'epsilon':eps,'completion':scalar,'full_frame':full,'count':count}
                             for (eps,scalar,full),count in sorted(census.items())]}

# The distinction is chosen-loop fixed points versus a common fixed point.
cycle=lambda x: ((x<<1)&7)|(x>>2)
h0=tuple(cycle(x) for x in range(8))
h1=tuple(cycle(x)^5 for x in range(8))
fixed0=[x for x in range(8) if h0[x]==x]
fixed1=[x for x in range(8) if h1[x]==x]
need(fixed0==[0,7] and fixed1==[3,4],'exact hostile fixed loci')
need(set(fixed0).isdisjoint(fixed1),'no common fixed point')
inverse0=tuple(h0.index(x) for x in range(8))
word=tuple(h1[inverse0[x]] for x in range(8))
need(word==tuple(x^5 for x in range(8)),'generated holonomy has translation5')
need(all(word[x]!=x for x in range(8)),'the generated word is fixed-point-free')

print(json.dumps({'status':'PASS independent finite checks', 'checks':checks,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'graph_universes':summary,'group_split':'48=2*24; all 2304 products checked',
 'hostile':{'fixed_sets':[fixed0,fixed1],'fixed_free_word':list(word)}},
 indent=2,sort_keys=True))
