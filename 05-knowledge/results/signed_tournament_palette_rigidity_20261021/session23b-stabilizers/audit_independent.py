#!/usr/bin/env python3
"""Second representation: vector-linear Fano action and set partitions.
No producer imports. Exhaustive finite controls do not establish general theorems.
"""
from itertools import permutations,product,combinations
from math import prod, comb
from collections import Counter
from pathlib import Path
import argparse,json,hashlib

def need(x,msg):
    if not x:raise RuntimeError(msg)
def moved(mask,p):return sum(1<<p[i] for i in range(len(p)) if mask>>i&1)
def regs(G):
    n=len(G[0]); return [s for s in range(1<<n) if sum(moved(s,p)==s for p in G)==1]
def partitions(n,k):
    def rec(word,m):
        if len(word)==n:
            yield tuple(word);return
        for a in range(min(m+1,k-1)+1):yield from rec(word+[a],max(m,a))
    yield from rec([0],0)
def regular_partitions(G,k):
    n=len(G[0]);out=Counter();first={}
    for c in partitions(n,k):
        blocks=[sum(1<<i for i in range(n) if c[i]==v) for v in range(max(c)+1)]
        if sum(all(moved(s,p)==s for s in blocks) for p in G)==1:
            b=len(blocks);out[b]+=1;first.setdefault(b,c)
    return out,first
def gf8_mul(a,b):
    z=0
    while b:
        if b&1:z^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return z

def main():
    affine={a:[tuple((a*x+b)%7 for x in range(7)) for b in range(7)] for a in range(1,7)}
    G21=sum((affine[a] for a in (1,2,4)),[]);G42=sum(affine.values(),[])
    C=[1,2,4,3,6,7,5];Ci={x:i for i,x in enumerate(C)}
    GL=[]
    for a,b,c in permutations(range(1,8),3):
        cols=(a,b,c)
        table=[0]+[0 for _ in range(7)]
        for x in range(1,8):
            for bit in range(3):
                if x>>bit&1:table[x]^=cols[bit]
        if len(set(table))==8:GL.append(tuple(Ci[table[x]] for x in C))
    need(len(GL)==168,'GL order')
    source=json.loads(Path(__file__).with_name('receipt.json').read_text())
    result={}
    for name,G in [('paley_aut',G21),('paley_aut_pm',G42),('fano_aut',GL)]:
        counts,first=regular_partitions(G,4)
        qc={q:sum(v*prod(range(q-b+1,q+1)) for b,v in counts.items() if b<=q) for q in range(1,5)}
        for q,record in source['groups'][name]['distinguishing_colors'].items():
            need(qc[int(q)]==record['regular_colorings'],name+' color count')
        rr=regs(G);need(len(rr)==source['groups'][name]['regular_subsets'],name+' subset count')
        result[name]={'group_order':len(G),'regular_set_partitions_by_blocks':dict(counts),'regular_ordered_colorings':qc}
    G8=[]
    for a in range(1,8):
        for b in range(8):
            for j in range(3):
                p=[]
                for x in range(8):
                    y=x
                    for _ in range(j):y=gf8_mul(y,y)
                    p.append(gf8_mul(a,y)^b)
                G8.append(tuple(p))
    need(len(set(G8))==168,'affine semilinear order')
    need(not regs(G8),'degree8 exceptional regular subsets absent')
    result['degree8_affine_semilinear']={'degree':8,'order':168,'regular_subsets':0,'max_constant_cardinality_layer':comb(8,4)}
    product_records=[]
    for H,t,u in [(G21,1,1),([tuple(range(0))],3,2),([tuple((i+a)%3 for i in range(3)) for a in range(3)],2,1)]:
        b=len(H[0]);n=b+3*t+u
        G=[]
        for h in H:
            for rotations in product(range(3),repeat=t):
                p=list(h)
                for j,r in enumerate(rotations):p.extend(b+3*j+(i+r)%3 for i in range(3))
                p.extend(range(b+3*t,n));G.append(tuple(p))
        rH=regs(H);rG=regs(G)
        rhoH=min(s.bit_count() for s in rH);rhoG=min(s.bit_count() for s in rG)
        need(len(rG)==len(rH)*6**t*2**u,'product regular count')
        need(rhoG==rhoH+t,'product minimum cost')
        need(len(rG)%len(G)==0,'regular orbit count')
        poly=Counter(x.bit_count() for x in rH)
        for _ in range(t):
            nxt=Counter()
            for k,v in poly.items():nxt[k+1]+=3*v;nxt[k+2]+=3*v
            poly=nxt
        for _ in range(u):
            nxt=Counter()
            for k,v in poly.items():nxt[k]+=v;nxt[k+1]+=v
            poly=nxt
        need(poly==Counter(x.bit_count() for x in rG),'regular subset polynomial')
        product_records.append({'B':b,'H_order':len(H),'t':t,'u':u,'R_H':len(rH),'R_T':len(rG),'rho_H':rhoH,'rho_T':rhoG,'regular_orbits':len(rG)//len(G)})
    result['disjoint_product_controls']=product_records
    # Every labelled tournament n<=4, every ordered directed nonempty partial path.
    cells=0
    for n in range(1,5):
        E=list(combinations(range(n),2));P=list(permutations(range(n)))
        for tm in range(1<<len(E)):
            arcs=set((i,j) if tm>>k&1 else (j,i) for k,(i,j) in enumerate(E))
            G=[p for p in P if all((p[i],p[j]) in arcs for i,j in arcs)]
            for k in range(1,n+1):
                for path in permutations(range(n),k):
                    if any((path[i],path[i+1]) not in arcs for i in range(k-1)):continue
                    U=[v for v in range(n) if v not in path];ui={v:i for i,v in enumerate(U)}
                    H=[p for p in G if all(p[v]==v for v in path)]
                    HU=[tuple(ui[p[v]] for v in U) for p in H]
                    need(len(set(HU))==len(H),'faithful residual action')
                    need(len(H)%2==1,'odd path stabilizer')
                    rr=regs(HU);need(bool(rr),'residual regular subset exists')
                    need(min(s.bit_count() for s in rr)<=len(U)//2,'residual minority bound')
                    cells+=1
    result['partial_path_control']={'all_labelled_tournaments_n_le_4':True,'ordered_path_objects_checked':cells}
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args=argparse.ArgumentParser();args.add_argument('--output',default='audit-receipt.json');a=args.parse_args()
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
