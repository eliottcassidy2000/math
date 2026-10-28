"""Bit squares, aligned consecutive tournaments, and odd self-converse doubles."""
from itertools import combinations,permutations
from pathlib import Path
import json
import reconstruct as r

HERE=Path(__file__).resolve().parent
refine=r.old.source('pair_refine','05-knowledge/results/root_ports_closed_carriers_20261028/check_carriers.py')
channels=r.old.source('pair_channels','05-knowledge/results/converse_tournament_count_channels_20261021/session22-converse-channels/channels.py')

def relation(n,fn):
    rows=[0]*n
    for i,j in combinations(range(n),2):
        if fn(i,j):rows[i]|=1<<j
        else:rows[j]|=1<<i
    return tuple(rows)

def decode(n,bits):
    bit=lambda i,j:(bits>>(i*n+j))&1
    a=relation(n,lambda i,j:bit(i,j))
    b=relation(n+1,lambda i,j:(1-bit(i,i)) if j==n else bit(j,i))
    return a,b

def encode(a,b):
    n=len(a);bits=0
    for i in range(n):
        for j in range(n):
            bit=r.core.edge(a,i,j) if i<j else r.core.edge(b,j,i) if i>j else r.core.edge(b,n,i)
            bits|=bit<<(i*n+j)
    return bits

def double(a,b):
    n=len(a);pole=2*n
    def arc(i,j):
        if j==pole:return i>=n
        if j<n:return r.core.edge(a,i,j)
        if i>=n:return 1-r.core.edge(a,i-n,j-n)
        k=j-n
        return r.core.edge(b,n,i) if i==k else r.core.edge(a,i,k)^r.core.edge(b,i,k)
    return relation(2*n+1,arc)

def undouble(t):
    n=(len(t)-1)//2
    a=tuple(row&((1<<n)-1) for row in t[:n])
    b=relation(n+1,lambda i,j:1-r.core.edge(t,i,n+i) if j==n else r.core.edge(t,i,n+j)^r.core.edge(a,i,j))
    return a,b

def main():
    out=[]
    for n in range(4):
        ps=list(permutations(range(n)));pairs=set();naive=set();doubles=set()
        engine=refine.Canon()
        for bits in range(1<<(n*n)):
            a,b=decode(n,bits);t=double(a,b);rho=tuple(range(n,2*n))+tuple(range(n))+(2*n,)
            r.need(encode(a,b)==bits and undouble(t)==(a,b),'two literal inverses')
            r.need(r.core.relabel(r.core.converse(t),rho)==t,'actual anti-involution')
            orbit=[];matrix_orbit=[]
            for p in ps:
                aa=r.core.relabel(a,p);bb=r.core.relabel(b,p+(n,))
                orbit.append(encode(aa,bb))
                pp=p+tuple(n+i for i in p)+(2*n,)
                r.need(double(aa,bb)==r.core.relabel(t,pp),'doubled equivariance')
                matrix_orbit.append(sum(((bits>>(p[i]*n+p[j]))&1)<<(i*n+j) for i in range(n) for j in range(n)))
            pairs.add(min(orbit));naive.add(min(matrix_orbit));doubles.add(engine.maps(t)[0])
        q,blue=r.old.native_classes(2*n+1);reps={}
        for x,c in enumerate(q):reps.setdefault(c,x)
        sc=set()
        for x in reps.values():
            t=r.core.presentation(2*n+1,x)
            if r.oracle.isomorphisms(t,r.core.converse(t)):
                sc.add(engine.maps(t)[0])
        r.need(doubles==sc and len(pairs)==len(sc),'full independent odd self-converse census')
        out.append(dict(n=n,bit_squares=1<<(n*n),synchronized_pair_classes=len(pairs),
                        naive_matrix_conjugacy_classes=len(naive),self_converse_classes=len(sc)))
    r.need(out[2]['synchronized_pair_classes']==8 and out[2]['naive_matrix_conjugacy_classes']==10,'non-equivariance hostile')
    values,stats=channels.channel_dp(100,a=2,b=1)
    r.need(all(x.denominator==1 for x in values),'paired integer counts')
    values=list(map(int,values))
    r.need(values[:4]==[x['synchronized_pair_classes'] for x in out],'independent channel check')
    paley=tuple(sum(1<<((i+d)%7) for d in (1,2,4)) for i in range(7))
    group=r.oracle.isomorphisms(paley,paley,True);rho=tuple((-i)%7 for i in range(7))
    central=[p for p in group if all(p[rho[i]]==rho[p[i]] for i in range(7))]
    phase=tuple((2*i)%7 for i in range(7))
    r.need(len(group)==21 and len(central)==3 and phase in central,'Paley pole centralizer retains Frobenius phase')
    result=dict(status='FINITE-EXACT',exhaustive=out,paired_counts_through_100=values,
                inherited_channel_stats=stats,Paley=dict(automorphisms=21,anti_involution_centralizer=central,
                anti_involution=rho,Frobenius_phase=phase),
                distinction='Class bijection is classical Eplett plus the inherited XOR map; it need not preserve full automorphism groups.')
    (HERE/'pair-square-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(exhaustive=out,first_paired_counts=values[:8],Paley_centralizer=3)))

if __name__=='__main__':main()
