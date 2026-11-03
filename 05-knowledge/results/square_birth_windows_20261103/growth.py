"""Exact prefix growth and typed prime comparisons; standard library only."""
from pathlib import Path
from math import isqrt
from collections import deque
from itertools import combinations
import json
from birth import bank, need

HERE=Path(__file__).resolve().parent

def square_edges(vertices):
    V=set(vertices)
    if not V:return []
    return sorted((a,s*s-a) for s in range(2,isqrt(2*max(V))+1)
                  for a in sorted(V) if a<s*s-a and s*s-a in V)

def components(V,E):
    adj={v:set() for v in V}
    for a,b in E:adj[a].add(b);adj[b].add(a)
    left=set(V);out=[]
    while left:
        seen={min(left)};todo=list(seen)
        while todo:
            for v in adj[todo.pop()]-seen:seen.add(v);todo.append(v)
        out.append(sorted(seen));left-=seen
    return sorted(out,key=lambda c:(-len(c),c))

def paths(n,E):
    if n==1:return [(1,)]
    if len(components(range(1,n+1),E))!=1:return []
    adj={v:[] for v in range(1,n+1)}
    for a,b in E:adj[a].append(b);adj[b].append(a)
    leaves=[v for v in adj if len(adj[v])==1]
    if len(leaves)>2:return []
    out=[]
    def walk(p,seen):
        if len(p)==n:out.append(p);return
        for v in adj[p[-1]]:
            if v not in seen:walk(p+(v,),seen|{v})
    for v in leaves[:1] if leaves else adj:walk((v,),{v})
    if leaves:out += [p[::-1] for p in out]
    return sorted(out)

def dag(E,x,n):
    adj=[[] for _ in range(n+1)];indeg=[0]*(n+1)
    for k,(a,b) in enumerate(E):
        if not(x>>k&1):a,b=b,a
        adj[a].append(b);indeg[b]+=1
    q=deque(v for v in range(1,n+1) if indeg[v]==0);count=0
    while q:
        v=q.popleft();count+=1
        for w in adj[v]:
            indeg[w]-=1
            if indeg[w]==0:q.append(w)
    return count==n

def main():
    ledger=[]
    for n in range(1,301):
        E=square_edges(range(1,n+1))
        # Independent direct pair enumeration and the old birth formula.
        direct=[e for e in combinations(range(1,n+1),2) if isqrt(sum(e))**2==sum(e)]
        need(E==direct,'independent square-edge construction')
        K=isqrt(n);M=isqrt(2*n-1)
        exact12=12*(M-K)*n-12*K+2*K*(K+1)*(2*K+1)-M*(M+1)*(2*M+1)+6*((M+1)//2)
        need(exact12==12*len(E),'incoming Session148 cumulative edge formula')
        roots=list(range(isqrt(n)+1,isqrt(2*n-1)+1))
        born=sorted(s*s-n for s in roots)
        need(born==sorted(a for a,b in E if b==n),'literal birth ports')
        C=components(range(1,n+1),E)
        row={'n':n,'roots':roots,'parents':born,'birth':len(roots),
             'edges':len(E),'component_sizes':[len(c) for c in C],
             'cycle_rank':len(E)-n+len(C),'excess_over_twice_path':len(E)-2*(n-1)}
        if n<=28:
            P=paths(n,E);row['directed_paths']=len(P)
            if P:
                used={tuple(sorted(e)) for p in P for e in zip(p,p[1:])}
                row['unused_in_every_path']=sorted(set(E)-used)
        ledger.append(row)
    need([r['n'] for r in ledger[:25] if r['directed_paths']]==[1,15,16,17,23,25],
         'inherited small Hamilton spectrum')
    threshold=next(r['n'] for r in ledger[1:] if r['excess_over_twice_path']>0)
    need(threshold==56,'first unused-edge majority')
    need(all(r['excess_over_twice_path']>0 for r in ledger[55:]),'finite tail check')
    witness56=json.loads((HERE/'path56.json').read_text(encoding='utf-8'))
    need(sorted(witness56)==list(range(1,57)),'nonvacuous threshold witness vertices')
    need(all(isqrt(a+b)**2==a+b for a,b in zip(witness56,witness56[1:])),
         'nonvacuous threshold witness edges')
    prime=[n for n in range(2,100) if all(n%d for d in range(2,isqrt(n)+1))]
    ogg=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
    prime_ledger=[]
    for k,p in enumerate(prime,1):
        V=prime[:k];E=square_edges(V);C=components(V,E)
        A=[v for v in V if v==2 or v%4==1];B=[v for v in V if v%4==3]
        need(all((a in A)!=(b in A) for a,b in E),'prime bipartition')
        prime_ledger.append({'index':k,'prime':p,'prime_graph_edges':len(E),
            'prime_graph_components':C,'bipartition_sizes':[len(A),len(B)],
            'integer_Q_index_edges':ledger[k-1]['edges'],
            'integer_Q_prime_edges':ledger[p-1]['edges'],
            'is_ogg_prime':p in ogg})
    need([len(c) for c in prime_ledger[-1]['prime_graph_components']]==[22,1,1,1],
         'inherited prime25 disconnected control')
    need(prime_ledger[-1]['prime_graph_edges']==24,'prime25 edge control')
    ogg_A=[v for v in ogg if v==2 or v%4==1];ogg_B=[v for v in ogg if v%4==3]
    need([len(ogg_A),len(ogg_B)]==[6,9],'Ogg square-sum bipartition obstruction')
    S=list(range(1,13))+[16,17,18];SE=square_edges(S);SC=components(S,SE)
    need([len(c) for c in SC]==[6,5,4],'user component control')
    orientation=[]
    for n in (25,26,27):
        g=bank(n);P=paths(n,square_edges(range(1,n+1)))
        need(P==g['P'],'independent path storage and traversal')
        dag_count=sum(dag(g['E'],x,n) for x in g['active'])
        need(dag_count==len(P),'unique-order orientation bijection')
        pos={v:i for i,v in enumerate(P[0])}
        forward=sum(1<<k for k,(a,b) in enumerate(g['E']) if pos[a]<pos[b])
        complemented=forward^(g['full']^g['patterns'][0][0])
        need(dag(g['E'],forward,n) and not dag(g['E'],complemented,n),
             'path-order orientations fail actual complement closure')
        orientation.append({'n':n,'active_guard_orientations':len(g['active']),
                            'acyclic_active_guard_orientations':dag_count})
    result={'status':'FINITE-EXACT; all-size formulas are proved in GROWTH.md',
        'prefix_ledger':ledger,'first25_prime_comparisons':prime_ledger,
        'first_unused_majority':threshold,'ogg_A':ogg_A,'ogg_B':ogg_B,
        'user_sets_components':SC,'user_sets_edges':SE,'orientation_controls':orientation}
    (HERE/'growth.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    table=['# First 25 births and three different prime-indexed graphs','',
           'Generated by `growth.py`. H counts directed paths; the singleton counts once.',
           'Prime edges/components refer to the graph on the actual first k prime VALUES.',
           'Q_k instead has integer vertices 1 through k. Q_(p_k) has vertices 1 through p_k.','',
           '| k | p_k | New parents in Q_k | b(k) | E(Q_k) | H(Q_k) | E(prime prefix) | Prime components | E(Q_(p_k)) |',
           '|---:|---:|:---|---:|---:|---:|---:|---:|---:|']
    for r,s in zip(ledger,prime_ledger):
        values=[r['n'],s['prime'],','.join(map(str,r['parents'])) or '-',r['birth'],
                r['edges'],r['directed_paths'],s['prime_graph_edges'],
                len(s['prime_graph_components']),s['integer_Q_prime_edges']]
        table.append('| '+' | '.join(map(str,values))+' |')
    (HERE/'LEDGER.md').write_text('\n'.join(table)+'\n',encoding='utf-8')
    print(json.dumps({'first_unused_majority':threshold,
        'windows':[ledger[n-1] for n in (8,9,10,15,24,25,26,27,28,55,56)],
        'prime25':prime_ledger[-1],'ogg_parts':[ogg_A,ogg_B],
        'acyclic':orientation},sort_keys=True))

if __name__=='__main__':main()
