"""Independent all-k DAG control and |Q|=2 bridge-sink boundary test."""
from check_independent import covers,need
from itertools import permutations,combinations
from pathlib import Path
import json

def topo(D):
    n=len(D);left=set(range(n));order=[]
    while left:
        v=next((v for v in left if not any(D[u][v] for u in left)),None)
        if v is None:return None
        order.append(v);left.remove(v)
    return {v:i for i,v in enumerate(order)}

def dag_tests():
    total=0;levels=0
    for n in (3,4):
        extra=[(u,v) for u in range(n) for v in range(n) if v!=u and v!=(u+1)%n]
        for bits in range(1<<len(extra)):
            A=[[int(v==(u+1)%n) for v in range(n)] for u in range(n)]
            for b,(u,v) in enumerate(extra):A[u][v]=(bits>>b)&1
            D=[[0]*n for _ in range(n)]
            for u,v in extra:
                if A[u][v]:D[(u+1)%n][v]=1
            f=topo(D)
            if f is None:continue
            total+=1
            for k in range(1,n+1):
                ports=list(permutations(range(n),k));match={}
                for s in ports:
                    # Stop just before the next selected start on the cycle.
                    end=[]
                    for a in s:
                        v=(a+1)%n
                        while v not in s:v=(v+1)%n
                        end.append((v-1)%n)
                    match[s]=tuple(end)
                need(len(set(match.values()))==len(ports),('bijection',n,k))
                inv={t:s for s,t in match.items()};M=covers(A,k)
                for s in ports:need(M[s,match[s]]==1,('unit',n,k,s))
                for (s,t),v in M.items():
                    if t!=match[s]:need(sum(f[u] for u in s)<sum(f[u] for u in inv[t]),('potential',n,k,s,t))
                levels+=1
    return {'digraphs':total,'all_k_observers':levels,'orders':[3,4]}

def sink_tests():
    total=0
    for n in (3,4):
        a,s=0,1;other=list(combinations(range(1,n),2))
        for bits in range(1<<len(other)):
            R=[[0]*n for _ in range(n)];R[s][a]=1
            for v in range(2,n):R[a][v]=1
            for b,(u,v) in enumerate(other):R[u][v]=(bits>>b)&1;R[v][u]=1-R[u][v]
            for q in (a,s):
                r,c=n,n+1;A=[row+[0,0] for row in R]+[[0]*(n+2) for _ in range(2)]
                for u in range(n):A[r][u]=A[c][u]=1
                A[q][r]=1;A[r][q]=0;A[c][r]=1
                M=covers(A,2)
                for st in permutations(range(n+2),2):need(M[st,(a,r)]-M[st,(r,a)]-M[st,(s,r)]+M[st,(r,s)]==0,('sinkQ2',n,bits,q,st))
                total+=1
    return {'Q_size':2,'R_sizes':[3,4],'cases':total}

if __name__=='__main__':
    out={'all_k_DAG':dag_tests(),'bridge_sink_boundary':sink_tests()}
    Path(__file__).with_name('all_k_dag_sink_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
