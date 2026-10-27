#!/usr/bin/env python3
"""Independent literal adjacency-matrix audit; producer programs are not imported."""
from pathlib import Path
from collections import Counter
import hashlib,json
ROOT=Path(__file__).resolve().parent
N=11
QR={1,3,4,5,9}

def require(p,msg):
    if not p: raise RuntimeError(msg)

def score_L(T):
    d=[sum(row) for row in T];s=sorted(d)
    cost=sum(abs(d-i) for i,d in enumerate(s))
    require(cost%2==0,'nonintegral score displacement')
    return cost//2,d

def apply(T,path):
    require(sorted(path)==list(range(N)),'path is not a permutation')
    require(all(T[u][v]==1 for u,v in zip(path,path[1:])),('illegal path',path))
    keep={frozenset((u,v)) for u,v in zip(path,path[1:])}
    U=[[0]*N for _ in range(N)]
    for u in range(N):
        for v in range(N):
            if u!=v:U[u][v]=T[u][v] if frozenset((u,v)) in keep else 1-T[u][v]
    require(all(U[u][v]+U[v][u]==1 for u in range(N) for v in range(u)), 'invalid resulting tournament')
    return U

def row_bits(T):return [sum(T[u][v]<<v for v in range(N)) for u in range(N)]

def parse(name):
    lines=(ROOT/name).read_text().splitlines();out=[]
    for i,line in enumerate(lines):
        if line.startswith('path '):
            path=tuple(map(int,line.split()[1:]));rows_line=lines[i+1] if lines[i+1].startswith('rows ') else lines[i+2]
            require(rows_line.startswith('rows '),('missing rows',name,i))
            out.append((path,list(map(int,rows_line.split()[1:]))))
    return out

initial=[[int(u!=v and (v-u)%N in QR) for v in range(N)] for u in range(N)]
require(score_L(initial)==(15,[5]*11),'P11 initial scores')
c1=parse('C1-six-moves.pinned.out');greedy=parse('extend-greedy.pinned.out');tail=parse('backtrack.pinned.out')
require((len(c1),len(greedy),len(tail))==(6,8,3),'input path counts')
word=c1+greedy[:6]+tail
T=initial;stages=[];states=[initial]
for k,(path,claimed_rows) in enumerate(word,1):
    T=apply(T,path);states.append(T)
    require(row_bits(T)==claimed_rows,('adjacency mismatch',k))
    L,d=score_L(T);require(L==15-k,('score descent failure',k,L))
    stages.append(dict(depth=k,path=path,adjacency_rows=row_bits(T),scores=d,sorted_scores=sorted(d),L=L,distance_from_P11_exact=k,distance_to_transitive_class_exact=15-k))
final_order=(6,0,2,9,7,4,8,5,10,3,1)
require(all(T[u][v]==1 for i,u in enumerate(final_order) for v in final_order[i+1:]),'final order is not transitive')
require(row_bits(T)==[1982,0,1978,2,1322,1034,1983,1338,1066,1466,10],'final row pin')
# Frozen C1 endpoint remains a forest frame; verify again with matrix operations.
forest_order=(1,10,4,6,8,7,5,2,3,9,0);F=states[6]
backward=[(u,v) for i,u in enumerate(forest_order) for v in forest_order[i+1:] if F[v][u]]
adj=[set() for _ in range(N)]
for u,v in backward:adj[u].add(v);adj[v].add(u)
component={0}
while True:
    new=component.union(*(adj[v] for v in component))
    if new==component:break
    component=new
require(len(backward)==10 and len(component)==11,'C1 endpoint lost its tree frame')
# Verify the discarded greedily chosen branch from the common L3 state.
G=states[12];trap_branch=[]
for k,(path,rows) in enumerate(greedy[6:],13):
    G=apply(G,path);L,d=score_L(G)
    require(row_bits(G)==rows and L==15-k,('discarded greedy branch mismatch',k))
    trap_branch.append(dict(depth=k,path=path,adjacency_rows=rows,L=L,scores=d))
require(score_L(G)[0]==1,'greedy trap not at L1')
# All Hamiltonian paths of the L1 trap, not a state-space BFS.
paths=[]
def extend(path,remaining):
    if not remaining:
        paths.append(tuple(path));return
    for v in sorted(remaining):
        if not path or G[path[-1]][v]:extend(path+[v],remaining-{v})
extend([],set(range(N)))
hist=Counter(score_L(apply(G,p))[0] for p in paths)
require(len(paths)==243 and hist.get(0,0)==0,('greedy trap census',len(paths),hist))
# Second way to count paths: subset DP, independent of the DFS.
D=[[0]*N for _ in range(1<<N)]
for v in range(N):D[1<<v][v]=1
for S in range(1,1<<N):
    for v in range(N):
        if S>>v&1:
            R=S^(1<<v)
            if R:D[S][v]=sum(D[R][u] for u in range(N) if R>>u&1 and G[u][v])
require(sum(D[-1])==len(paths),'independent path-count DP mismatch')
inputs={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('C1-six-moves.pinned.out','extend-greedy.pinned.out','backtrack.pinned.out')}
report=dict(status='PASS',input_hashes=inputs,initial_scores=[5]*11,initial_L=15,stages=stages,final_transitive_forward_order=final_order,final_adjacency_rows=row_bits(T),exact_distance_to_transitive=15,exact_distance_to_forest=6,discarded_greedy_branch=trap_branch,greedy_L1_trap=dict(hamiltonian_paths=len(paths),independent_subset_DP_path_count=sum(D[-1]),successor_L_histogram=dict(sorted(hist.items())),transitive_one_step_successors=0),scope='Fixed eleven vertices; legal Hamiltonian-path complements; optional free relabellings/converse. No state-space BFS.')
print(json.dumps(report,indent=2,sort_keys=True))
