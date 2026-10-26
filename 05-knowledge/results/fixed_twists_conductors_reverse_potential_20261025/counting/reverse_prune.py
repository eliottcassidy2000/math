# Fixed-precision pruning: modulo 2^K an edge with reduced cost J' >= K contributes nothing; drop it and all states that no longer
# reach the terminal. Measure the surviving graph at larger N and verify the K-digit result is unchanged.
import json, time, sys
from reverse_potential import build_graph, ordinary_direct, oddfact
def analyze(n, Ks=(1,4,8,16,32,64)):
    t0=time.time()
    vertices, edges, initial, terminal = build_graph(n)
    out={}
    for (s,d,T,J) in edges: out.setdefault(s,[]).append((d,T,J))
    V={terminal:0}
    for v in sorted(vertices, key=lambda v:v[0]):
        if v!=terminal: V[v]=min(J+V[w] for (w,T,J) in out[v])
    B=V[initial]
    red=[(v,w,T,J+V[w]-V[v]) for v in vertices if v!=terminal for (w,T,J) in out[v]]
    naive_negative = sum(1 for (w,T,J) in out[initial] if J - B < 0)   # hostile: subtracting B at the initial edge
    a=ordinary_direct(n); ON=oddfact(n); exact=(ON*a)>>B
    res={'n':n,'vertices':len(vertices),'edges':len(edges),'B':B,'naive_initial_negative_edges':naive_negative,'per_K':{}}
    for K in Ks:
        mod=1<<K
        kept=[e for e in red if e[3]<K]
        # forward reachability from initial over kept edges, backward reachability to terminal
        fw={initial}; changed=True
        adj={}
        for (v,w,T,Jp) in kept: adj.setdefault(v,[]).append((w,T,Jp))
        order=sorted(vertices,key=lambda v:-v[0])
        for v in order:
            if v in fw:
                for (w,T,Jp) in adj.get(v,[]): fw.add(w)
        bw={terminal}
        radj={}
        for (v,w,T,Jp) in kept: radj.setdefault(w,[]).append(v)
        for v in sorted(vertices,key=lambda v:v[0]):
            if v in bw:
                for u in radj.get(v,[]): bw.add(u)
        live=fw&bw
        live_edges=[e for e in kept if e[0] in live and e[1] in live]
        Z={initial:1}
        for v in order:
            if v in live and v!=terminal and Z.get(v,0):
                for (w,T,Jp) in adj.get(v,[]):
                    if w in live: Z[w]=(Z.get(w,0)+T*pow(2,Jp,mod)*Z[v])%mod
        ok = Z.get(terminal,0)==exact%mod
        res['per_K'][K]={'live_vertices':len(live),'live_edges':len(live_edges),'correct':ok}
    return res
for n in [60,80,100]:
    print(json.dumps(analyze(n)))
