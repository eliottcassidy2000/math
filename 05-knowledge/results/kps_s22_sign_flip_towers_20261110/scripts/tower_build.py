"""tower_build.py - explicit tower tournaments and an independent polynomial-time path counter (kind-pasteur-2026-11-10-S22).

A tower is given by a start (n0, P0): base TT_n0 with out-set P0 (a list of positions 1..n0, possibly empty), and levels
(delta, i, m): tail TT_a, a = i + delta + m; the next closer beats exactly the tail vertices i and i + delta.
The final closing tail has t vertices and the final closer beats the last out-set.
Vertices: the chain (base, tails, final tail, in transitive order) and the closers (apexes); an earlier closer beats a
later one; a closer loses to its own tail and to all earlier chain vertices outside its out-set, and beats all later
chain vertices.
count_apex() counts Hamiltonian paths by a transfer matrix over the chain: a path is R_0 a_1 R_1 ... a_q R_q with the
apexes in some order and R_i increasing runs of the chain. It shares no code with the subset programme of towers.py."""
import itertools
def build(start, levels, t):
    n0,P0=start
    chain=[]; apex=[]; beats={}; lev=[]
    cid=n0; base=list(range(n0)); chain+=base
    O=[base[p-1] for p in P0]
    for seg in list(levels)+[None]:
        a=t if seg is None else seg[1]+seg[0]+seg[2]
        tail=list(range(cid,cid+a)); cid+=a
        v=len(apex); apex.append(v); beats[v]=set(O); lev.append((v,len(chain)+a))
        chain+=tail
        if seg is not None: d,i,m=seg; O=[tail[i-1],tail[i-1+d]]
    N=len(chain)
    for v,frm in lev: beats[v]|=set(range(frm,N))
    return N,len(apex),beats
def adjacency(N,q,beats):
    A=[[0]*(N+q) for _ in range(N+q)]
    for i in range(N):
        for j in range(i+1,N): A[i][j]=1
    for a in range(q):
        for c in range(N):
            if c in beats[a]: A[N+a][c]=1
            else: A[c][N+a]=1
        for b in range(a+1,q): A[N+a][N+b]=1
    return A
def count_apex(N,q,A):
    total=0
    for perm in itertools.permutations(range(q)):
        ap=[N+a for a in perm]
        states={tuple([0]*(q+1)):1}     # per run: 0 empty, 1 nonempty and its last vertex beats the next apex (run q: no condition), 2 nonempty and it does not
        for c in range(N):
            new={}
            for stt,cnt in states.items():
                for r in range(q+1):
                    if stt[r]==0 and r>=1 and not A[ap[r-1]][c]: continue
                    good=1 if (r==q or A[c][ap[r]]) else 2
                    s2=stt[:r]+(good,)+stt[r+1:]
                    new[s2]=new.get(s2,0)+cnt
            states=new
        for stt,cnt in states.items():
            ok=True
            for r in range(q+1):
                if stt[r]==2: ok=False; break
                if stt[r]==0 and 1<=r<=q-1 and not A[ap[r-1]][ap[r]]: ok=False; break
            if ok: total+=cnt
    return total
def start_state(start):
    """state (r, M, h, E) of (TT_n0, P0) from the closed forms: empty set, one vertex, or two vertices at distance 1 or 2."""
    n0,P0=start
    if len(P0)==0: return (1,0,1,0)
    if len(P0)==1: p=P0[0]; return (1,2**(n0-p),1,1 if p==n0 else 0)
    i,j=P0; d=j-i; assert d in (1,2); m=n0-j
    return (-1,(2**d+1)*2**m,1,1 if m==0 else 0)
def step(s,level):
    """one tower level on a state (r, M, h, E); level = (delta, i, m), delta in {1, 2}."""
    r,M,h,E=s; d,i,m=level; al=2**d+1; ga=3**d+2**d; a=i+d+m
    return (-r, al*2**m*(r+h)+2*ga*2**(i-1)*3**m*M-(E if i==1 else 0), r+2**a*M, al*2**(i-1)*M+((r-h) if m==0 else 0))
def tower_state(start, levels):
    s=start_state(start)
    for l in levels: s=step(s,l)
    return s
