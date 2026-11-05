#!/usr/bin/env python3
"""Independent exact certificate for the Session195 rerooted seam grammar.
Standard library only. No repository imports. Use --bound 1000 (default), or
--bound 20000 to replay the larger target-bounded arbitrary-depth closure.
"""
import argparse, hashlib, itertools, json
from collections import deque
from pathlib import Path

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def tournaments(n):
    for mask in range(1 << (n*(n-1)//2)):
        a=[0]*n; k=0
        for i in range(n):
            for j in range(i+1,n):
                if mask>>k&1: a[i]|=1<<j
                else: a[j]|=1<<i
                k+=1
        yield tuple(a)

def validate(a):
    n=len(a)
    for i in range(n):
        require(0 <= a[i] < 1<<n and not(a[i]>>i&1),'invalid row')
        for j in range(i): require(((a[i]>>j)^(a[j]>>i))&1,'invalid pair')

def dp(a):
    validate(a); n=len(a); N=1<<n
    H=[0]*N; H[0]=1
    E=[[0]*n for _ in range(N)]; S=[[0]*n for _ in range(N)]
    for v in range(n): E[1<<v][v]=S[1<<v][v]=1
    for mask in range(1,N):
        for v in range(n):
            if not mask>>v&1:continue
            rest=mask^(1<<v)
            if rest:
                E[mask][v]=sum(E[rest][u] for u in range(n) if rest>>u&1 and a[u]>>v&1)
                S[mask][v]=sum(S[rest][u] for u in range(n) if rest>>u&1 and a[v]>>u&1)
        H[mask]=sum(E[mask])
        require(H[mask]==sum(S[mask]),'forward/reverse HP mismatch')
    full=N-1; profiles=[]
    for b in range(n):
        L=sum(E[u][b]*H[full^u] for u in range(N) if u>>b&1)
        R=sum(S[u][b]*H[full^u] for u in range(N) if u>>b&1)
        deleted=full^(1<<b)
        B=sum(H[u]*H[deleted^u] for u in range(N) if u&deleted==u)
        profiles.append((H[full],S[full][b],L,E[full][b],R,B,int(n==1),H[deleted]))
    return H,profiles

def literal(a):
    n=len(a); h=0; s=[0]*n; e=[0]*n; L=[0]*n; R=[0]*n
    def path(p):return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
    for p in itertools.permutations(range(n)):
        if path(p):h+=1;s[p[0]]+=1;e[p[-1]]+=1
        for k in range(n+1):
            if path(p[:k]) and path(p[k:]):
                if k:L[p[k-1]]+=1
                if k<n:R[p[k]]+=1
    out=[]
    for b in range(n):
        vertices=[v for v in range(n) if v!=b]; B=0; d=0
        for p in itertools.permutations(vertices):
            if path(p):d+=1
            for k in range(n):
                if path(p[:k]) and path(p[k:]):B+=1
        out.append((h,s[b],L[b],e[b],R[b],B,int(n==1),d))
    return out

def signature(p):return p[:7]
def triple(p):return p[0],p[3],p[4]
def transfer(x,y):
    h,e,R=x;g,s,L,E,Q,B,delta=y[:7]
    return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R

def seam(a,r,g,b):
    n=len(a);m=len(g);G=((1<<m)-1)<<n
    out=[row|G for row in a]+[row<<n for row in g]
    out[r]^=1<<(n+b);out[n+b]|=1<<r
    return tuple(out),n+b

def sink(a):
    n=len(a);return tuple(row|(1<<n) for row in a)+(0,)

def census(bound,bank):
    seed=(1,1,1);seen={seed};q=deque([seed]);hist={};values={1};witness={1:seed}
    terminal={}
    while q:
        x=q.popleft()
        for y in bank:
            z=transfer(x,y);h,e,R=z
            require(h>=max(x[0],x[2]) and R>x[2], 'monotonicity violation')
            if h>bound:continue
            values.add(h)
            if h not in witness:witness[h]=z;terminal[z]=(x,y)
            if R<=bound and z not in seen:
                seen.add(z);q.append(z);hist[z]=(x,y)
    return seen,values,hist,witness,terminal

def witness_graph(N,hist,witness,terminal,bank):
    z=witness[N];steps=[]
    while z!=(1,1,1):
        prev,y=hist[z] if z in hist else terminal[z]
        steps.append((y,z));z=prev
    a=(0,);r=0;records=[]
    for y,z in reversed(steps):
        g,b,p=bank[y];a,r=seam(a,r,g,b)
        records.append({'gadget_rows':list(g),'gadget_root':b,'output_triple':list(z)})
    _,p=dp(a);require(triple(p[r])[0]==N,'witness target fails')
    return {'H':N,'rows':list(a),'root':r,'triple':list(triple(p[r])),'steps':records}

def test_endpoint_word(a,u,w):
    # Peer theorem: Z -> C; v beats C-u, loses to u, and beats z_i iff w_i=1.
    n=len(a);m=len(w);v=n+m;out=list(a)+[0]*(m+1)
    for i in range(m):
        z=n+i;out[z]|=(1<<n)-1
        for j in range(i+1,m):out[z]|=1<<(n+j)
        if w[i]:out[v]|=1<<z
        else:out[z]|=1<<v
    for i in range(n):
        if i==u:out[i]|=1<<v
        else:out[v]|=1<<i
    return tuple(out),v

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=1000);args=ap.parse_args()
    require(args.bound>=1,'target bound must be a positive integer')
    bank={};graphs=0;rooted=0;small=[]
    for n in range(1,6):
        for a in tournaments(n):
            H,ps=dp(a);ls=literal(a)
            require(ps==ls,'DP/literal gadget inventory mismatch')
            graphs+=1;rooted+=n
            for b,p in enumerate(ps):
                require(p[2]>=p[0] and p[4]>=p[0] and p[0]>=p[1] and p[0]>=p[3], 'profile inequalities')
                if signature(p) in bank:
                    require(bank[signature(p)][2][7]==p[7], 'deletion coordinate lost on matrix deduplication')
                bank.setdefault(signature(p),(a,b,p))
                if n<=3:small.append((a,b,p))
    require((graphs,rooted,len(bank))==(1099,5405,65),'inventory sizes differ')
    checks=0
    for a,r,p in small:
        for g,b,q in small:
            c,t=seam(a,r,g,b);_,cp=dp(c)
            require(triple(cp[t])==transfer(triple(p),q),'small exhaustive seam mismatch')
            require(cp[t][7]==p[0]*q[7],'root deletion product mismatch')
            checks+=1
    # Check every distinct matrix against nontrivial inputs and output graphs.
    controls=[((0,),0),((0,1),1),((4,1,2),0),((2,12,9,17,7),0),((24,9,3,20,6),0)]
    for a,r in controls:
        H,ps=dp(a);p=ps[r]
        for y,(g,b,q) in bank.items():
            c,t=seam(a,r,g,b);_,cp=dp(c)
            require(triple(cp[t])==transfer(triple(p),y),'bank seam mismatch')
            require(cp[t][7]==p[0]*q[7],'deletion product mismatch')
            checks+=1
    # Independently cross-check peer near-sink K and D laws, retaining old root.
    retained_checks=0
    for a,r,p in small:
        h,e,R=triple(p);d=p[7]
        c,t=seam(a,r,(0,),0);_,cp=dp(c)
        require((*triple(cp[r]),cp[r][7])==(h+R-e,d,2*R-d,d),'retained root K law')
        _,sp=dp(sink(a))
        require((*triple(sp[r]),sp[r][7])==(h,0,2*R,d),'retained sink D law')
        retained_checks+=2
    # Independently cross-check peer endpoint-one word transformation on C3.
    peer_checks=0
    a=(4,1,2);u=0;_,ps=dp(a);p=ps[u];h,s,L=p[0],p[1],p[2]
    for m in range(5):
        for w in itertools.product((0,1),repeat=m):
            if m and w[0]!=1:continue
            A=C=0;g=1;last=0
            for i,bit in enumerate(w):
                if not bit:g+=A+last
                A=2*A+bit;last=bit;C+=(1-bit)*(1<<i)
            c,v=test_endpoint_word(a,u,w);_,qs=dp(c);q=qs[v]
            if m:expected=(g*h-(1-last)*s+(A+1)*L,h,(1+C)*h+(1<<m)*L)
            else:expected=(h+L-s,h-s,h+L)
            require((q[0],q[1],q[2])==expected,'peer endpoint-word formula mismatch')
            peer_checks+=1
    # Minimal linear-state Hankel determinant, computed without external algebra.
    X=[(1,1,1),(1,0,2),(3,1,3)]
    O=[lambda x:x[0],lambda x:x[0]-x[1]+x[2],lambda x:x[0]+2*x[2]]
    M=[[f(x) for x in X] for f in O]
    det=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    require(det==-4,'minimal-state minor mismatch')
    hostile=(28,17,26,2,8);_,hostile_profiles=dp(hostile)
    require((*triple(hostile_profiles[3]),hostile_profiles[3][7])==(9,3,9,3),'content hostile fails')
    for y in bank:
        require(all(v%3==0 for v in transfer((9,3,9),y)),'content divisibility fails')
    a=(0,);r=0;triangle_states=[]
    for t in range(4):
        _,ps=dp(a);triangle_states.append(list(triple(ps[r])))
        if t:
            for v in range(len(a)):
                reach=1<<v
                while True:
                    new=reach
                    for u in range(len(a)):
                        if reach>>u&1:new|=a[u]
                    if new==reach:break
                    reach=new
                require(reach==(1<<len(a))-1,'triangle chain not strong')
        a,r=seam(a,r,(4,1,2),0)
    require(triangle_states==[[1,1,1],[5,1,7],[35,5,43],[229,35,277]],'triangle recurrence fails')
    states,values,hist,wit,terminal=census(args.bound,bank)
    misses=sorted(set(range(1,args.bound+1,2))-values)
    if args.bound>=1000:require([v for v in misses if v<=1000]==[7,21,67,163,173,247,511],'small exact closure mismatch')
    selected=[v for v in [11,13,41,43,47,127,137,139,299,983] if v<=args.bound and v in values]
    witnesses={str(v):witness_graph(v,hist,wit,terminal,bank) for v in selected}
    inventory=[{'signature':list(y),'rows':list(a),'root':b,'H_deleted_root':p[7]} for y,(a,b,p) in bank.items()]
    serialized=json.dumps(inventory,sort_keys=True,separators=(',',':')).encode()
    receipt={'status':'PASS','gadget_graphs':graphs,'rooted_gadgets':rooted,'distinct_matrices':len(bank),'graph_composition_checks':checks,'near_sink_peer_checks':retained_checks,'endpoint_word_peer_checks':peer_checks,'minimal_linear_dimension':3,'triangle_chain_triples':triangle_states,'content_hostile_state':[9,3,9,3],'Hankel_minor':M,'Hankel_determinant':det,'target_bound':args.bound,'retained_states':len(states),'attained_counts':len(values),'missing_odd_counts':misses,'inventory_sha256':hashlib.sha256(serialized).hexdigest(),'coverage_scope':'Complete arbitrary-depth closure of singleton seed under every rooted gadget of order at most five; not all tournaments.'}
    here=Path(__file__).resolve().parent
    (here/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    (here/'witnesses.json').write_text(json.dumps(witnesses,indent=2)+'\n')
    (here/f'receipt-{args.bound}.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
