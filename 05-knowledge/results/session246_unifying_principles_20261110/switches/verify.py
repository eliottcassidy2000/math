"""Standalone exact shared-core cut-exchange certificate; no repo imports."""
import itertools,json,sys

def need(ok,label):
    if not ok: raise RuntimeError(label)

def reverse(a):
    return [((1<<len(a))-1)^row^(1<<i) for i,row in enumerate(a)]

def tables(a):
    n=len(a); full=1<<n
    end=[[0]*full for _ in a]; h=[0]*full; h[0]=1
    for mask in range(1,full):
        for j in range(n):
            if not(mask>>j&1):continue
            old=mask^(1<<j)
            end[j][mask]=1 if old==0 else sum(end[i][old] for i in range(n) if old>>i&1 and a[i]>>j&1)
            h[mask]+=end[j][mask]
    return h,end

def coordinates(a):
    h,end=tables(a);_,start=tables(reverse(a));n=len(a);full=(1<<n)-1
    av=[sum(h[m]*start[i][full^m] for m in range(full+1)) for i in range(n)]
    W=[[0]*n for _ in a]
    for i in range(n):
        for j in range(i):
            W[i][j]=W[j][i]=sum(end[i][m]*start[j][full^m]+end[j][m]*start[i][full^m] for m in range(full+1))
    return h[-1],av,W

def closure(a,cut):
    n=len(a);g=[row|(1<<n) for row in a]+[0]
    return [row|((1<<(n+1)) if not(cut>>i&1) else 0) for i,row in enumerate(g)]+[cut]

def strong(a):
    for i in range(len(a)):
        seen=1<<i
        while True:
            new=seen
            for j in range(len(a)):
                if seen>>j&1:new|=a[j]
            if new==seen:break
            seen=new
        if seen!=(1<<len(a))-1:return False
    return True

def scalar(h,av,W,cut):
    ids=[i for i in range(len(av)) if cut>>i&1]
    return h+2*sum(av[i] for i in ids)-sum(W[i][j] for i,j in itertools.combinations(ids,2))

def complete(weights):
    s=0
    for w in sorted(w for w in weights if w):
        if w>s+1:return False
        s+=w
    return True

def select(weights,target):
    bits=[0]*len(weights)
    for i in sorted(range(len(weights)),key=lambda i:weights[i],reverse=True):
        if weights[i] and target>=weights[i]:bits[i]=1;target-=weights[i]
    need(target==0,'bounded selector remainder')
    return bits

def value(u,v,b):
    return sum(x*y for x,y in zip(u,b))+sum(c*b[i]*b[j] for (i,j),c in v.items())

def lip(u,v):
    return all(u[i]+sum(min(0,c) for e,c in v.items() if i in e)>=-1 and u[i]+sum(max(0,c) for e,c in v.items() if i in e)<=1 for i in range(len(u)))

def extremum(u,v,maximize):
    """Path/cycle component DP, returning value and an attaining bit vector."""
    n=len(u);adj=[[] for _ in u]
    for (i,j),c in v.items():
        if c:adj[i].append(j);adj[j].append(i)
    need(all(len(z)<=2 for z in adj),'interaction degree')
    unseen=set(range(n));answer=[0]*n;total=0
    pick=max if maximize else min
    while unseen:
        root=min(unseen);comp={root};stack=[root]
        while stack:
            for j in adj[stack.pop()]:
                if j not in comp:comp.add(j);stack.append(j)
        unseen-=comp
        ends=[i for i in comp if len(adj[i])<2]
        first=min(ends or comp);order=[first];previous=None;at=first
        while True:
            nxt=[j for j in adj[at] if j!=previous and j!=first]
            if not nxt:break
            previous,at=at,nxt[0];order.append(at)
        need(set(order)==comp,'component traversal')
        cyc=not ends
        options=[]
        for initial in (0,1):
            dp={initial:(u[first]*initial,[initial])}
            for prev,cur in zip(order,order[1:]):
                edge=v[tuple(sorted((prev,cur)))];nd={}
                for b in (0,1):
                    nd[b]=pick((old+u[cur]*b+edge*a*b,bits+[b]) for a,(old,bits) in dp.items())
                dp=nd
            for last,(old,bits) in dp.items():
                if cyc:old+=v[tuple(sorted((order[-1],first)))]*last*initial
                options.append((old,bits))
        best,bits=pick(options);total+=best
        for i,b in zip(order,bits):answer[i]=b
    need(value(u,v,answer)==total,'DP attained extremum')
    return total,answer

def path_select(u,v,target):
    lo,b=extremum(u,v,False);hi,end=extremum(u,v,True)
    need(lo<=target<=hi,'path target range')
    if value(u,v,b)==target:return b
    for i in range(len(u)):
        if b[i]!=end[i]:
            old=value(u,v,b);b[i]=end[i];new=value(u,v,b)
            need(abs(new-old)<=1,'unit path step')
            if new==target:return b
    raise RuntimeError('path failed to hit target')

def main():
    if '--fail' in sys.argv:need(False,'intentional sentinel')
    core=[10,28,49,36,41,3];h,av,W=coordinates(core)
    need((h,av)==(43,[55,64,72,55,78,57]),'literal core scalar and fields')
    expected=[[0,26,14,20,26,28],[26,0,20,18,20,20],[14,20,0,26,28,24],[20,18,26,0,24,22],[26,20,28,24,0,22],[28,20,24,22,22,0]]
    need(W==expected,'literal interaction matrix')
    for cut in range(64):need(tables(closure(core,cut))[0][-1]==scalar(h,av,W,cut),'all-cut independent polynomial')
    # Check the positive lower bound using independent deletion graphs.
    for i,j in itertools.combinations(range(6),2):
        ids=[k for k in range(6) if k not in(i,j)]
        deleted=[sum((core[u]>>w&1)<<k for k,w in enumerate(ids)) for u in ids]
        need(W[i][j]>=2*tables(deleted)[0][-1]>0,'OCF endpoint lower bound')
    pairs=[(2,4),(5,3)];base=37;rows=[]
    need(W[4][3]-W[4][5]-W[2][3]+W[2][5]==0,'shared-core orthogonality')
    need([-2*W[p][q] for p,q in pairs]==[-56,-44],'negative Gram diagonal')
    for bits in itertools.product((0,1),repeat=2):
        cut=base
        for bit,(p,q) in zip(bits,pairs):
            if bit:cut^=(1<<p)|(1<<q)
        graph=closure(core,cut);H,end=tables(graph);_,start=tables(reverse(graph));full=(1<<len(graph))-1
        perm=sum(all(graph[i]>>j&1 for i,j in zip(order,order[1:])) for order in itertools.permutations(range(len(graph))))
        U=sum(H[m]*start[-1][full^m] for m in range(full+1))
        need(H[-1]==perm==345+2*sum(bits),'literal DP and permutations')
        need(U==411+12*bits[0]-4*bits[1],'retained response sidecar')
        need(end[-1][-1]==43 and strong(graph),'root endpoint and strongness')
        rows.append(dict(bits=bits,cut=cut,count=H[-1],response=U,root_end=end[-1][-1],rows=graph))
    weight_cases=0
    for k in range(5):
        for ws in itertools.product(range(5),repeat=k):
            sums={sum(w*b for w,b in zip(ws,bs)) for bs in itertools.product((0,1),repeat=k)}
            need(complete(ws)==(sums==set(range(sum(ws)+1))),'complete sequence iff')
            if complete(ws):
                for target in sums:need(sum(w*b for w,b in zip(ws,select(ws,target)))==target,'complete selector')
            weight_cases+=1
    need(not complete([1,4]),'primitive gcd hostile')
    lip_cases=accepted=0;edges=list(itertools.combinations(range(3),2))
    for us in itertools.product(range(-2,3),repeat=3):
        for cs in itertools.product(range(-2,3),repeat=3):
            vs=dict(zip(edges,cs));vals={bs:value(us,vs,bs) for bs in itertools.product((0,1),repeat=3)}
            brute=all(abs(f-vals[tuple(1-b if i==j else b for i,b in enumerate(bs))])<=1 for bs,f in vals.items() for j in range(3))
            need(lip(us,vs)==brute,'exact marginal criterion')
            if brute:
                accepted+=1
                need(all(sum(abs(c) for e,c in vs.items() if i in e)<=2 for i in range(3)),'weighted degree capacity')
                lo,_=extremum(us,vs,False);hi,_=extremum(us,vs,True)
                need((lo,hi)==(min(vals.values()),max(vals.values())),'component DP versus exhaustive')
                need(set(vals.values())==set(range(lo,hi+1)),'no gaps')
                for target in range(lo,hi+1):need(value(us,vs,path_select(us,vs,target))==target,'physical-coordinate path selector')
            lip_cases+=1
    # Growing formal path family checks the all-size DP. Physical embedding NOT claimed.
    u=[1]*100;v={(i,i+1):-1 for i in range(99)}
    need(lip(u,v),'formal run chain criterion')
    need(extremum(u,v,False)[0]==0 and extremum(u,v,True)[0]==50,'formal run chain extrema')
    for target in range(51):need(value(u,v,path_select(u,v,target))==target,'formal run chain selector')
    print(json.dumps(dict(status='PASS',physical_cube=rows,all_cuts=64,permutations_per_graph=40320,complete_sequence_cases=weight_cases,lipschitz_coefficient_cases=lip_cases,lipschitz_accepted=accepted,formal_chain_length=100,scope='One literal shared-core cube; proved conditional interval and capacity criteria; no unrestricted coverage'),indent=2))

if __name__=='__main__':main()
