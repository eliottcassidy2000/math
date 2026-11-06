import json
from itertools import combinations
from pathlib import Path

def need(x,msg):
    if not x: raise RuntimeError(msg)

def induced(a,vs):
    return [sum(((a[v]>>w)&1)<<j for j,w in enumerate(vs)) for v in vs]

def valid(a):
    n=len(a)
    need(all(not (a[i]>>i&1) and 0<=a[i]<1<<n for i in range(n)), 'diagonal/range')
    need(all(((a[i]>>j)&1)+((a[j]>>i)&1)==1 for i in range(n) for j in range(i)), 'tournament')

def reachable(a,r):
    s=1<<r
    while True:
        z=s
        for i in range(len(a)):
            if s>>i&1: z|=a[i]
        if z==s:return s
        s=z

def starts(a):
    n=len(a);N=1<<n;D=[[0]*n for _ in range(N)]
    for v in range(n):D[1<<v][v]=1
    for S in range(1,N):
        for v in range(n):
            if not S>>v&1:continue
            R=S^(1<<v)
            if R:D[S][v]=sum(D[R][w] for w in range(n) if R>>w&1 and a[v]>>w&1)
    return D[-1]

def graph(k,seed):
    core=[102,0,10,83,135,30,182,47];m=5+2*k
    slots=[];at=0
    for i in range(8):
        size=m if i==5 else 1
        slots.append(list(range(at,at+size)));at+=size
    a=[0]*at
    for i in range(8):
        for j in range(8):
            if core[i]>>j&1:
                for p in slots[i]:
                    for q in slots[j]:a[p]|=1<<q
    for i,p in enumerate(slots[5]):
        for q in slots[5][i+1:]:a[p]|=1<<q
    bits=(seed<<k)+(1<<k)-1
    for i,p in enumerate(slots[5]):
        if bits>>i&1:a[3]|=1<<p;a[p]&=~(1<<3)
        else:a[3]&=~(1<<p);a[p]|=1<<3
    return a,slots

def controls():
    checks=0
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            a=[0]*n
            for i,(x,y) in enumerate(pairs):
                if mask>>i&1:a[x]|=1<<y
                else:a[y]|=1<<x
            ss=starts(a)
            for subset in range(1,1<<n):
                vs=[i for i in range(n) if subset>>i&1];b=induced(a,vs);sb=starts(b)
                for r in range(n):
                    if reachable(a,r)!=(1<<n)-1:continue
                    # A shortest first-entry path to the kernel.
                    queue=[[r]];P=None
                    while queue:
                        p=queue.pop(0)
                        if p[-1] in vs:P=p;break
                        for j in range(n):
                            if j not in p and a[p[-1]]>>j&1:queue.append(p+[j])
                    need(P is not None,'entry path')
                    v=P[-1];need(ss[r]>=sb[vs.index(v)],'rooted insertion bound')
                    checks+=1
    return checks

rows=[];embeddings=0
for seed,want in zip([10,14,17,21],[119,127,133,141]):
    a,_=graph(0,seed);valid(a)
    need(a[1]==0,'global sink')
    b=induced(a,[i for i in range(len(a)) if i!=1]);valid(b)
    ss=starts(b)
    need(min(ss)==want,'literal minimum')
    need(all(reachable(b,r)==(1<<len(b))-1 for r in range(len(b))), 'strong kernel')
    full=starts(a)
    need(full[1]==0 and [x for i,x in enumerate(full) if i!=1]==ss,'sink path bijection')
    for k in range(21):
        big,slots=graph(k,seed);valid(big)
        vs=[]
        for i in range(8):
            if i==1:continue
            if i==5:vs.extend(slots[i][k:k+5])
            else:vs.extend(slots[i])
        need(induced(big,vs)==b,'literal retained strong kernel')
        embeddings+=1
    rows.append({'seed':seed,'rows':b,'starts':ss,'H':sum(ss),'minimum_start':min(ss)})
HERE=Path(__file__).resolve().parent
payload={'kernels':rows}
(HERE/'KERNELS.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','kernels':rows,'induced_kernel_checks':embeddings,'scales_checked':[0,20], 'rooted_insertion_controls':controls(),'scope':'Specified four-choice hereditary class only; F(8,137) remains open.'},indent=2,sort_keys=True))
