"""Independent arithmetic/physical audit; the imported constructor is the SUT.
The evaluations below do not import its recurrence, graph or profile counter.
"""
import hashlib
import importlib.util
import itertools
import json
import pathlib
import random

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('candidate', HERE/'source'/'verify.py')
candidate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(candidate)

def require(ok, why):
    if not ok:
        raise RuntimeError(why)

def evaluate(n):
    # Direct sum over the zero positions of the final integer. No recursion.
    return 1 + sum((n >> (i+1)) + ((n >> (i+1)) & 1)
                   for i in range(n.bit_length()) if not (n >> i & 1))

def Q(n):
    s = n.bit_length()
    return 1 + sum(1 << (s-i-1) for i in range(s) if not (n >> i & 1))

def actual_graph(L, strong=False, sinks=0):
    n=L-1; s=n.bit_length(); root=s
    rows=[0]*(s+1)
    for i in range(s):
        for j in range(i+1,s): rows[i] |= 1 << j
        if n >> i & 1: rows[i] |= 1 << root
        else: rows[root] |= 1 << i
    if strong:
        # Put the exceptional new vertex last, unlike the source constructor.
        a=len(rows)
        rows[root] |= 1 << a
        rows.append(((1 << a)-1) ^ (1 << root))
    for _ in range(sinks):
        a=len(rows)
        rows=[row | (1 << a) for row in rows]+[0]
    return rows,root

def literal_profile(rows,r):
    # Enumerate every simple directed word, including proper subsets.
    size=1 << len(rows); full=size-1
    h=[0]*size; starts=[0]*size; ends=[0]*size
    cycle=0
    def walk(mask,last,first):
        nonlocal cycle
        h[mask]+=1
        starts[mask]+=first==r
        ends[mask]+=last==r
        if mask==full and first==r and rows[last] >> r & 1: cycle+=1
        unused=rows[last] & (full^mask)
        while unused:
            bit=unused & -unused; unused-=bit
            walk(mask|bit,bit.bit_length()-1,first)
    for a in range(len(rows)): walk(1 << a,a,a)
    h[0]=1
    R=sum(starts[m]*h[full^m] for m in range(size))
    L=sum(ends[m]*h[full^m] for m in range(size))
    return (h[full],ends[full],starts[full],R,L,h[full^(1<<r)],cycle)

def slot_count(n):
    # A Hamilton word in a one-apex graph is an increasing base subset,
    # the apex, then the increasing complementary subset.
    s=n.bit_length(); full=(1<<s)-1; result=0
    for before in range(1<<s):
        after=full^before
        if before and not (n >> (before.bit_length()-1) & 1): continue
        if after and (n >> ((after & -after).bit_length()-1) & 1): continue
        result+=1
    return result

def main():
    rng=random.Random(1972026); digest=hashlib.sha256()
    counts=dict(joint=0,neutral=0,loop=0,suffix=0,slot=0,physical=0,sinks=0,boundary=0)
    for n in range(4096):
        require(evaluate(n)==slot_count(n),'literal two-slot count')
        counts['slot']+=1
    for k in range(1,11):
        M=1<<k
        triples=itertools.product(range(M),range(1,M,2),range(1,M,2)) if k<=4 else (
            (rng.randrange(M),rng.randrange(M)|1,rng.randrange(M)|1) for _ in range(500))
        for u,z,q in triples:
            w=candidate.construct(k,u,z,q)
            n=int(w,2); f=evaluate(n)
            require(w[0]=='1' and len(w)==n.bit_length(),'minimal width')
            require((n+1)%M==u and f%M==z and Q(n)%M==q,'requested triple')
            require(len(w)<=4*k*k+3*k and (1<<len(w))%M==0,'order bound')
            digest.update(f'{k},{u},{z},{q},{n},{f},{Q(n)}\n'.encode())
            counts['joint']+=1
            # Insert neutral loops into the already constructed word, without
            # asking the constructor to recompute its correction.
            for j in (1,3):
                ww=w[:2*k]+('1'+'0'*k)*j+w[2*k:]
                nn=int(ww,2)
                require(nn>n and ((nn+1)%M,evaluate(nn)%M,Q(nn)%M)==(u,z,q),'unchanged neutral insertion')
                counts['neutral']+=1
        for t in range(8):
            n=t*M
            for a in range(1,k+4):
                nn=((n<<a)+(1<<a)-1)<<k
                require(nn%M==0 and (evaluate(nn)-evaluate(n))%M==(2-(1<<a))%M,'loop from arbitrary reset state')
                counts['loop']+=1
            for r in range(M):
                nn=(n<<k)+r
                require((evaluate(nn)-evaluate(n))%M==(evaluate(r)-1)%M,'suffix from arbitrary reset state')
                counts['suffix']+=1
    for k in (16,32,64,128):
        M=1<<k
        for u,z,q in ((0,1,1),(M-1,M-1,M-1),(M//2,12345,54321)):
            w=candidate.construct(k,u,z,q); n=int(w,2)
            require(((n+1)%M,evaluate(n)%M,Q(n)%M)==(u,z,q),'large integer triple')
            require(len(w)<=4*k*k+3*k,'large order bound')
            counts['joint']+=1
    for L in range(2,258):
        n=L-1; g=evaluate(n); q=Q(n); B=1<<n.bit_length()
        rows,r=actual_graph(L)
        p=literal_profile(rows,r)
        require(p[:6]==(g,1,L%2,q,L,1),'actual gadget paths and covers')
        rows,r=actual_graph(L,True)
        p=literal_profile(rows,r)
        require(p[:4]==(g+L-L%2,1,1,q+B) and p[5:]==(1,1),'actual strong source and unique cycle')
        counts['physical']+=2
        if L<=65:
            for t in (1,2):
                rows,r=actual_graph(L,True,t)
                pp=literal_profile(rows,r)
                require(pp[:4]==(p[0],0,1,(1<<t)*p[3]) and pp[5:]==(1,0),'actual weak sink source')
                counts['sinks']+=1
    # Explicitly reject extending the strong corollary to the excluded L=1.
    rows,r=actual_graph(1,True)
    require(literal_profile(rows,r)==(1,0,1,2,1,1,0),'excluded transitive-edge boundary')
    counts['boundary']+=1
    print(json.dumps({'status':'PASS','counts':counts,'joint_digest':digest.hexdigest(),
        'scope':'Independent proof audit plus direct integer evaluation, literal base slots, and directed-word counting; constructor imported only as the system under test.'},sort_keys=True,indent=2))

if __name__=='__main__': main()
