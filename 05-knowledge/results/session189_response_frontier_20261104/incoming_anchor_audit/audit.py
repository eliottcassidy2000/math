#!/usr/bin/env python3
"""Independent cycle enumeration and literal Hamilton DP, no incoming code imported."""
from pathlib import Path
from itertools import permutations
from functools import lru_cache
import json, random, time, hashlib
P=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def tournament(n,code):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:rows[i]|=1<<j
            else: rows[j]|=1<<i
            k+=1
    return tuple(rows)
def op(rows):
    mask=(1<<len(rows))-1
    return tuple(x^(mask^(1<<i)) for i,x in enumerate(rows))
def tt(k):return tuple(sum(1<<j for j in range(i+1,k)) for i in range(k))
def permtt(order):
    rows=[0]*len(order)
    for i,u in enumerate(order):
        for v in order[i+1:]: rows[u]|=1<<v
    return tuple(rows)
def mul(p,var,complement=False):
    out=set(p) if complement else set()
    for m in p:
        q=m|var
        if q in out:out.remove(q)
        else:out.add(q)
    return out

@lru_cache(maxsize=4096)
def cycle_basis(rows):
    """Enumerate simple cycles anchored at x, then at y without x.
    Store their exact ANF parity as an affine function of old vertex weights.
    Edges involving a port are independent literal incidence polynomials.
    """
    n=len(rows);x=n;y=n+1;basis=[0]*(n+1)
    def edge(p,u,v):
        if u<n and v<n:return p if rows[u]>>v&1 else set()
        if u==x and v==y:return p
        if u==y and v==x:return set()
        if v==x:return mul(p,1<<u)
        if v==y:return mul(p,1<<(n+u))
        if u==x:return mul(p,1<<v,True)
        return mul(p,1<<(n+v),True)
    def walk(anchor,u,visited,pathlen,poly):
        if pathlen>=3:
            closed=edge(poly,u,anchor)
            if closed:
                value=0
                for m in closed:value ^= 1<<m
                old=visited&((1<<n)-1)
                ports=((visited>>n)&3).bit_count()
                if ports%2:basis[0]^=value
                for i in range(n):
                    if old>>i&1:basis[i+1]^=value
        for v in range(n+2):
            if visited>>v&1 or (anchor==y and v==x):continue
            nxt=edge(poly,u,v)
            if nxt:walk(anchor,v,visited|1<<v,pathlen+1,nxt)
    walk(x,x,1<<x,1,{0});walk(y,y,1<<y,1,{0})
    return tuple(basis)
def weight_poly(basis,w):
    p=basis[0]
    for i in range(len(basis)-1):
        if w>>i&1:p^=basis[i+1]
    return p
def delta(rows,w):return weight_poly(cycle_basis(rows),w)^weight_poly(cycle_basis(op(rows)),w)
def functional_mask(rows,v):
    n=len(rows);s=set()
    for u in range(n):
        if u==v:continue
        s.add((1<<u)|(1<<(n+v)))
        if rows[v]>>u&1:s.add((1<<u)|(1<<v)|(1<<(n+v)))
    return sum(1<<m for m in s)
def L(rows,p,v):return (p&functional_mask(rows,v)).bit_count()%2

def inflate(Q,children):
    blocks=[];start=0
    for child in children:blocks.append(list(range(start,start+len(child))));start+=len(child)
    rows=[0]*start
    for i,A in enumerate(children):
        for u,v in enumerate(blocks[i]):
            for a,b in enumerate(blocks[i]):
                if A[u]>>a&1:rows[v]|=1<<b
            for j,B in enumerate(blocks):
                if Q[i]>>j&1:rows[v]|=sum(1<<b for b in B)
    return tuple(rows),blocks

def ports(rows,a,c):
    n=len(rows);r=list(rows)+[1<<(n+1),0]
    for v in range(n):
        if a>>v&1:r[v]|=1<<n
        else:r[n]|=1<<v
        if c>>v&1:r[v]|=1<<(n+1)
        else:r[n+1]|=1<<v
    return tuple(r)
@lru_cache(maxsize=30000)
def H(rows):
    n=len(rows);dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,count in dp[mask].items():
            bits=rows[v]&~mask
            while bits:
                bit=bits&-bits;bits-=bit;u=bit.bit_length()-1;t=dp[mask|bit];t[u]=t.get(u,0)+count
    return sum(dp[-1].values())
def endpoints(Q):
    n=len(Q);start=[0]*n;end=[0]*n
    for p in permutations(range(n)):
        if all(Q[p[i]]>>p[i+1]&1 for i in range(n-1)):start[p[0]]+=1;end[p[-1]]+=1
    return [(a+b)%2 for a,b in zip(start,end)]
def swap_poly(p,n,u,v):
    out=0
    while p:
        b=p&-p;p-=b;m=b.bit_length()-1
        for offset in [0,n]:
            i=u+offset;j=v+offset
            if ((m>>i)^(m>>j))&1:m^=(1<<i)|(1<<j)
        out^=1<<m
    return out

def weighted_checks():
    counts={};raw=[];start=time.time()
    for n in range(3,6):
        cnt=0
        for code in range(1<<(n*(n-1)//2)):
            rows=tournament(n,code);B=cycle_basis(rows);C=cycle_basis(op(rows));F=[functional_mask(rows,v) for v in range(n)]
            for w in range(1<<n):
                d=weight_poly(B,w)^weight_poly(C,w)
                need(all((d&F[v]).bit_count()%2==(w>>v&1) for v in range(n)),f'E1 {n} {code} {w}')
                cnt+=n
        counts[str(n)]=cnt
        print('weighted',n,cnt,round(time.time()-start,2),flush=True)
    rng=random.Random(1894864)
    for n in [6,7]:
        for z in range(40):
            rows=tournament(n,rng.randrange(1<<(n*(n-1)//2)));w=rng.randrange(1<<n);d=delta(rows,w)
            need(all(L(rows,d,v)==(w>>v&1) for v in range(n)),'random E1')
            raw.append({'rows':rows,'weight_mask':w,'L':[L(rows,d,v) for v in range(n)]})
    return {'exhaustive_functional_checks':counts,'random_records':raw}

def kernel_checks():
    counts={};witnesses=[]
    for m in range(2,6):
        count=0
        for code in range(1<<(m*(m-1)//2)):
            Q=tournament(m,code);p=endpoints(Q)
            for j in range(m):
                for wb in [0,1]:
                    # Inflation keeps source u first and beta next. All others singleton symbols.
                    Z,blocks=inflate(Q,[tt(2) if i==j else tt(1) for i in range(m)])
                    u,b=blocks[j];n=m+1
                    unchanged,_=inflate(op(Q),[tt(2) if i==j else tt(1) for i in range(m)])
                    reverse=op(Z)
                    for outside in range(1<<(m-1)):
                        w=1<<u|(wb<<b);t=0;outsideones=[]
                        for i in range(m):
                            if i==j:continue
                            if outside>>t&1:w|=1<<blocks[i][0];outsideones+=blocks[i]
                            t+=1
                        d=weight_poly(cycle_basis(Z),w)^weight_poly(cycle_basis(unchanged),w)
                        f=delta(Z,w)
                        need(d!=0 and f!=0,'kernel nonzero')
                        need(L(Z,f,u)==1,'reverse L')
                        if outsideones:
                            need(all(L(Z,d,z)==1 for z in outsideones),'outside singleton L')
                        elif wb==1:need(L(Z,d,u)==1,'pair unchanged L')
                        else:
                            extra=(1<<u)|(1<<b)|(1<<(n+b))
                            val=L(Z,d,u)^L(Z,d,b)^((d>>extra)&1)
                            need(val==(1^p[j]),'hard unchanged lambda')
                        correction=d^f
                        need(correction==swap_poly(correction,n,u,b),'swap orbit identity')
                        count+=2
        counts[str(m)]=count
        print('kernel',m,count,flush=True)
    return {'exhaustive_kernel_comparisons':counts}

def evaluate(poly,assignment):
    parity=0
    while poly:
        b=poly&-poly;poly-=b;m=b.bit_length()-1
        if m&assignment==m:parity^=1
    return parity

def physical_checks():
    rng=random.Random(4864189);records=[]
    Q=(2,4,1)
    specs=[('singleton',[tt(1),tt(3),(2,4,1)],[tt(1),(2,4,1),tt(3)],None),
           ('unchanged_TT3',[tt(3),tt(2),tt(2)],[tt(3),op(tt(2)),tt(2)],0),
           ('reversed_TT3',[tt(3),tt(2),tt(2)],[op(tt(3)),op(tt(2)),tt(2)],0),
           ('different_TT4_shared_end',[tt(4),tt(2),tt(2)],[permtt([0,2,1,3]),tt(2),tt(2)],0),
           ('unchanged_TT2',[tt(2),tt(3),tt(2)],[tt(2),(2,4,1),tt(2)],0),
           ('nontransitive_common_source',[(14,4,8,2),tt(2),tt(2)],[tt(4),tt(2),tt(2)],0)]
    for name,A,B,j in specs:
        W,blocks=inflate(Q,A);U,_=inflate(op(Q),B)
        groups=[]
        for i,block in enumerate(blocks):
            if i==j:groups.extend([[block[0]],block[1:]])
            else:groups.append(block)
        r=len(groups)
        # Literal restriction graph and correct contraction weights.
        Z=tuple(sum(1<<b for b,g in enumerate(groups) if W[group[0]]>>g[0]&1) for group in groups)
        Zp=tuple(sum(1<<b for b,g in enumerate(groups) if U[group[0]]>>g[0]&1) for group in groups)
        w=sum(1<<i for i,g in enumerate(groups) if len(g)==1)
        poly=weight_poly(cycle_basis(Z),w)^weight_poly(cycle_basis(Zp),w)
        base=(H(ports(W,0,0))-H(ports(U,0,0)))//2%2
        contrasts=[];bank=[]
        for z in range(1<<(2*r)):
            a=sum(1<<v for i,g in enumerate(groups) if z>>i&1 for v in g)
            c=sum(1<<v for i,g in enumerate(groups) if z>>(r+i)&1 for v in g)
            h1=H(ports(W,a,c));h2=H(ports(U,a,c));d=h1-h2
            need(d%2==0 and d//2%2==(base^evaluate(poly,z)),'physical contraction')
            bank.append([a,c,h1,h2]);contrasts.append(d%4)
        need(2 in contrasts,'physical separator')
        first=next(row for row in bank if (row[2]-row[3])%4==2)
        records.append({'name':name,'source':W,'target':U,'groups':groups,'weights':w,'kernel_source':Z,'kernel_target':Zp,'base_half_contrast':base,'query_bank':bank,'first_separator':first})
        print('physical',name,len(bank),first,flush=True)
    # Ambient lifts use a directed triangle quotient containing the whole old core as a module.
    ambient=[]
    for record in records[:4]:
        W=tuple(record['source']);U=tuple(record['target']);a,c,_,_=record['first_separator'];n=len(W)
        AW,_=inflate(Q,[W,tt(1),tt(1)]);AU,_=inflate(Q,[U,tt(1),tt(1)])
        fixed=(1<<(n+1)) # last outside vertex points to the old module and both ports
        aa=a|fixed;cc=c|fixed
        x=H(ports(AW,aa,cc));y=H(ports(AU,aa,cc));need((x-y)%4==2,'ambient mod4 lift')
        ambient.append({'name':record['name'],'source':AW,'target':AU,'a':aa,'c':cc,'H_source':x,'H_target':y,'difference':x-y})
    # Exact fatal missing-hypothesis counterexample.
    failure={'m':1,'source':tt(2),'target':tt(2),'all_query_differences':[H(ports(tt(2),a,c))-H(ports(tt(2),a,c)) for a in range(4) for c in range(4)]}
    return {'physical_records':records,'ambient_controls':ambient,'m1_counterexample':failure}

def main():
    start=time.time();out={'status':'PASS','weighted':weighted_checks(),'kernels':kernel_checks(),'physical':physical_checks()}
    out['elapsed_seconds']=round(time.time()-start,3);out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (P/'controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','seconds':out['elapsed_seconds'],'weighted':out['weighted']['exhaustive_functional_checks'],'kernel':out['kernels'],'physical_query_comparisons':sum(len(r['query_bank']) for r in out['physical']['physical_records']),'ambient':out['physical']['ambient_controls']}),flush=True)
if __name__=='__main__':main()
