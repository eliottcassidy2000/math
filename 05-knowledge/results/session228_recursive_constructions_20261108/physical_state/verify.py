#!/usr/bin/env python3
"""Standalone exact rooted-state SCC/compression checks. No repository imports."""
import argparse, hashlib, itertools, json

def check(ok, msg):
    if not ok: raise RuntimeError(msg)

def graph(n, code):
    a = [[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        a[i][j]=bool(code>>k&1); a[j][i]=not a[i][j]
    return a

def counts(a):
    n=len(a); lim=1<<n
    end=[[0]*n for _ in range(lim)]
    start=[[0]*n for _ in range(lim)]
    h=[0]*lim; h[0]=1
    for s in range(1,lim):
        for v in range(n):
            if not s>>v&1: continue
            t=s^(1<<v)
            end[s][v]=1 if not t else sum(end[t][u] for u in range(n) if a[u][v])
            start[s][v]=1 if not t else sum(start[t][u] for u in range(n) if a[v][u])
        h[s]=sum(end[s])
    return h,end,start

def state(a,r=0):
    h,e,s=counts(a); f=(1<<len(a))-1
    u=sum(s[t][r]*h[f^t] for t in range(1,f+1) if t>>r&1)
    return (h[f],e[f][r],u,h[f^(1<<r)])

def literal(a,r=0):
    n=len(a); f=(1<<n)-1; paths=[]
    for mask in range(f+1):
        paths.append([p for p in itertools.permutations([v for v in range(n) if mask>>v&1]) if all(a[x][y] for x,y in zip(p,p[1:]))])
    return (len(paths[f]),sum(p[-1]==r for p in paths[f]),sum(sum(bool(p) and p[0]==r for p in paths[t])*len(paths[f^t]) for t in range(f+1)),len(paths[f^1]))

def components(a,r=0):
    vs=[v for v in range(len(a)) if v!=r]
    reach={(i,j): (i==j or a[i][j]) for i in vs for j in vs}
    for k in vs:
        for i in vs:
            for j in vs: reach[i,j]=reach[i,j] or (reach[i,k] and reach[k,j])
    blocks=[]; unseen=set(vs)
    while unseen:
        v=min(unseen); b=sorted(w for w in unseen if reach[v,w] and reach[w,v]); blocks.append(b); unseen.difference_update(b)
    old=blocks[:]
    blocks.sort(key=lambda b: -sum(a[b[0]][c[0]] for c in old if c!=b))
    check(all(a[x][y] for i,b in enumerate(blocks) for c in blocks[i+1:] for x in b for y in c),'SCC order')
    return blocks

def induced(a,vs): return [[a[i][j] for j in vs] for i in vs]

def profile(a):
    bs=components(a); d=1; v=0; rows=[]; remove=[]; activated=False
    for b in bs:
        z=induced(a,b); h,_,starts=counts(z); f=(1<<len(b))-1
        q=sum(h[t]*h[f^t] for t in range(f+1))
        alpha=sum(starts[t][i]*h[f^t] for t in range(1,f+1) for i,x in enumerate(b) if a[0][x])
        check(q>=2*h[f] and q>=1<<len(b), 'Q floor')
        check((alpha>0)==any(a[0][x] for x in b),'alpha positivity')
        if not activated and not alpha and len(b)==1: remove+=b
        activated=activated or alpha>0
        v=q*v+alpha*d; d=h[f]*d
        rows.append([len(b),h[f],q,alpha,d,v])
    return rows,remove,d,v

def bound(u,d):
    if u==d: return 1 if d==1 else d+1
    return (0 if d==1 else d)+(u-d).bit_length()+1

def find_state(x):
    """Exact terminating decision plus labeled witness. Impractical for large bounds.
    Returns None iff no nonempty rooted tournament realizes (D,E,U,d).
    """
    if len(x)!=4 or any(type(t) is not int for t in x): return None
    D,E,U,d=x
    if min(D,U,d)<1 or E<0 or E>D or U<d: return None
    if U==d and x!=(d,d,d,d): return None
    for n in range(1,bound(U,d)+1):
        for code in range(1<<(n*(n-1)//2)):
            if state(graph(n,code))==x: return dict(n=n,code=code,root=0)
    return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--max-n',type=int,default=6); args=ap.parse_args()
    digest=hashlib.sha256(); summary=[]
    for n in range(1,args.max_n+1):
        total=1<<(n*(n-1)//2); removed=0; literal_cases=0; states=set(); active=0; worst=0
        for code in range(total):
            a=graph(n,code); x=state(a); rows,rm,d,v=profile(a)
            check((d,v)==(x[3],x[2]-x[3]),('recurrence',n,code,x,rows))
            check(v>=0,'nonnegative V')
            if v:
                first=next(i for i,z in enumerate(rows) if z[3])
                suffix=sum(z[0] for z in rows[first+1:])
                check((1<<suffix)<=v,'suffix vertex bound')
            if not v: check(x==(d,d,d,d),'sink state')
            if n<=5:
                check(x==literal(a),('literal',n,code)); literal_cases+=1
            kept=[i for i in range(n) if i not in rm]; y=state(induced(a,kept))
            check(x==y,('compression',n,code,x,y))
            check(len(kept)<=bound(x[2],d),('bound',n,code,x,len(kept)))
            check(sum(z[0] for z in rows if z[0]>1)<=(d if d>1 else 0),'nontrivial budget')
            removed+=bool(rm); active+=v>0; worst=max(worst,len(kept)); states.add(x)
            digest.update((json.dumps([n,code,x,rows,rm],separators=(',',':'))+'\n').encode())
        summary.append(dict(n=n,labeled_root0_cases=total,literal_cases=literal_cases,distinct_states=len(states),compressed_cases=removed,active_cases=active,max_compressed_order=worst))
    print(json.dumps(dict(status='FINITE-EXACT PASS; theorem pending independent audit',coverage='All labeled tournaments with root fixed at 0, hence every rooted isomorphism type through max_n',checks=['literal path/covers through n=5','SCC recurrence','alpha positivity','Q floors and suffix vertex bound','four-coordinate compression','Camion nontrivial budget','representative order bound','sink diagonal'],rows=summary,sha256=digest.hexdigest()),indent=2))
if __name__=='__main__':main()
