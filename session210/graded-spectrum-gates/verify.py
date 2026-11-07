#!/usr/bin/env python3
"""Standalone exact checks; no repository imports or external dependencies."""
import itertools, json, math

def demand(ok, msg):
    if not ok: raise RuntimeError(msg)

def tournament(n, mask):
    a=[0]*n
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        u,v=(i,j) if mask>>k&1 else (j,i)
        a[u]|=1<<v
    return tuple(a)

def paths(a):
    n=len(a); out=[0]*(1<<n); out[0]=1
    dp=[[0]*n for _ in out]
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not s>>v&1: continue
            t=s^(1<<v)
            if t: dp[s][v]=sum(dp[t][u] for u in range(n) if a[u]>>v&1)
        out[s]=sum(dp[s])
    return out

def covers(h, q):
    f=[0]*len(h); f[0]=1
    for _ in range(q):
        g=[0]*len(h)
        for s in range(len(h)):
            t=s
            while True:
                g[s]+=f[s^t]*h[t]
                if not t: break
                t=(t-1)&s
        f=g
    return f[-1]

def close(a):
    n=len(a); tail=n; head=n+1
    return tuple(x|(1<<tail) for x in a)+(1<<head,(1<<n)-1)

def cycles(a):
    # Enumerate each simple cycle at its smallest vertex: no factor-of-length ambiguity.
    n=len(a); c=[0]*(n+1)
    for k in range(3,n+1):
        for s in itertools.combinations(range(n),k):
            for rest in itertools.permutations(s[1:]):
                p=(s[0],)+rest
                if all(a[u]>>v&1 for u,v in zip(p,p[1:]+p[:1])): c[k]+=1
    return c

def gates(a):
    n=len(a)
    return [(u,v) for u in range(n) for v in range(n) if a[u]==1<<v and sum(x>>v&1 for x in a)==1]

def strong(a):
    n=len(a)
    for u in range(n):
        s=1<<u
        while True:
            t=s
            for v in range(n):
                if s>>v&1: t |= a[v]
            if t==s: break
            s=t
        if s!=(1<<n)-1: return False
    return True

def perm_h(a):
    n=len(a)
    return sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(range(n)))

def apply_L(f):
    # Monomial coefficients, exact binomial expansion.
    g=[0]*(len(f)+2)
    for j,c in enumerate(f):
        g[j+2]+=c
        for k in range(j+1):
            g[k+1]+=c*math.comb(j,k)*(1-(-1)**(j-k))
    return g

def evaluate(f,x): return sum(c*x**i for i,c in enumerate(f))

def recover_monic_parity(n, values):
    from fractions import Fraction as F
    d=(n-1)//2
    demand(len(values)>=d,'too few samples')
    # Convert closure samples to P(1)..P(d) by triangular Laurent-step rows.
    evals=[]
    rows=[{1:1}]
    for r in range(1,d):
        nxt={}
        for k,c in rows[-1].items():
            for j,z in ((k+1,k),(k,k*k),(k-1,-k)):
                if j>0: nxt[j]=nxt.get(j,0)+c*z
        rows.append(nxt)
    for r in range(d):
        row=rows[r]
        demand(row[r+1]==math.factorial(r),'triangular coefficient')
        evals.append(F(values[r]-sum(c*evals[k-1] for k,c in row.items() if k<=r),row[r+1]))
    exps=list(range(1 if n%2 else 2,n,2))
    mat=[[F((j+1)**e) for e in exps]+[evals[j]-F((j+1)**n)] for j in range(d)]
    for i in range(d):
        p=next(k for k in range(i,d) if mat[k][i]); mat[i],mat[p]=mat[p],mat[i]
        z=mat[i][i]; mat[i]=[v/z for v in mat[i]]
        for k in range(d):
            if k!=i:
                z=mat[k][i]; mat[k]=[u-z*v for u,v in zip(mat[k],mat[i])]
    f=[F(0)]*(n+1); f[n]=F(1)
    for e,row in zip(exps,mat): f[e]=row[-1]
    demand(all(x.denominator==1 for x in f),'noninteger recovered polynomial')
    return [int(x) for x in f]

def main():
    report={'schema':'gate-correspondence-v1','construction_universe':{},'target_universe':{}}
    collision=None; seen={}; checks=0
    for n in range(1,6):
        spectrum=set()
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask); h=paths(a); b=close(a); hb=paths(b)
            demand(h[-1]==perm_h(a),'path counters')
            cb=cycles(b); ca=cycles(a)
            demand(cb[-1]==h[-1],'cycle/path correspondence')
            demand(gates(b)==[(n,n+1)] if n>=2 else len(gates(b))==3,'gate uniqueness')
            demand(strong(b),'closure strong')
            for k in range(3,n+3):
                demand(cb[k]==(ca[k] if k<=n else 0)+sum(h[s] for s in range(1,1<<n) if s.bit_count()==k-2),'cycle deck law')
            for q in range(1,5):
                lhs=covers(hb,q); rhs=q*(covers(h,q+1)+q*covers(h,q)-covers(h,q-1))
                demand(lhs==rhs,'cover operator'); checks+=1
            # Alternating convolution proves parity; checked on every induced subset.
            for s in range(1,1<<n):
                demand(sum((-1)**t.bit_count()*h[t]*h[s^t] for t in range(1<<n) if t&s==t)==0,'sign involution identity')
            d=(n-1)//2; vals=[]; z=a
            for _ in range(d): vals.append(paths(z)[-1]); z=close(z)
            f=recover_monic_parity(n,vals)
            demand(all(evaluate(f,q)==covers(h,q) for q in range(6)),'response reconstruction')
            # Independent permutation forward-adjacency distribution and its binomial transform.
            dist=[0]*n
            for word in itertools.permutations(range(n)):
                fwd=sum(a[u]>>v&1 for u,v in zip(word,word[1:]))
                dist[fwd]+=1
            demand(all(covers(h,q)==sum(c*math.comb(q+fwd,n) for fwd,c in enumerate(dist) if q+fwd>=n) for q in range(6)),'binomial descent transform')
            demand(hb[-1]==(n+2)*h[-1]+(dist[n-2] if n>=2 else 0),'one-defect closure bijection')
            spectrum.add(h[-1])
            key=(n,h[-1]); data=(covers(h,2),mask,a,f)
            if key in seen and seen[key][0]!=data[0] and collision is None:
                prev=seen[key]
                collision={'order':n,'H':h[-1],'first':{'mask':prev[1],'rows':prev[2],'P2':prev[0],'polynomial':prev[3]},'second':{'mask':mask,'rows':a,'P2':data[0],'polynomial':f}}
            seen[key]=data
        report['construction_universe'][str(n)]={'labelled':1<<(n*(n-1)//2),'spectrum':sorted(spectrum)}
    report['cover_operator_checks']=checks; report['descent_transform_cores']=1099; report['same_H_context_hostile']=collision
    # Independently exhaust target tournaments and derive their local gate class.
    hostile=None
    for N in range(3,7):
        count=0; spec=set()
        for mask in range(1<<(N*(N-1)//2)):
            a=tournament(N,mask); gs=gates(a)
            if gs:
                count+=1; spec.add(cycles(a)[-1])
            if hostile is None and strong(a) and any(x.bit_count()==1 for x in a) and any(sum(x>>v&1 for x in a)==1 for v in range(N)) and not gs:
                hc=cycles(a)[-1]
                if hc%2==0: hostile={'order':N,'mask':mask,'rows':a,'HC':hc,'outdegrees':[x.bit_count() for x in a],'indegrees':[sum(x>>v&1 for x in a) for v in range(N)]}
        expected=2 if N==3 else N*(N-1)*(1<<((N-2)*(N-3)//2))
        demand(count==expected,'target labelled enumeration')
        demand(spec==set(report['construction_universe'][str(N-2)]['spectrum']),'target spectrum')
        report['target_universe'][str(N)]={'labelled':count,'HC_spectrum':sorted(spec)}
    report['unlinked_extreme_degree_hostile']=hostile
    # Universal monomial reconstruction, beyond actual tournament realizations.
    universal=0
    for n in range(1,21):
        f=[0]*(n+1); f[n]=1
        for e in range(1 if n%2 else 2,n,2): f[e]=(-1)**e*(e+3)**2
        seq=[]; g=f
        for _ in range((n-1)//2): seq.append(evaluate(g,1)); g=apply_L(g)
        demand(recover_monic_parity(n,seq)==f,'universal polynomial reconstruction'); universal+=1
    report['universal_polynomial_reconstruction_orders']=universal
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__': main()
