from itertools import combinations, product
from random import Random

def bits(S):
    while S:
        b = S & -S
        yield b.bit_length()-1
        S ^= b

def tournament(n, k):
    a = [0]*n
    for j,(u,v) in enumerate(combinations(range(n),2)):
        if k >> j & 1: u,v = v,u
        a[u] |= 1 << v
    return a

def leaf(a, p):
    return [r | ((i!=p)<<len(a)) for i,r in enumerate(a)] + [1<<p]

def data(a):
    n = len(a); M = 1 << n
    P = [[[0]*n for _ in range(n)] for _ in range(M)]
    H = [0]*M; H[0] = 1
    for s in range(1,M):
        for u in bits(s):
            t = s ^ (1<<u)
            if not t: P[s][u][u] = 1
            else:
                for v in bits(t & a[u]):
                    for w,x in enumerate(P[t][v]): P[s][u][w] += x
        H[s] = sum(map(sum,P[s]))
    return P,H

def check(x, y, tag):
    if x != y: raise RuntimeError((tag,x,y))

rng = Random(645)
controls = 0
for n in range(2,6):
    width = n*(n-1)//2
    seeds = range(1<<width) if n < 5 else rng.sample(range(1<<width), 24)
    for code in seeds:
        A = tournament(n,code)
        P,H = data(A)
        f = len(H)-1
        M = [[sum(row) for row in z] for z in P]
        E = [sum(M[s][q]*H[f^s] for s in range(1,f+1) if not s.bit_count()%2) for q in range(n)]
        L = [data(leaf(A,q))[1][-1] for q in range(n)]
        for q,p in product(range(n),repeat=2):
            if p == q: continue
            tag = n,code,q,p
            W = sum((M[s][p]-P[s][p][q])*H[f^s] for s in range(1,f+1) if s.bit_count()%2)
            D = [[0,0],[0,0]]
            for S in range(1,f+1):
                if not M[S][q]: continue
                U=f^S; s=U
                while s:
                    D[S.bit_count()%2][s.bit_count()%2] += M[S][q]*M[s][p]*H[U^s]
                    s = (s-1)&U
            B=A
            g=[]; d=[]
            for k in range(1,5):
                B=leaf(B,q)
                g.append(data(leaf(B,p))[1][-1])
                d.append(g[-1]-2*data(B)[1][-1])
            check(g[0],H[-1]+2*(E[p]+E[q]+W+D[0][1]+D[1][0]+2*D[0][0]),tag)
            b,rem=divmod(g[1]-g[0],2)
            c,rem2=divmod(g[2]-g[0]-6*b,6)
            check((rem,rem2), (0,0),tag)
            check(c>0, bool(A[q] & ~(1<<p)),tag)
            check(b>=0 and c>=0,True,tag)
            check(g[3],g[0]+14*b+36*c,tag)
            check(d[2]-3*d[1]+2*d[0],6*c,tag)
            controls+=1
print(f"Exact cross-response, third-visit test, four-height reconstruction: {controls} ordered pairs.")

def stable(a, P, p):
    return all(not x for v in bits(a[p]) for x in [sum(row[v] for row in P[((1<<len(a))-1)^(1<<p)])])

stable_controls = switches = 0
for n in range(2,5):
    for code in range(1<<(n*(n-1)//2)):
        A=tournament(n,code)
        P,H=data(A)
        for q in range(n):
            B=leaf(A,q)
            R, J=data(B)
            all_end_q=sum(row[q] for row in P[-1])
            predicted={q}
            if A[q].bit_count()==1: predicted.update(bits(A[q]))
            if not all_end_q: predicted.add(n)
            check({p for p in range(n+1) if stable(B,R,p)}, predicted, (n,code,q))
            if stable(A,P,q) and J[-1]-2*H[-1] in (1,-1):
                if A[q].bit_count()==1:
                    p=next(bits(A[q]))
                    for k in range(1,4):
                        d2=data(leaf(B,p))[1][-1]-2*data(B)[1][-1]
                        old_d=A[p].bit_count()
                        check(d2,2**(k+old_d)-3,(n,code,q,p,k))
                        check(d2 in (1,-1),k==old_d==1,(n,code,q,p,k))
                        if d2==1: switches+=1
                        B=leaf(B,q)
            stable_controls+=1
print(f"Direct stable-menu replay: {stable_controls} moves; {switches} labelled stable-unit switches.")
