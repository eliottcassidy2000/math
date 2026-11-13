"""Literal controls for the one-return-star transfer and binary seam."""
from itertools import permutations, combinations, product


def ck(ok, receipt):
    if not ok:
        raise AssertionError(receipt)


def matrices(n):
    edges = list(combinations(range(n), 2))
    for z in product(range(2), repeat=len(edges)):
        yield {(u,v) if flag else (v,u) for flag,(u,v) in zip(z,edges)}


def covers(n, a, two):
    for p in permutations(range(n)):
        for i in range(n+1) if two else (n,):
            if all(e in a for q in (p[:i],p[i:]) for e in zip(q,q[1:])):
                yield p[:i], p[i:]


def literal(ni, nu, nw, i, u, w, s, b):
    a = set(i)
    a |= {(x+ni,y+ni) for x,y in u}
    a |= {(x+ni+nu,y+ni+nu) for x,y in w}
    v = ni+nu+nw
    bounds = (0,ni,ni+nu,v)
    for l in range(3):
        for r in range(l+1,3):
            a |= {(x,y) for x in range(bounds[l],bounds[l+1]) for y in range(bounds[r],bounds[r+1])}
    a |= {(v,x) if x in s else (x,v) for x in range(ni)}
    a |= {(x,v) for x in range(ni,ni+nu)}
    a |= {(x+ni+nu,v) if x in b else (v,x+ni+nu) for x in range(nw)}
    end = [[0]*(v+1) for _ in range(1<<(v+1))]
    for j in range(v+1):
        end[1<<j][j] = 1
    for mask in range(1,1<<(v+1)):
        for j in range(v+1):
            if end[mask][j]:
                for k in range(v+1):
                    if not (mask>>k&1) and (j,k) in a:
                        end[mask|1<<k][k] += end[mask][j]
    return sum(end[-1])


def L(b):
    v = 0
    while b:
        b, bit = divmod(b,2)
        if bit:
            v += b+(b%2)
    return v


def run():
    cases = 0
    for ni,nu,nw in product(range(1,4),range(1,4),range(3)):
        for i,u,w in product(matrices(ni),matrices(nu),matrices(nw)):
            hi = sum(1 for _ in covers(ni,i,False))
            hw = sum(1 for _ in covers(nw,w,False))
            hu = sum(1 for _ in covers(nu,u,False))
            ku = sum(1 for _ in covers(nu,u,True))
            for sb,bb in product(range(1,2**ni),range(2**nw)):
                s = {j for j in range(ni) if sb>>j&1}
                b = {j for j in range(nw) if bb>>j&1}
                ac = [(p,q) for p,q in covers(ni,i,True) if q and q[0] in s]
                alpha = len(ac)
                lam = sum(1 for p,_ in ac if p and p[-1] in s)
                beta = sum(1 for p,_ in covers(nw,w,True) if not p or p[-1] in b)
                wb = w | {(x,nw) if x in b else (nw,x) for x in range(nw)}
                d = sum(1 for _ in covers(nw+1,wb,False))
                actual = literal(ni,nu,nw,i,u,w,s,b)
                ck(actual == alpha*beta*ku+(hi*d-hw*lam)*hu,(ni,nu,nw,i,u,w,s,b))
                cases += 1
    print("literal factorization controls:",cases)
    pairs = 0
    for n in range(1,8):
        a = set(combinations(range(n),2))
        for b in range(1,2**n):
            s = {j for j in range(n) if b>>(n-1-j)&1}
            ac = [(p,q) for p,q in covers(n,a,True) if q and q[0] in s]
            ck(len(ac)==b,(n,b,"alpha"))
            ck(sum(1 for p,_ in ac if p and p[-1] in s)==L(b),(n,b,"lambda"))
            pairs += 1
    print("literal transitive mask controls:",pairs)
    for b in range(4096):
        for r in range(1,10):
            q,l = divmod(b,2**r)
            rev = int(f'{l:0{r}b}'[::-1],2)
            ck(L(b)==L(q)+q*rev+L(l)+(q%2)*(l//(2**(r-1))),(b,r))
    print("exact radix seam controls:",4096*9)


if __name__ == "__main__":
    run()
