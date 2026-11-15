"""Small exact controls, not a proof by growing the census."""
from itertools import permutations

def need(ok):
    if not ok:
        raise RuntimeError("control failed")

def step(f, q):
    return tuple((t*t*f[t]+t*(f[(t+1)%q]-f[t-1]))%q for t in range(q))

def table(f,q,k):
    out=[tuple(f)]
    for _ in range(k):
        out.append(step(out[-1],q))
    return out

def literal(r):
    n=len(r)
    return sum(all(r[a]>>b&1 for a,b in zip(p,p[1:]))
               for p in permutations(range(n)))

def guard(r):
    n=len(r)
    return tuple(x|(1<<(n+1)) for x in r)+((1<<n)-1,1<<n)

def poly_value(rows,t):
    if t==0:
        return int(not rows)
    n=len(rows)
    H=[1]+[0]*((1<<n)-1)
    D={}
    for m in range(1,1<<n):
        for i in range(n):
            if not m>>i&1:
                continue
            p=m^(1<<i)
            D[m,i]=sum(D[p,j] for j in range(n) if p>>j&1 and rows[j]>>i&1) if p else 1
            H[m]+=D[m,i]
    V=H
    for _ in range(1,t):
        U=[0]*len(H)
        for m in range(len(H)):
            a=m
            while True:
                U[m]+=V[a]*H[m^a]
                if a==0:
                    break
                a=(a-1)&m
        V=U
    return V[-1]

def main():
    for p,want in zip((1,3,5,7),((7,3),(5,9),(13,1),(13,1))):
        F=table([pow(t,p,16) for t in range(16)],16,6)
        need((F[4][1],F[5][1])==want)
        need(tuple((a-b)%16 for a,b in zip(F[6],F[4]))==tuple(8*(t%2) for t in range(16)))
        print("seed",p,want)
    B=(24,1,3,6,14)
    need(literal(B)==poly_value(B,1)==11)
    for core in [(0,),B]:
        for k in range(3):
            next_core=guard(core)
            need(poly_value(next_core,1)==poly_value(core,1)+poly_value(core,2))
            print("literal",len(core),k,poly_value(core,1),poly_value(next_core,1))
            core=next_core
    for r in range(4,10):
        q=1<<r; P=q//4
        for p in (1,3,5,7,9,75):
            F=table([pow(t,p,q) for t in range(q)],q,r+2*P)
            A=[v[1] for v in F[r:r+P]]
            need(sorted(A)==list(range(A[0]%4,q,4)))
            need(A==[v[1] for v in F[r+P:r+2*P]])
        q*=9
        found=set()
        for core in [(0,), B]:
            Pvals=[poly_value(core,t) for t in range(len(core)+1)]
            # Pad directly by multiplying the recovered count at each
            # residue. Newton forward interpolation uses exact binomials.
            from math import comb
            coeff=[]
            while Pvals:
                coeff.append(Pvals[0])
                Pvals=[b-a for a,b in zip(Pvals,Pvals[1:])]
            def eval_p(t):
                return sum(c*comb(t,i) for i,c in enumerate(coeff) if i<=t)%q
            values=[eval_p(t) for t in range(q)]
            for j in (1,2,3):
                init=[p*pow(t,4*j,q)%q for t,p in enumerate(values)]
                F=table(init,q,r+P)
                found.update(v[1] for v in F[r:r+P])
        need(found==set(range(3,q,6)))
        print("depth modulus",1<<r,"joint residue count",len(found))
    print('PASS')
if __name__=='__main__':
    main()
