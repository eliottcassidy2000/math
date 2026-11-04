import itertools, subprocess, random
from lib import subst, TT, conv, random_tournament, coeffs, is_prime
def ham(A,S):
    S=list(S)
    return sum(1 for p in itertools.permutations(S) if all(A[p[i]][p[i+1]] for i in range(len(S)-1))) if S else 1
def profile(A):
    k=len(A); q=[0]*(k+1)
    for labels in itertools.product(range(k),repeat=k):
        h=max(labels)+1
        if set(labels)!=set(range(h)): continue
        pr=1
        for b in range(h): pr*=ham(A,[i for i in range(k) if labels[i]==b])
        q[h]+=pr
    return q[1:]
def F_runword(jobs):
    inp=[str(len(jobs))]
    for (Q,j,named,blocks,abits,cbits) in jobs:
        m=len(Q); bl=[i for i in range(m) if i!=j]; B=len(bl); k=len(named); S=k+2
        inp.append('%d %d'%(B,S))
        for i in bl:
            pr=profile(blocks[i]); inp.append(' '.join(map(str,[len(blocks[i])]+pr)))
        # atoms: blocks bl[0..B-1], named 0..k-1, x, y
        A=B+S; X=B+k; Y=B+k+1
        def quot(a):
            return bl[a] if a<B else j
        M=[['0']*A for _ in range(A)]
        for a in range(A):
            for b in range(A):
                if a==b: continue
                if a<B+k and b<B+k:
                    qa,qb=quot(a),quot(b)
                    if qa!=qb: M[a][b]='1' if Q[qa][qb] else '0'
                    else: M[a][b]='1' if named[a-B][b-B] else '0'
                elif a==X and b==Y: M[a][b]='1'
                elif a==Y and b==X: M[a][b]='0'
                else:
                    # port vs old atom
                    port,old,out_of_port=(a,b,True) if a in (X,Y) else (b,a,False)
                    bits=abits if port==X else cbits
                    beats_port=bits[old]   # old -> port
                    M[a][b]='1' if (beats_port and not out_of_port) or ((not beats_port) and out_of_port) else '0'
        inp+= [''.join(r) for r in M]
    out=subprocess.run(['./runword.exe'],input='\n'.join(inp)+'\n',capture_output=True,text=True,check=True).stdout.split()
    return list(map(int,out))
if __name__=='__main__':
    # validation vs direct named DP on small instances
    rng=random.Random(2); ok=0; tot=0
    for t in range(12):
        m=4; Q=random_tournament(m,rng); j=rng.randrange(m)
        named=random_tournament(rng.choice([2,3]),rng)
        blocks={i:random_tournament(rng.choice([1,2,3]),rng) for i in range(m) if i!=j}
        kids=[named if i==j else blocks[i] for i in range(m)]
        T,blk=subst(Q,kids); n=len(T)
        bl=[i for i in range(m) if i!=j]; B=len(bl)
        abits=[rng.randint(0,1) for _ in range(B+len(named))]; cbits=[rng.randint(0,1) for _ in range(B+len(named))]
        # expand to physical query and count directly
        a=[0]*n; c=[0]*n
        for ai,i in enumerate(bl):
            for v in blk[i]: a[v]=abits[ai]; c[v]=cbits[ai]
        for r,v in enumerate(blk[j]): a[v]=abits[B+r]; c[v]=cbits[B+r]
        Cd=coeffs([T])[0]
        sa=[v for v in range(n) if a[v]]; sc=[v for v in range(n) if c[v]]
        Ss=[()]+[(v,) for v in sa]+list(itertools.combinations(sa,2)); Rs=[()]+[(v,) for v in sc]+list(itertools.combinations(sc,2))
        direct=sum(Cd[(S,R)] for S in Ss for R in Rs)
        rwv=F_runword([(Q,j,named,blocks,abits,cbits)])[0]
        tot+=1; ok+= direct==rwv
    print('runword validation',ok,'/',tot)
