from functools import lru_cache
from itertools import combinations
from random import Random

rng=Random(31545)
cases=0
for n in range(1,7):
    pairs=list(combinations(range(n),2))
    for rep in range(16):
        rows=[0]*n
        for i,j in pairs:
            if rng.randrange(2):i,j=j,i
            rows[i]|=1<<j
        rev=[sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n)]
        @lru_cache(None)
        def paths(s,v,conv):
            a=s^(1<<v)
            return sum(paths(a,w,conv) for w in range(n)
                       if a & (rev if conv else rows)[v] & (1<<w)) if a else 1
        p=[1]+[sum(paths(s,v,False) for v in range(n) if s>>v&1)
               for s in range(1,1<<n)]
        def subs(s):
            a=s
            while a:
                yield a
                a=(a-1)&s
            yield 0
        def m(A,B):
            return [sum(A[a]*B[s^a] for a in subs(s)) for s in range(1<<n)]
        inv=[1]
        for s in range(1,1<<n):
            inv.append(-sum(p[a]*inv[s^a] for a in subs(s) if a))
        Q={0:[1]+[0]*((1<<n)-1)}
        for j in range(1,8):
            Q[j]=m(Q[j-1],p)
            Q[-j]=m(Q[-j+1],inv)
        def read(s,v,k,conv):
            return sum(paths(a,v,conv)*Q[k][s^a] for a in subs(s) if a>>v&1)
        for s in range(1,1<<n):
            sign=(-1)**(s.bit_count()+1)
            for v in range(n):
                if not s>>v&1:continue
                for k in range(7):
                    S,T=read(s,v,k,False),read(s,v,k,True)
                    if S!=sign*read(s,v,-k-1,True) or (T-sign*S)%(2*k+1):
                        raise AssertionError((rows,s,v,k))
                t,s1=read(s,v,1,True),read(s,v,1,False)
                if (t-sign*s1-read(s,v,4,True)+t)%9:
                    raise AssertionError((rows,s,v))
                cases+=1
print('all induced-root checks of reflection, odd modulus, and first mod-9 defect:',cases)
