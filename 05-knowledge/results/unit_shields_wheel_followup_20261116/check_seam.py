from itertools import combinations
from functools import cache

SEED=(2,4,25,19,35,15)
def rev(x,n):return int(f'{x:0{n}b}'[::-1],2) if n else 0
def J(x):
    v=[i for i in range(x.bit_length()) if (x>>i)&1]
    return sum(2**(j-i-1)+(j==i+1) for i,j in combinations(v,2))
def H(b,n):return -11+82*(2**(n+1)-rev(b,n))+9*(2*b-J(b))
def pair(u,U):
    assert u>=5 and u%6==5 and 0<=U<2**u
    V=2**u;i=u+3;n=u+i+3
    K,r=divmod(62*V-13,9)
    assert r==0 and K%2==1
    X=K-U
    assert 0<=X<2**i
    L=rev(X,i)
    lo=(U<<(i+3))|(2<<i)|L
    hi=lo^(5<<i)
    return n,lo,hi
def literal(b,n):
    rows=SEED
    for e in '0'+format(b,f'0{n}b')+'0':
        k=len(rows)
        rows=tuple(x|((int(j!=3 or e=='1'))<<k) for j,x in enumerate(rows))+((0 if e=='1' else 1<<3),)
    return rows
def count_by_literal_split(rows, root=3, old=6):
    """Count from the actual adjacency, splitting a path at root."""
    from itertools import permutations
    seed=[i for i in range(old) if i!=root]
    def paths(a):
        return [w for w in permutations(a) if all(rows[x]>>y&1 for x,y in zip(w,w[1:]))]
    total=[]
    for has_before in (False,True):
        for out_allowed in (False,True):
            score=0
            for m in range(1<<len(seed)):
                B=[i for j,i in enumerate(seed) if m>>j&1]
                A=[i for j,i in enumerate(seed) if not(m>>j&1)]
                before=sum(not w or has_before or (rows[w[-1]]>>root&1) for w in paths(B))
                after=sum((rows[root]>>w[0]&1) if w else out_allowed for w in paths(A))
                score+=before*after
            total.append(score)
    F=(1<<(len(rows)-old))-1
    val=0
    for m in range(F+1):
        if m and not(rows[old+m.bit_length()-1]>>root&1):
            continue
        other=F^m
        can_out=not other or bool(rows[root]>>(old+(other&-other).bit_length()-1)&1)
        val+=total[2*bool(m)+can_out]
    return val

if __name__=='__main__':
    for u in range(5,30,6):
        V=2**u
        for U in (range(V) if V<4096 else [0,1,V//3,V//2,V-2,V-1]):
            n,b,c=pair(u,U)
            assert H(c,n)-H(b,n)==2
            a,z=literal(b,n),literal(c,n)
            assert sum((x^y).bit_count() for x,y in zip(a,z))//2==2
    n,b,c=pair(5,0)
    assert count_by_literal_split(literal(b,n))==H(b,n)
    assert count_by_literal_split(literal(c,n))==H(c,n)
    print('smallest supplied order',n+8,'addresses',b,c,'counts',H(b,n),H(c,n))
    print('verified exact +2 for the dyadic exterior examples')
