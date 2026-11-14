from itertools import combinations

def req(ok):
    if not ok:raise RuntimeError('replay mismatch')

A=(238,220,216,16,225,206,8,72)
W=((0,5),(1,2),(6,7))

def ends(A):
    n=len(A);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,1<<n):
        for i,c in enumerate(D[S]):
            if c:
                a=A[i]&~S
                while a:
                    t=a&-a;a-=t;D[S|t][t.bit_length()-1]+=c
    return D

def strong(A):
    for i in range(len(A)):
        x=1<<i
        for _ in A:
            for j,a in enumerate(A):
                if x>>j&1:x|=a
        if x!=(1<<len(A))-1:return False
    return True

def check(A,h):
    for i,j in combinations(range(len(A)),2):
        req((A[i]>>j&1)+(A[j]>>i&1)==1)
    req(strong(A) and sum(ends(A)[-1])==h)

check((36,109,8,81,103,76,5),49)
check(A,63)
check((250,136,99,148,70,90,138,52),343)
req(pow(4,21,49)==1 and (15*4**3-1)==7*137)

def cuts(A):
    n=len(A);allmask=(1<<n)-1
    E=ends(A);S=ends(tuple(allmask^a^(1<<i) for i,a in enumerate(A)))
    h=sum(E[-1]);M=[[0]*n for _ in A]
    for T in range(1,allmask):
        e=E[T];s=S[allmask^T]
        for i,x in enumerate(e):
            if x:
                for j,y in enumerate(s):
                    M[i][j]+=x*y
    ret=set()
    for U in range(1,allmask):
        val=0
        for i in range(n):
            if U>>i&1:val+=S[-1][i]
            else:
                val+=E[-1][i]
                val+=sum(z for j,z in enumerate(M[i]) if U>>j&1)
        if val==3*h:ret.add(U)
    return ret

def add(A,U):
    n=len(A)
    return tuple(a|(int(not(U>>i&1))<<n) for i,a in enumerate(A))+(U,)

def pair_cut(A,i,j):
    req((A[i]^A[j])&~((1<<i)|(1<<j))==0)
    if not(A[i]>>j&1):i,j=j,i
    return (A[i]&~(1<<j))|(1<<i)

# Interface checks in the old seed (all listed contact vertices stay singletons).
for i,j in W:req((A[i]^A[j])&~((1<<i)|(1<<j))==0)
def outs(i,j):return {z for z in range(8) if z not in (i,j) and A[i]>>z&1}
def ins(i,j):return set(range(8))-{i,j}-outs(i,j)
req(ins(0,5)=={4} and outs(6,7)=={3})
req({(x,y) for x in outs(1,2) for y in ins(1,2) if A[x]>>y&1}=={(4,0),(4,5)})

def traverse(A,remaining,depth=0):
    check(A,63*3**depth)
    expected={pair_cut(A,*p):p for p in remaining}
    req(set(expected)==cuts(A))
    return sum(traverse(add(A,u),tuple(q for q in remaining if q!=p),depth+1)
               for u,p in expected.items()) if expected else 1

req(traverse(A,W)==6)
print('49, 63, 343 strong witnesses; all six exact strong tripling histories; no fourth step')
