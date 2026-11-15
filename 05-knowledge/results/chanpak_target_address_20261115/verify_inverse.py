from itertools import product
from inverse import *
from packets import packets, rows, dp

def trans(n):
    return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def graph(x,y,n,c,m):
    left=[int(z) for z in bin(x)[2:]]
    right=[1-int(z) for z in bin(y)[:1:-1]]
    bits=left+[m>>i&1 for i in range(n)]+right
    N=len(bits)
    a=[]
    for B in (trans(len(left)),rows(n,c),trans(len(right))):
        s=len(a);e=s+len(B);tail=((1<<N)-1)^((1<<e)-1)
        a.extend((v<<s)|tail for v in B)
    for i,t in enumerate(bits):
        if not t:a[i]|=1<<N
    a.append(sum(t<<i for i,t in enumerate(bits)))
    return a

physical=[(6,8850,54),(6,8881,54),(6,8881,53)]
for j,(n,c,m) in enumerate(physical):
    data={mm:(k,N) for k,N,mm in packets(rows(n,c))}
    assert data[m]==((45,432,258,30,129,15),K[j])
    for x,y in ((1,1),(2,3),(7,1),(1,5),(7,6)):
        assert dp(graph(x,y,n,c,m))[0][-1]==value(x,y,j)
for a,d in product((1,5,11,78,395,732),(0,15,30,45)):
    rev={}
    for x in range(1499//a+1):
        rev.setdefault(scalar(a,d,x),[]).append(x)
    for N in range(1500):
        assert sorted(solve_row(a,d,N))==rev.get(N,[])
print('Physical graph counts: 15. Scalar inverses: 36000.')
hits={}
cap=12345
for x,y,j in product(range(1,cap//Q+1),range(1,cap//Q+1),range(3)):
    N=value(x,y,j)
    if N<=cap:hits.setdefault(N,[]).append((x,y,j))
for N in range(1,cap+1,2):
    assert decode(N)==sorted(hits.get(N,[]))
print('Checked all', (cap+1)//2,'odd target inverses through',cap, 'hits',len(hits))
for N in [999,1001,1003,35099,36833,1000000001]:
    print(N,decode(N))
