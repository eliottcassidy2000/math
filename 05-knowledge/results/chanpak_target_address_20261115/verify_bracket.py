from random import Random
from bracket import bracket, low_signature
from inverse import p
for s in range(7):
    g=1<<s
    for b in range(g):
        A,z=low_signature(g,b)
        for u in range(50):
            assert p(g*u+b)==p(u)+A*u+p(b)+z*(u&1)
print('Exact decomposition checked for 6350 low-digit/source pairs.')
R=Random(231515)
steps=0
worst=(0,0)
for bits in [*range(14,39),48,64,96,128,192,256,384,512]:
    for _ in range(33):
        N=R.randrange(max(10000,1<<bits),1<<(bits+1))
        a,b=bracket(N)
        ga,gb=N-a[2][0],b[2][-1]-N
        assert a[2][1]-a[2][0]==a[2][2]-a[2][1]==2
        assert b[2][1]-b[2][0]==b[2][2]-b[2][1]==2
        assert len({v%3 for v in a[2]})==len({v%3 for v in b[2]})==3
        rat=max(ga,gb)/N**.375
        worst=max(worst,(rat,N))
        steps+=1
print('Verified',steps,'two-sided brackets; worst tested max-error/N^.375:',worst[0])
for N in (35099,10**6+1,10**12+1,10**36+1):
    print(N,bracket(N))
