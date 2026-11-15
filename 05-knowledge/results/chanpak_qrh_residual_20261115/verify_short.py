from short_guard import bracket, certify
from random import Random
R=Random(921015)
c=0
for m,expected in [(7,13),(11,12),(13,13)]:
    n=certify(m)
    assert n==expected
    print("Shortest universal exact-length guard",m,n)
    for k in [100,127,255,511]:
        for i in range(13):
            N=R.randrange(1<<k,1<<(k+1))
            S,a,b=bracket(N,m,n)
            assert a[0]+a[1]==b[0]+b[1]==S
            assert a[-1][a[2]]%(3*m)==b[-1][b[2]]%(3*m)==N%(3*m)
            c+=1
print('Joint fixed-sum full brackets:',c)
