from itertools import product
from replay_intercyclic import literal
def word(j,a,b):
    return ('1' if j==0 else '2')+'0'*a+'30'+'2'*(b+1)+('1' if j==0 else '12')
def terms(j,b):
    t,r=2**b,3**b
    return (1+18*t+12*r,-6+18*t+42*r) if j==0 else (5+36*t+24*r,8+48*t+36*r)
def solve(j,N):
    b=1
    while 3**b<N:
        c,d=terms(j,b)
        A,q=divmod(N-c,d)
        if q==0 and A>=2 and A&(A-1)==0:yield A.bit_length()-1,b
        b+=1
for j in (0,1):
    grid=[]
    for a,b in product((1,2,3),repeat=2):
        c,d=terms(j,b);N=c+d*2**a
        assert literal(word(j,a,b))==N
        assert (a,b) in list(solve(j,N))
        grid.append(N)
    print('family',j,'all-length coefficient grid',grid)
