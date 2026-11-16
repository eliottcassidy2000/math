from functools import cache
FIVE=((2,260,273,7,43,15,319,383,57),(0,257,267,35,271,87,415,319,41),
      (0,257,323,23,39,271,315,383,25),(272,5,9,259,14,31,319,383,54))
def join(a,b,root):
    n=len(a);m=len(b);ar=((1<<m)-1)<<n
    return tuple(v|(0 if i==root else ar) for i,v in enumerate(a))+tuple((v<<n)|(1<<root) for v in b)
def H(a):
    n=len(a);full=(1<<n)-1
    @cache
    def f(M,i):
        old=M^(1<<i)
        return sum(f(old,j) for j in range(n) if (old&a[i])>>j&1) if old else 1
    z=sum(f(full,i) for i in range(n))
    f.cache_clear()
    return z
for i,a in enumerate(FIVE):
    b=join(a,(0,),8);c=join(b,(2,4,1),8)
    assert H(b)==259+2*i and H(c)==3783+6*i
# One backward arc with exactly four vertices between its endpoints.
u=[((1<<6)-1)^((1<<(i+1))-1) for i in range(6)]
u[0]^=1<<5;u[5]^=1
c=tuple(x|(1<<6) for x in u)+(1<<7,63)
assert H(tuple(u))==17 and H(c)==241
print('Literal rooted triangle amplifiers and 17→241 checked.')
