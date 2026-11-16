from functools import lru_cache

def H(a):
    n=len(a)
    for i in range(n):
        assert a[i]>>i&1==0
        for j in range(i):
            assert (a[i]>>j&1)+(a[j]>>i&1)==1
    @lru_cache(None)
    def path(mask,j):
        m=mask^(1<<j)
        return sum(path(m,i) for i in range(n) if m>>i&1 and a[i]>>j&1) if m else 1
    return sum(path((1<<n)-1,j) for j in range(n))

def step(a,p,mask):
    a=list(a);n=len(a)
    if p is not None:
        i,j=p
        a[i]^=1<<j
        a[j]^=1<<i
    for i in range(n):
        if not mask>>i&1:a[i]|=1<<n
    return tuple(a+[mask])

old=(62,124,120,112,32,64,17)
assert H(old)==27
for h,p,m in [(53,None,16),(105,(0,1),181)]:
    old=step(old,p,m)
    assert H(old)==h
assert old==(188,509,248,496,32,192,401,16,181)
for h,p,m in [(209,(2,3),208),(417,(3,8),656),(833,(4,7),692)]:
    old=step(old,p,m)
    assert H(old)==h
assert old==(3772,4093,1776,3828,160,1728,3473,0,3773,208,2704,692)
late=(6,12,24,49,99,71,15,383,127)
assert H(late)==105
for h,p,m in [(209,(0,1),31),(417,(1,3),2),(833,(1,10),1088),(1665,(3,5),3087)]:
    late=step(late,p,m)
    assert H(late)==h
print('Original linked moves through 833; separate 105 core through 1665.')
