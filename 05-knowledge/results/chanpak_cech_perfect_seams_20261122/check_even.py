from collections import Counter, defaultdict
from math import prod
def strong(rows):
    allv=(1<<len(rows))-1
    for i in range(len(rows)):
        z=1<<i
        while True:
            a=z
            for j in range(len(rows)):
                if z>>j&1:a|=rows[j]
            if a==z:break
            z=a
        if z!=allv:return False
    return True


def cap_two(rows,marks):
    for a in marks:
        E=[(o,i) for o in range(len(rows)) if rows[a]>>o&1
                 for i in range(len(rows)) if (rows[i]>>a&1) and (rows[o]>>i&1)]
        assert all(not ((u!=v or u in marks) and (x!=y or x in marks))
                   for u,x in E for v,y in E)


def signature(rows,marks):
    cap_two(rows,marks)
    n=len(rows)
    U=((1<<n)-1)^sum(1<<j for j in marks)
    dp={(0,(0,)*len(marks),-1):1}
    data=Counter()
    for _ in range(n+len(marks)+1):
        m=defaultdict(int)
        for (state,counts,end),ways in dp.items():
            if state==U and all(counts):
                data[counts]+=ways
            for j in range(n):
                if end>=0 and not rows[end]>>j&1:continue
                if j in marks:
                    i=marks.index(j)
                    if counts[i]==2:continue
                    v=counts[:i]+(counts[i]+1,)+counts[i+1:]
                    m[state,v,j]+=ways
                elif not state>>j&1:
                    m[state|(1<<j),counts,j]+=ways
        dp=m
    return dict(data)


def power_coeffs(sig):
    r=[0]*8
    for s,v in sig.items():
        y=sum(1<<j for j in range(3) if s[j]==2)
        q=y
        while True:
            r[q]+=v*(-2)**((q^y).bit_count())
            if not q:break
            q=(q-1)&y
    return tuple(r)


def H(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for j in range(n):dp[1<<j][j]=1
    for k in range(1,1<<n):
        for j,v in enumerate(dp[k]):
            if not v:continue
            z=rows[j]&~k
            while z:
                b=z&-z;t=b.bit_length()-1
                dp[k|b][t]+=v
                z^=b
    return sum(dp[-1])


def inflate(rows,ports,widths):
    values=dict(zip(ports,widths))
    seq=[i for i in range(len(rows)) for _ in range(values.get(i,1))]
    return tuple(sum(1<<j for j,w in enumerate(seq)
                  if (w!=v and rows[v]>>w&1) or (w==v and i<j))
                 for i,v in enumerate(seq))


def cube(sig,widths):
    return sum(c*prod((1<<w)-2 if k==2 else 1 for w,k in zip(widths,r))
               for r,c in sig.items())


rows=(2,16,11,19,133,31,319,111,191)
marks=(8,5,1)
s=signature(rows,marks)
assert strong(rows)
assert power_coeffs(s)==(-1,-1,-5,-5,5,5,8,8)
for p in (6,8,10,24,62,120):
    w=(p,p//2-2,p//2-2)
    v=2**(p-1)*(2**p-1)-1
    assert cube(s,w)==v
    if p==6: assert H(inflate(rows,marks,w))==v
print(s)
print('PASS')

upper=(2,16,11,19,133,31,63,111)
upper_marks=(1,3,5)
su=signature(upper,upper_marks)
assert strong(upper)
assert tuple(v for _,v in sorted(su.items()))==(53,22,27,12,37,16,18,8)
assert power_coeffs(su)==(1,1,-1,2,-2,0,-4,8)
for p in (6,8,10,24,62,120):
    w=(p//2-1,p-1,p//2-2)
    v=2**(p-1)*(2**p-1)+1
    assert cube(su,w)==v
    if p==6:
        assert H(inflate(upper,upper_marks,w))==v
print('Even upper rail PASS')
