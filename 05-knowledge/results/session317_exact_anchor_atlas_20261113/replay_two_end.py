from itertools import permutations
def good(row,p): return all(row[i]>>j&1 for i,j in zip(p,p[1:]))
def add_rows(n,arc): return [sum(1<<j for j in range(n) if i!=j and arc(i,j)) for i in range(n)]
def htable(row):
    n=len(row)
    d=[[0]*n for _ in range(1<<n)]
    h=[1]+[0]*((1<<n)-1)
    for i in range(n): d[1<<i][i]=1
    for m in range(1,1<<n):
        h[m]=sum(d[m])
        for i,u in enumerate(d[m]):
            if u:
                a=row[i]&~m
                while a:
                    b=a&-a;a-=b;d[m|b][b.bit_length()-1]+=u
    return h
row=(10,28,25,16,1)
R=htable(row); assert R[-1]==9
pp=list(permutations(range(5)))
assert [sum(good(row,p) and p[0]==b for p in pp) for b in (0,2)]==[1,3]
assert [sum(p[j-1]==b and good(row,p[:j]) and good(row,p[j:]) for p in pp for j in range(1,6)) for b in (0,2)]==[11,11]
for L in range(1,5):
    n=L+7; vals=[]
    for b in (0,2):
        def arc(i,j):
            if i==n-1: return j==0 or L<j<n-1 and j!=L+1+b
            if j==n-1: return not arc(j,i)
            if i>L and j>L: return row[i-L-1]>>(j-L-1)&1
            return i<j
        tab=htable(add_rows(n,arc))
        vals.append((tab[-1],sum(x*y for x,y in zip(tab,reversed(tab)))))
    assert [p[0] for p in vals]==[20*(1<<L)+19,20*(1<<L)+17]
    F2=2**L-2; F3=3**L-3*(2**L)+3
    assert vals[1][1]-vals[0][1]==-128-88*F2-8*F3
print("absolute two-ended gate and two-slot loss through four M-orders")
