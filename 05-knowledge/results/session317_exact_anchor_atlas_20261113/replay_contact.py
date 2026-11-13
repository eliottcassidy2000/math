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
def require(ok):
    if not ok:
        raise RuntimeError("contact replay mismatch")

def val(n,p=3):
    e=0
    while n%p==0: e+=1;n//=p
    return e

def nine(i,j):
    return (j//3-i//3)%3==1 if i//3!=j//3 else (j-i)%3==1

def hq(rows):
    tab=htable(rows)
    return tab[-1],sum(x*y for x,y in zip(tab,reversed(tab)))

h,q=hq(add_rows(9,nine)); S=(q-h)//9
require((h,q,S)==(3159,48924,5085))
for r in range(2,40):
    f=S*h**(r-1)+h**(r-2)*(h-h//9)*(h-(r==2)*(h//9))+S*(S+h)*q**(r-2)
    require(val(f)==4*r-4)
print("ternary nine-block",h,q,S)

row=(14,28,24,16,1)
require(hq(row)==(9,84))
pp=list(permutations(range(5)))
require(sum(good(row,p) and p[0]==1 for p in pp)==4)
require(sum(good(row,p) and p[-1]==1 for p in pp)==1)
require(sum(p[j]==1 and good(row,p[:j]) and good(row,p[j:]) for p in pp for j in range(5))==25)
require(sum(p[j-1]==1 and good(row,p[:j]) and good(row,p[j:]) for p in pp for j in range(1,6))==10)
def rigid(r):
    return 10*9**(r-1)+5*(9-(r==2))*9**(r-2)+475*84**(r-2)
for r in range(2,40): require(val(rigid(r))==r-2)
for r in (2,3):
    n=5*r+1
    win={1,*[j for j in range(5*(r-1),5*r) if j%5!=1]}
    def arc(i,j):
        if i==n-1: return j in win
        if j==n-1: return i not in win
        if i//5==j//5: return row[i%5]>>(j%5)&1
        return i<j
    h=htable(add_rows(n,arc))[-1]
    require(h==rigid(r))
    print("rigid literal aperture",r,h)
