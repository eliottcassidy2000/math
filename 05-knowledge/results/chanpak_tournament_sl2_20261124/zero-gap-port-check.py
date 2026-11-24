from itertools import product
def rows(gaps):
    n=1+len(gaps)+sum(gaps)
    a=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
    i=0
    for g in gaps:
        j=i+1+g;a[i]^=1<<j;a[j]^=1<<i;i=j
    return a
def hp(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i,v in enumerate(D[s]):
            if not v:continue
            w=a[i]&~s
            while w:
                b=w&-w;w-=b;D[s|b][b.bit_length()-1]+=v
    return sum(D[-1])
def state(gaps):
    A,O,E,P=1,0,0,0
    for g in gaps:
        if not g:A,O,E,P=A+2*E+P,A+E,O,P
        else:A,O,E,P=A+2*E+P,0,0,(1<<g)*(A+O+E+P)
    return A+2*E+P
for m in range(7):
 for g in product(range(3),repeat=m):
    if 1+len(g)+sum(g)<=10:assert hp(rows(g))==state(g),(g,hp(rows(g)),state(g))
def K(t,a,b):
 z=0
 for x in product(range(2),repeat=t):
    s=(a,)+x+(b,);v=1;j=0
    while j<len(s):
      if not s[j]:j+=1;continue
      i=j
      while j<len(s) and s[j]:j+=1
      if i and j<len(s):v*=2 if (j-i)%2==0 else 0
    z+=v
 return z
assert [[K(3,a,b) for b in (0,1)] for a in (0,1)]==[[5,6],[6,7]]
print('verified exact strip height and unimodular bare port')
