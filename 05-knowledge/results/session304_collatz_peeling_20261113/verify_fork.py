def hp(rows):
    n=len(rows)
    dp=[{} for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i,c in dp[s].items():
            a=rows[i]&~s
            while a:
                b=a&-a; a-=b; j=b.bit_length()-1
                dp[s|b][j]=dp[s|b].get(j,0)+c
    return sum(dp[-1].values())

def word(w):
    n=len(w)
    return tuple(sum(1<<j for j in range(i+1,n))+(int(w[i]=='0')<<n)
                 for i in range(n))+(sum(1<<j for j,e in enumerate(w) if e=='1'),)

def extend(r,loses):
    n=len(r);mask=sum(1<<i for i in loses)
    return tuple(x|(((mask>>i)&1)<<n) for i,x in enumerate(r))+(((1<<n)-1)^mask,)

def crossing(r,i):
    O=[j for j in range(len(r)) if r[i]>>j&1]
    I=set(range(len(r)))-set(O)-{i}
    return [(j,k) for j in O for k in sorted(I) if r[j]>>k&1]

def eq(a,b):
    if a!=b: raise AssertionError((a,b))

eq(crossing(extend(word("101010"),[0,2]),3),[(6,0),(6,2)])
eq(crossing(extend(word("100100"),[1,6]),3),[(4,6),(5,6)])
for k in range(1,9):
    z=2**k
    r=word("101"+"0"*k+"10")
    eq(hp(r),15*z+1)
    eq(hp(word("1001"+"0"*k+"110")),45*z+1)
    eq(hp(extend(r,[0,2])),45*z+5)
    if k>=2:
        loses=[i for i in range(len(r)) if i!=4 and (i==3 or r[i]>>3&1)]
        loses.remove(3); loses.append(4)
        eq(hp(extend(r,loses)),45*z+3)
    r=word("100"+"1"*k+"00")
    eq(hp(r),15*z+3)
    eq(hp(extend(r,[1,len(r)-1])),45*z+7)
print("55 direct subset-DP controls and two crossing signatures verified")
