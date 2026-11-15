# Self-contained direct endpoint-DP check of the physical witnesses.
def dp(rows, h_all=False):
    n=len(rows);h=[1]+[0]*((1<<n)-1);a=[[0]*n for i in h]
    for m in range(1,len(h)):
        for i in range(n):
            if (m>>i)&1:
                s=m^(1<<i)
                a[m][i]=1 if not s else sum(a[s][j] for j in range(n) if s&rows[i]&(1<<j))
        h[m]=sum(a[m])
    return (h,a) if h_all else h[-1]

def apex(r,mask):
    return tuple(x+(1<<len(r) if not(mask>>i&1) else 0) for i,x in enumerate(r))+(mask,)

def append(r,mask,w):
    r=list(r)
    for v in w:
        k=len(r)
        r=[x|(1<<k) for x in r]+[0]
        if v=='1':mask|=1<<k
    return r,mask

def state(r,m):
    H,a=dp(r,True);f=len(H)-1;p=q=0
    for s in range(1,f+1):
        x=H[f^s]*sum(a[s][i] for i in range(len(r)) if m>>i&1)
        if s.bit_count()%2:p+=x
        else:q+=x
    return H[-1],dp(apex(r,m)),p+q,p-q

def p(t):
    if not t:return 0
    u,b=divmod(t,2)
    return p(u)+(0 if b else u+(u%2))

def tail(t):
    return '' if not t else ''.join('1' if c=='0' else '0' for c in bin(t)[2:][::-1])

P=[((2,12,49,21,67,27,47,127),232),((2,36,25,67,11,29,55,127),228)]
T=[((4,5,56,195,75,27,167,55),40),((12,33,50,22,35,201,31,95),167),((2,12,33,85,135,155,55,79),136)]
W=[((10,20,33,38,13,19,63,127),9),((10,20,33,22,37,11,63,127),9),((16,1,3,7,78,31,175,63),18),((10,20,33,38,13,19,63,127),10),((10,20,33,22,37,11,63,127),10),((10,20,33,22,37,11,63,127),12)]
V=[((0,1,43,83,135,155,55,79),39),
   ((0,1,43,83,135,91,151,47),39),
   ((16,1,3,7,78,31,175,63),18),
   ((0,1,43,83,135,155,55,79),43),
   ((0,1,43,83,135,91,151,47),43),
   ((0,1,43,83,135,91,151,47),51)]
def strong(rows):
    for i in range(len(rows)):
        q=1<<i
        for _ in rows:
            for j,r in enumerate(rows):
                if q>>j&1:q|=r
        if q!=(1<<len(rows))-1:return False
    return True

def main():
    for (r,m),N in zip(P,(1267,1269)):
        if state(r,m)!=(81,N,1717,23):raise AssertionError(state(r,m))
        for t in (0,1,3,5):
            H=dp(apex(*append(r,m,'0'+tail(t))))
            if H!=N+1740+3434*t+81*p(t):raise AssertionError(H)
    for (r,m),N,K in zip(T,(1261,1263,1265),(1395,2790,1395)):
        h,val,S,d=state(r,m)
        if (h,val,S+h)!=(239,N,K):raise AssertionError((h,val,S,d))
        for g in (2,3):
            for t in (0,1,3):
                H=dp(apex(*append(r,m,'1'*(g-(K==2790))+'0'+tail(t))))
                if H!=N+(1395<<g)*(2*t+1)-478*t+239*p(t):raise AssertionError(H)
    for (r,m),N in zip(W,range(129,140,2)):
        h,val,S,d=state(r,m)
        if (h,val,S)!=(45,N,129):raise AssertionError((h,val,S,d))
        for t in (0,1,3,5):
            H=dp(apex(*append(r,m,'10'+tail(t))))
            if H!=N+348+606*t+45*p(t):raise AssertionError(H)
    print('Matched pair, aligned triple and aligned six all checked by direct path DP (50 literal extensions)')
    for (r,m),N,K in zip(V,range(129,140,2),(696,696,174,696,696,696)):
        h,val,S,d=state(r,m)
        if (h,val,S+h)!=(45,N,K):raise AssertionError((h,val,S,d))
        for t in (0,1):
            ww='1'*(1 if K==696 else 3)+'0'+tail(t)
            T0=apex(*append(r,m,ww))
            H=dp(T0)
            if not strong(T0) or H!=N+1392*(2*t+1)-90*t+45*p(t):raise AssertionError(H)
    print('Independent literal strength and count check for six strong arms (12 more extensions)')
if __name__=='__main__':main()
