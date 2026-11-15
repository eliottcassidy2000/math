from replay_carry_inverse import dp,p

def apex(r,mask):
    n=len(r)
    return tuple(x|(1<<n if not mask>>i&1 else 0) for i,x in enumerate(r))+(mask,)

def append(rows,mask,w):
    rows=list(rows)
    for x in w:
        n=len(rows)
        for i in range(n):rows[i]|=1<<n
        rows.append(0)
        if x=='1':mask|=1<<n
    return tuple(rows),mask

def state(r, mask):
    n=len(r); f=(1<<n)-1
    a=[[0]*n for _ in range(f+1)]
    hh=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if s>>i&1:
                z=s^(1<<i)
                a[s][i]=sum(a[z][j] for j in range(n) if (r[i]&z)>>j&1) if z else 1
        hh[s]=sum(a[s])
    E=O=0
    for k in range(1,f+1):
        v=sum(a[k][i] for i in range(n) if mask>>i&1)*hh[f^k]
        if k.bit_count()%2:O+=v
        else:E+=v
    return hh[-1],dp(apex(r,mask)),E+O,O-E

def test(arms, ready):
    for (r,mask,tail),(h,N,S,d) in zip(arms,ready):
        before=state(r,mask)
        a,b=append(r,mask,tail)
        got=state(a,b)
        print('initial and aligned',r,mask,tail,before,got)
        if got!=(h,N,S,d):raise AssertionError(got)
        for more in (0,2):
         for t in (0,1,2,3,5,7):
            w=''.join('0' if c=='1' else '1' for c in bin(t)[2:][::-1]) if t else ''
            a,b=append(r,mask,tail+'1'*more+'0'+w)
            nS=(1<<more)*(S+h)-h
            want=N+nS+d+2*nS*t+h*p(t)
            actual=dp(apex(a,b))
            if actual!=want:raise AssertionError((w,actual,want))

def main():
    r=(48,29,25,1,8,30)
    s=(18,20,17,7,8)
    test([(r,25,'1'),(r,16,'11'),(s,1,'11'),(r,1,'11')],[(15,N,121,15) for N in (27,29,31,33)])
    r=(6,8,10,1)
    test([(r,8,'11'),(r,6,'1')],[(5,N,35,5) for N in (9,11)])
    print('72 literal matched-window whole-glue computations')
if __name__=='__main__':main()
