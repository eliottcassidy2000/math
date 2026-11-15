from replay_turnback_focus import fields, slots, turnback


def wrap(r):
    n=len(r)
    return tuple(v|(1<<n) for v in r)+(1<<(n+1), (1<<n)-1)


def strong(r):
    for j in range(len(r)):
        x=1<<j
        for _ in r:
            for i,a in enumerate(r):
                if x>>i&1: x|=a
        if x!=(1<<len(r))-1:
            return False
    return True


def main():
    A=(172,225,2,70,79,92,133,60)
    B=(86,92,152,1,232,143,172,11)
    for r,m,v in [(A,4,-1508),(B,8,-1484)]:
        F=fields(r,m,2)
        h,u,_,K=F[0]; q,u2,v2,J=F[1]
        if (h,q,u,u2,K,J,v2)!=(389,6032,406,3196,0,0,v) or not strong(r):
            raise AssertionError((r,F))
        print('literal old kernel',m,F)
    for j in range(2,5):
        ans=[]
        for r,m in [(A,4),(B,8)]:
            w=turnback(r,m,j,j)
            _,H,Q=slots(w,2)
            got=slots(wrap(w),1)[1][-1]
            if got!=H[-1]+Q[-1]:
                raise AssertionError((j,m,got,H[-1],Q[-1]))
            ans.append((H[-1],Q[-1],got))
        a,b=ans
        if (a[0],b[0],b[1]-a[1],b[2]-a[2])!=(389+406*2**(j-1),389+406*2**(j-1),48,48):
            raise AssertionError((j,ans))
        print('literal returned and wrapped',j,a,b)

if __name__=='__main__':
    main()
