from replay_back_cross import ordinal,tt,rev,state

def main():
    T=(18,20,17,7,8);R=(48,29,25,1,8,30)
    U=rev(rev(ordinal(ordinal(T,(0,)),(0,)),0,6),4,6)
    V=rev(rev(rev(ordinal(R,(0,)),0,6),1,6),4,6)
    if (state(U,97),state(V,1))!=((77,173,285,51),(77,175,104,6)):raise AssertionError()
    for j in (3,4,5):
        for i,(r,k) in enumerate(((T,j-2),(R,j-1))):
            n=len(r)
            # First j-1 ones, the guard mark 0, then the original first fold vertex.
            z=ordinal(r,tt(j+1))
            x=n+j
            m=1|(((1<<(j-1))-1)<<n)|((1-i)<<x)
            dest=([0,4]+list(range(n+1,x))) if i==0 else ([0,1,4]+list(range(n,x)))
            for a in dest:z=rev(z,a,x)
            expected=(77,173+2*i+181*(1<<(j-1)),181*(1<<j)-154,0)
            a,b=(U,97) if i==0 else (V,1)
            q=ordinal(a,tt(k+1));s=b|(((1<<k)-1)<<7)
            if state(q,s)!=expected or state(z,m)!=expected:raise AssertionError((i,j,state(z,m),expected))
        print('Literal opposite-mark return at j=',j,':',173+181*(1<<(j-1)),175+181*(1<<(j-1)),'shared h=77 S=',181*(1<<j)-154)
if __name__=='__main__':main()
