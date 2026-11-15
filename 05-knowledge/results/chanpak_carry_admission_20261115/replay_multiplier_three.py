from replay_matched_windows import state, append, apex
from replay_carry_inverse import dp
from itertools import product

def trans(st,w,h):
    N,S,d=st
    for c in w:
        b=int(c)
        N,S,d=N+(1-b)*(S+d),2*S+b*h,b*h
    return N,S,d

def direct(st,w,pos,h):
    u,v=w[:pos],w[pos+2:]
    i=int(w[pos:pos+2],2)
    N,S,d=trans(st,u,h)
    mid=3*N+2*(i in (1,2))*h+(9-3*i)*S+(3-i)*d,12*S+3*h*i,h*i
    return trans(mid,v,3*h)

def main():
    n=0
    for h,N,S in [(5,49,70),(15,163,242),(1,3,2),(7,101,7)]:
        for l in range(2,10):
            for chars in product('01',repeat=l):
                w=''.join(chars)
                NN,SS,_=trans((N,S,0),w,h)
                gaps=set()
                for pos in range(l-1):
                    i=int(w[pos:pos+2],2)
                    _,s,d=trans((N,S,0),w[:pos],h)
                    e=w[pos+2:pos+3]=='0'
                    J=0 if i in (0,3) else 3*s-d+2*h-2*h*e if i==1 else -3*s+d-4*h+2*h*e
                    x,y,_=direct((N,S,0),w,pos,h)
                    if (x,y)!=(3*NN+J,3*SS):raise AssertionError((w,pos,x,y))
                    gaps.add(J);n+=1
                a=sorted(gaps)
                if any(y-x<2*h for x,y in zip(a,a[1:])):raise AssertionError((w,a))
    print('Effective-word position and offset checks:',n)
    tri=(2,4,1)
    for w,pos in [('00101',1),('11010',2)]:
        R,mask=append((6,8,10,1),6,'10'+w[:pos])
        n=len(R)
        m=int(w[pos:pos+2],2)
        R=tuple(r|(7<<n) for r in R)+tuple(r<<n for r in tri)
        mask|=((1<<m)-1)<<n
        R,mask=append(R,mask,w[pos+2:])
        actual=dp(apex(R,mask))
        exp=direct((51,70,0),w,pos,5)[0]
        if actual!=exp:raise AssertionError((w,pos,actual,exp))
    print('Two independent literal marked-tournament multiplier-three checks')
if __name__=='__main__':main()
