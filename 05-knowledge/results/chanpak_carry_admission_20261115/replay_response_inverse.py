from replay_carry_inverse import p, sig

def v2(a):
    return (a&-a).bit_length()-1

bank={1:(1,),3:(12,),9:(84,144)}
factors={1:((1,1),),3:((1,3),(3,1)),9:((1,9),(3,3),(9,1))}

def inverse_d(L,t,k):
    q=1<<k
    target=(pow(q,-1,t)*L)%t if t>1 else 0
    d=(L//q)-(L//q-target)%t
    z,rem=divmod(L-d*q,t)
    return (d,q-1-z) if d>0 and d%2==0 and not rem and 0<=z<q else None

def recover(L,t):
    out=set()
    for k in range(L.bit_length()):
        rec=inverse_d(L,t,k)
        if rec is None:continue
        d,r=rec
        v=f'{r:0{k}b}'[::-1] if k else ''
        for u,h in factors[t]:
            for Qu in bank[u]:
                for Qw in bank[h]:
                    s=v2(d)-v2(Qw)
                    if s<0 or h==1 and s==0:continue
                    qw=Qw<<s
                    if d%qw:continue
                    qu=d//qw+u
                    z,rem=divmod(qu,Qu)
                    if rem or z<=0 or z&(z-1) or u==1 and z==1:continue
                    A,c,t_=sig((u,qu),(h,qw),v)
                    if A-t_!=L or t_!=t:raise AssertionError((L,t,d,r))
                    out.add(c)
    return out

def direct(limit,t):
    ans={}
    for u,h in factors[t]:
        for Qu in bank[u]:
            for Qw in bank[h]:
                for s in range(u==1,limit.bit_length()+1):
                    for j in range(h==1,limit.bit_length()+1):
                        qu,qw=Qu<<s,Qw<<j
                        d=(qu-u)*qw
                        if d>limit:continue
                        for k in range((limit//d).bit_length()):
                            for b in range(1<<k):
                                v=f'{b:0{k}b}' if k else ''
                                D,c,t_=sig((u,qu),(h,qw),v)
                                L=D-t_
                                if L<=limit:ans.setdefault(L,set()).add(c)
    return ans

def main():
    for t in range(1,20,2):
        for L in range(1,257):
            for k in range(6):
                q=1<<k
                brute=[(d,q-1-(L-d*q)//t) for d in range(2,L//q+1,2) if (L-d*q)%t==0 and (L-d*q)//t<q]
                got=inverse_d(L,t,k)
                if brute!=([got] if got else []):raise AssertionError((L,t,k))
    print('15360 width/congruence inverse controls')
    for t in factors:
        exp=direct(1000,t)
        for L in range(1,1001):
            if recover(L,t)!=exp.get(L,set()):raise AssertionError((L,t,recover(L,t),exp.get(L)))
    print('3000 complete fixed-response lists compared with direct finite-bank forward controls')
    got=recover(252,9)
    if got!={207,243}:raise AssertionError(got)
    print('Exact full constant list at L=252,t=9:',sorted(got))
if __name__=='__main__':main()
