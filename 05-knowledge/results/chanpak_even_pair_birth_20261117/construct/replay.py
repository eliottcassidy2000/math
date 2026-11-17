from math import isqrt

def ck(a,b):
    if a!=b:raise AssertionError((a,b))
def TT(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def table(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*(len(d)-1)
    for i in range(n):d[1<<i][i]=1
    for m in range(1,len(d)):
        h[m]=sum(d[m])
        for i,c in enumerate(d[m]):
            if not c:continue
            z=a[i]&~m
            while z:
                b=z&-z;z-=b;d[m|b][b.bit_length()-1]+=c
    return h
def profile(a):
    h=table(a);f=len(h)-1
    return h[-1],sum(h[m]*h[f^m] for m in range(f+1))
def strong(a):
    for i in range(len(a)):
        s=1<<i
        for _ in range(len(a)):
            for j,r in enumerate(a):
                if s>>j&1:s|=r
        if s!=(1<<len(a))-1:return False
    return True
def ext(R,mask):
    m=len(R)
    return tuple(a|(((mask>>i&1)^1)<<m)|(1<<(m+1)) for i,a in enumerate(R))+(mask,1<<m)
def subst(a,slot,B):
    d=len(B)-1;n=len(a)
    def old(i):return i+(d if i>slot else 0)
    b=[]
    for i in range(n):
        c=sum(1<<old(j) for j in range(n) if j!=slot and a[i]>>j&1)
        if i==slot:b.extend(c|(r<<slot) for r in B)
        else:b.append(c|((((1<<len(B))-1)<<slot) if a[i]>>slot&1 else 0))
    return tuple(b)
R=(
((118,56,18,37,8,20,62),52,18,13),
((16,181,33,71,44,9,55,125),12,67,9),
((350,68,16,6,42,207,28,95,254),14,109,7),
)
T=(782,780,776,768,751,975,783,847,16,256)
S=(782,652,776,784,679,975,543,845,82,256)
def prime(n):
    ck(all(n%d for d in range(2,isqrt(n)+1)),True)
def main():
    for a,q,alpha,beta in R:
        slot=len(a)-1;z=ext(a,q)
        ck(strong(a) and strong(z),True)
        ck(sum(r>>slot&1 for r in a),1)
        ck(sum(r>>slot&1 for r in z),1)
        for b in (TT(1),TT(2),TT(3),(2,4,1)):
            h,p=profile(b)
            U=subst(a,slot,b);V=subst(z,slot,b)
            ck(V[:len(U)],tuple((r&((1<<len(U))-1))|(v&~((1<<len(U))-1)) for r,v in zip(U,V)))
            ck((table(U)[-1],table(V)[-1]),((alpha*p+beta*h),7*(alpha*p+beta*h)))
        print('induced',alpha,beta)
    ck(T[-1],S[-1]);ck(T[-1].bit_count(),1)
    ck(strong(T) and strong(S),True)
    for b in (TT(1),TT(2),TT(3),(2,4,1)):
        h,p=profile(b);x=128*p+h
        ck((table(subst(T,9,b))[-1],table(subst(S,9,b))[-1]),(x,7*x))
    b=(90,40,59,64,42,9,54);ck(strong(b),True)
    h,P=profile(b);ck((h,P),(71,1012))
    p,q=872*P+7*h,128*P+h
    ck((p,q,7*p,7*q),(882961,129607,6180727,907249))
    ck(((7*q-p)//24,(16*p-109*q)//3),(P,h))
    prime(p);prime(q);print('prime pair',p,q,'sevenfold values',7*p,7*q)
if __name__=='__main__':main()
