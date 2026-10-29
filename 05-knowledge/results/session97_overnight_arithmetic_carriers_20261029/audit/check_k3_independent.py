"""Independent stdlib reconstruction of three good reductions at p=233."""
import json
P=233

def need(b,m):
    if not b:raise ArithmeticError(m)
def trim(a):
    while a and a[-1]==0:a.pop()
    return a

def add(*args):
    out=[0]*max(map(len,args))
    for a in args:
        for i,v in enumerate(a):out[i]=(out[i]+v)%P
    return trim(out)
def scale(a,k):return trim([v*k%P for v in a])
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,u in enumerate(a):
        for j,v in enumerate(b):out[i+j]=(out[i+j]+u*v)%P
    return trim(out)
def power(a,n):
    r=[1]
    for _ in range(n):r=mul(r,a)
    return r

def mod(a,b):
    a=a[:]
    while len(a)>=len(b):
        k=len(a)-len(b);v=a[-1]*pow(b[-1],-1,P)%P
        for i,c in enumerate(b):a[k+i]=(a[k+i]-v*c)%P
        trim(a)
    return a

def evaluate(a,t):
    r=0
    for v in reversed(a):r=(r*t+v)%P
    return r

def raw(i):
    if i==273:
        L=[-4436170,48716518,-10340482];p=[-43997,6071604,1345905];q=scale([10560,169808,1757],169)
        D=[-5344262336985,1291548172478,94651330477309,100399938810750,11340016558200]
        E=scale([5489594880,37829024608,-72530579204,59220491056,6409666485],5408)
    elif i==302:
        L=[-239031,-471466,-446667];p=[-72570,368554,318552];q=[-14960,-45082,733413]
        D=scale([-205134150,3909144679,-22069002960,-7114589515,7174492962],5)
        E=[-2386325360,-32493137198,-1174040743,811447034567,882769396002]
    else:
        L=[32416,-828066,466548];p=[103008,-1076788,-478380];q=[81622,-180205,82612]
        D=[-53374213248,239872553168,11566653808,-91335649060,-26666192028]
        E=[-24372268764,76235905664,-104259187849,118525002574,-47423161625]
    L,p,q,D,E=[scale(a,1) for a in [L,p,q,D,E]]
    B=add(mul(p,E),mul(q,D),mul(mul(p,q),add(L,scale(p,-1),scale(q,-1))))
    return L,add(D,E),B,mul(D,E),[]

def main():
    freq=[]
    for b in range(P):
        f=[0]*P
        for y in range(P):f[(y*y+b*y)%P]+=1
        freq.append(f)
    ans={}
    for i in [273,302,4239]:
        aa=raw(i);a1,a2,a3,a4,a6=aa
        b2=add(power(a1,2),scale(a2,4));b4=add(mul(a1,a3),scale(a4,2));b6=add(power(a3,2),scale(a6,4))
        b8=add(mul(power(a1,2),a6),scale(mul(a2,a6),4),scale(mul(mul(a1,a3),a4),-1),mul(a2,power(a3,2)),scale(power(a4,2),-1))
        delta=add(scale(mul(power(b2,2),b8),-1),scale(power(b4,3),-8),scale(power(b6,2),-27),scale(mul(mul(b2,b4),b6),9))
        need(len(delta)==25,'degree24')
        u=delta[:];v=[j*delta[j]%P for j in range(1,len(delta))]
        while v:u,v=v,mod(u,v)
        need(len(u)==1,'squarefree discriminant')
        total=0
        for t in range(P+1):
            vals=[evaluate(a,t) if t<P else (a[w] if len(a)>w else 0) for a,w in zip(aa,[2,4,6,8,12])]
            A1,A2,A3,A4,A6=vals;total+=1
            for x in range(P):total+=freq[(A1*x+A3)%P][(x**3+A2*x*x+A4*x+A6)%P]
        ans[i]=total
    need(ans=={273:58509,302:58404,4239:58374},'three distinct surface counts')
    print(json.dumps({'p':P,'squarefree_degree24':True,'smooth_proper_surface_counts':ans},sort_keys=True))
if __name__=='__main__':main()
