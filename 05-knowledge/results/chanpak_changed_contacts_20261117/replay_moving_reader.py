from math import factorial, isqrt
from replay_mersenne import need, ordinal, strong, table

SCAFFOLD=(8,1,3,6)
def converse(rows):
    return tuple(sum((row>>i&1)<<j for j,row in enumerate(rows))
                 for i in range(len(rows)))

def inflate(block):
    labels=[(0,i) for i in range(len(block))]+[(i,0) for i in (1,2,3)]
    def arc(i,j):
        x,u=labels[i]; y,v=labels[j]
        return (block[u]>>v&1) if x==y==0 else (SCAFFOLD[x]>>y&1)
    return tuple(sum(arc(i,j)<<j for j in range(len(labels))) for i in range(len(labels))),len(block)

def profile(rows,v):
    d,H=table(rows)
    b=sum(R[v]*H[(len(H)-1)^S] for S,R in enumerate(d) if S>>v&1)
    c=table(converse(rows))[0][-1][v]
    need(d[-1][v]>0)
    return H[-1],b,c

def join(left,v,right_reverse,u):
    n=len(left); mask=(1<<n)-1
    result=list(left)+[(row<<n)|mask for row in converse(right_reverse)]
    result[v]|=1<<(n+u)
    result[n+u]^=1<<v
    return tuple(result)


def transitive(n):
    M=(1<<n)-1
    return tuple(M^((1<<(i+1))-1) for i in range(n))

for p,m,reader,u,expected in [
    (257,1,(34,124,9,17,37,76,29),6,(101,146,20)),
    (65537,4,(34,460,9,177,167,326,157,293,93),2,(1907,5711,139)),
]:
    left,v=inflate(transitive(m))
    a,b,c=profile(left,v)
    P=1<<m
    need((a,b,c)==(3*P-1,4*P+1,P-1))
    d,e,f=profile(reader,u)
    need((d,e,f)==expected)
    child=join(left,v,reader,u)
    need(a*d+b*e-c*f==7*p and strong(child))
    need(table(child)[1][-1]==7*p)
    need(all(p%i for i in range(2,isqrt(p)+1)))
    print(p,len(child),3*d+4*e-f,-d+e+f,child)

# Necessary endpoint profiles left after the elementary small-target screening.
for N,expected in [(35,[]),(119,[(1,2,5,11,5),(1,2,9,9,7)])]:
    survivors=[]
    for n in range(1,8):
        for h in range(1,min(factorial(n),N//(6*(n+1)+1))+1,2):
            for P in range(max(2**n,(n+1)*h),N//6+1):
                if P==2 and (n,h)!=(1,1):
                    continue
                if n==1 and (h,P)!=(1,2):
                    continue
                if n>=2 and P%4:
                    continue
                for d in range(1,N//(6*P+h)+1,2):
                    if d in (7,21):
                        continue
                    if d==1:
                        actual=[(1,1)]+[(1<<j,0) for j in range(1,6)]
                    else:
                        actual=[(e,f) for e in range(d,N//(4*P+h)+d+1) for f in range(d+1)]
                    for e,f in actual:
                        if (3*P-h)*d+(4*P+h)*e-(P-h)*f==N:
                            survivors.append((h,P,d,e,f))
    need(survivors==expected)
    print('necessary profiles',N,survivors)

need(27*factorial(25)-1 < 2**89)
need(pow(2,18,27)==1 and [s for s in range(18) if pow(2,s,27)==26]==[9])
need([s+7 for s in range(89) if pow(2,s,27)==26]==[16,34,52,70,88])
