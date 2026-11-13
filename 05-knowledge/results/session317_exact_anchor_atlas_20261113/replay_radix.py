from math import gcd

def descend(N,n,q,h,D,A,B):
    if not n:
        if not 0<=N<=D:raise ValueError(N)
        return [N]
    Q=q**n
    d=(N*pow(Q,-1,h))%h
    lo=max(0,-((h*B[n-1]-N)//Q))
    hi=min(D,(N-h*A[n-1])//Q)
    d+=max(0,(lo-d+h-1)//h)*h
    if not lo<=d<=hi:raise ValueError((N,n,q,h,d,lo,hi))
    return descend((N-d*Q)//h,n-1,q,h,D,A,B)+[d]

done=0
for q in range(3,12):
    for h in range(2,q):
        if gcd(q,h)!=1:continue
        D=(q+1)*(h-1);A=[0];B=[D]
        for n in range(1,12):
            A.append(h*A[-1]+(h-1)*q**n)
            B.append(h*B[-1]+(D-h+1)*q**n)
        for N in range(12000):
            n=next(i for i in range(len(A)) if A[i]<=N<=B[i])
            digits=descend(N,n,q,h,D,A,B)
            if sum(d*q**i*h**(n-i) for i,d in enumerate(digits))!=N:raise ValueError(digits)
            done+=1
print('literal mixed-power integer decodes:',done)
