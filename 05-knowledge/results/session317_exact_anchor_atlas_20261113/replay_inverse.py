from bisect import bisect_left

def K(B):
    return K(B//2)+(B%2)*(B//2+(B//2)%2) if B else 0

def F(B,h,E,s,c):
    return c+(h+E)*B-h*K(B)-s*(B%2)

def inverse(N,h,E,s,c):
    L=N//E
    m=1
    while E*2**m<h*(L//2**m+2):m+=1
    T=2**m; found=[]
    for Y in range(min(T,L+1)):
        rev=int(f'{Y:0{m}b}'[::-1],2)
        eta=Y//(T//2)
        alpha=(h+E)*T-h*rev
        beta=F(Y,h,E,s,c)
        def f(X):return beta+alpha*X-h*K(X)-h*eta*(X%2)
        a,b=0,(L-Y)//T+1
        top=b
        while a<b:
            mid=(a+b)//2
            if f(mid)<N:a=mid+1
            else:b=mid
        if a<top and f(a)==N:found.append(T*a+Y)
    return sorted(found)

for h,E,s,c in [(9,245,0,91),(3,12,3,13),(1,1,0,1),(33,1,12,33)]:
    bucket={}
    for B in range(12001//E+1):
        bucket.setdefault(F(B,h,E,s,c),[]).append(B)
    for N in range(0,12001,3):
        if inverse(N,h,E,s,c)!=bucket.get(N,[]):
            raise ValueError((h,E,s,c,N))
print('exact target inverse agrees with exhaustive prefix values')
