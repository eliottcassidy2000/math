from binary_cycle import create,literal

def delete(A,k):
    return tuple((r&((1<<k)-1))|((r>>(k+1))<<k) for i,r in enumerate(A) if i!=k)

def calc(K,n):
    T,_=create(K,False,spine_bits=n)
    return literal(T)[-1]-literal(delete(T,n))[-1]

if __name__=='__main__':
    for K in range(1,31):
        n=K.bit_length()
        digits=tuple(map(int,f'{K:0{n}b}'))
        X=sum(d*(2**j)*(3**(n-j-1)) for j,d in enumerate(digits))
        Y=sum(d*(2**j) for j,d in enumerate(digits))
        assert calc(K,n)==2*K+2*X+Y-digits[0]
        a,b=calc(K,n+1),calc(K,n+2)
        assert a==2*K+2*(2*X+Y) and 2*a-b==2*K
    print('literal path gap and four-count extraction for 30 cycle targets')

def full_formula(K,n):
    d=tuple(map(int,f'{K:0{n}b}'))
    Y=sum(b*(1<<i) for i,b in enumerate(d))
    B=sum(b*(1<<i)*3**(n-i-1) for i,b in enumerate(d))
    p0=p2=0
    F2=1<<(n+1)
    for i in range(n):
        for j in range(i+1,n):
            if d[i] and not d[j]:
                m=j-i-1
                p0+=m==0
                p2+=1<<m
                F2+=(1<<(n-m-1))*(3**m+1)
    return 1+p0+p2+F2+3*K+2*B+2*Y+d[-1]-(1<<n)-2*d[0]

if __name__=='__main__':
    for K in range(1,30):
        for n in (K.bit_length(),K.bit_length()+1):
            T,_=create(K,False,spine_bits=n)
            assert full_formula(K,n)==literal(T)[-1]
    for K,A,B,res in ((11,37,152,0),(13,45,176,2),(26,97,564,1)):
        for j in range(1,6):
            n=K.bit_length()+j
            h=full_formula(K,n)
            assert h==A+B*(1<<j) and (h%7==0)==(j%3==res)
    assert full_formula(13,5)==397 and full_formula(13,6)==749
    print('whole ordinary spectrum of binary cycle words and three exact mod-seven padding phases')
