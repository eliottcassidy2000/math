"""New literal one-return certificates and exact triangular target checks."""
from check_leech_labels import paths, bits, reaches, assemble, strong, tt
DYADIC = {-15: ((18, 12, 145, 37, 170, 71, 31, 107), 28, 411, 2048, 426), -11: ((4, 33, 18, 7, 11, 29), 44, 31, 128, 42), -7: ((64, 69, 49, 7, 11, 27, 60), 70, 121, 512, 128), -5: ((24, 5, 193, 22, 166, 15, 59, 107), 21, 375, 2048, 380), -3: ((4, 1, 82, 7, 11, 31, 59), 40, 31, 128, 34), 3: ((2, 24, 3, 5, 13), 20, 11, 32, 8), 5: ((0, 17, 3, 7, 13), 1, 5, 32, 0), 7: ((8, 17, 3, 6, 13), 24, 13, 32, 6), 9: ((16, 1, 3, 7, 14), 4, 9, 16, 0), 13: ((16, 1, 67, 7, 14, 31, 59), 48, 41, 128, 28)}
COXETER = {1: ((4, 1, 82, 7, 11, 31, 59), 40, 31, 128, 34), 3: ((8, 65, 3, 70, 15, 31, 437, 63, 191), 288, 227, 768, 230), 5: ((32, 129, 67, 23, 391, 286, 59, 109, 207), 208, 1193, 5120, 1196), 7: ((32, 65, 3, 7, 15, 30, 61), 38, 85, 448, 88), 9: ((4, 1, 2), 7, 3, 9, 6), 13: ((16, 1, 3, 71, 78, 31, 423, 63, 191), 384, 155, 416, 158), 25: ((8, 17, 3, 70, 13, 31, 55), 36, 49, 200, 52)}

def require(b, context):
    if not b:
        raise AssertionError(context)


def profile(data, beta):
    D,I,h,G,J=data
    U,end,start,H=paths(D)
    g=j=0
    for S in range(1,U+1):
        T=U^S
        e=sum(end(S,q) for q in bits(S & I))
        g+=e*H(T)
        j+=e*sum(start(T,q) for q in bits(T & I))
    require((H(U),g,j)==(h,G,J) and h-J==beta and reaches(D,I), data)
    host=assemble(D,I,tt(2))
    W,_,_,Z=paths(host)
    require(strong(host) and Z(W)==4*G+beta, data)


def factor_offsets(p, r):
    t=1<<(p-1)
    return ((t-(r+1)//2,2*t+r) if r%4==1
            else (t+r//2,2*t-(r+1)))


def check():
    for beta,data in DYADIC.items():
        require(data[3] & (data[3]-1)==0, data)
        profile(data,beta)
    for a,data in COXETER.items():
        g=data[3]
        require(g//(g&-g)==a, data)
        profile(data,-3)
    for p in [13,17,19,31,61]:
        t=1<<(p-1)
        E=t*(2*t-1)
        for r in [2,5,6,9,10,13,14]:
            x,y=factor_offsets(p,r)
            require(x*y==E-r*(r+1)//2 and x>0 and y>0,(p,r))
            for z,pow2 in ((x,t),(y,2*t)):
                beta=z-pow2
                if beta!=1:
                    g=DYADIC[beta][3]
                    require(pow2%g==0 and pow2//g>=2,(p,r,beta))
    for a in [2,3,4,5,7,9,13,25]:
        odd=a//(a&-a)
        g=COXETER[odd][3]
        for k in (10,18,30,60):
            A=a<<k
            require((2*A)%g==0 and (2*A)//g>=2,(a,k))
            require((A+1)*(2*A-3)==2*A*A-A-3,(a,k))
    E=(1<<30)*((1<<31)-1)
    require(factor_offsets(31,2)==(1073741825,2147483645),E)
    print("17 profile checks (16 distinct cores); fixed + signed factors give")
    print("E_31-3 =",E-3,"from",factor_offsets(31,2))


if __name__=='__main__':
    check()

