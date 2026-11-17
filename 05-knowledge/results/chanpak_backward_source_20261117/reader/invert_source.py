from replay import U,rev,pair,ck
def invert_source(N,j,A,r):
    q=1<<j
    if j<1 or A<q or A%2 or not 0<=r<q:return None
    rho=rev(r,j);G=A+rho
    if rho%2 or G%4!=2:return None
    L=2*q*A*A+r*G+q*(1-U(A)-U(r))
    if L<=0:return None
    h,z=divmod(q*N-4*q*A-2*r,L)
    if z or h<1 or not (h%2):return None
    p,z=divmod(h*G+2,q)
    if z or p<2*h or p%2:return None
    return h,p
def main():
    ck(invert_source(158860293255395167777,25,634979302,0),(197,3728))
    cases=0
    for h,p in ((1,2),(1,4),(5,32),(197,3728)):
        for j in range(1,35):
            q=1<<j
            for rho in (0,min(2,q-2),q-2):
                v=pair(h,p,j,rho)
                if not v:continue
                N,u=v;A,r=divmod(u,2*q)
                ck(invert_source(N,j,A,r),(h,p))
                cases+=1
    print(cases,'source-unaware inverse equalities, plus explicit crossover.')
if __name__=='__main__':main()
