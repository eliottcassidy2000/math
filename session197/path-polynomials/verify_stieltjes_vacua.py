#!/usr/bin/env python3
from verify import *

def ab(par,n):return n*n+6*n+(3 if par else 5),2*n*(n if par else n+1)*(2*n+1)
def moments(par,N):
    p=mon(1 if par else 2);out=[]
    for _ in range(N+1):out.append(value(p,1));p=D(p)
    return out

def run():
    reports=[]
    for par in (1,0):
        ell=[Fraction(ab(par,0)[0])]
        for n in range(1,101):
            a,b=ab(par,n);ell.append(Fraction(a)-b/ell[-1])
            require(ell[-1]>0,'positive pivot')
            if par:require(ell[-1]>=(n+1)**2,'odd induction bound')
            if not par and n>=8:require(ell[-1]>=n*n+2*n-4,'even induction bound')
        if not par:
            require(ell[:9]==list(map(Fraction,[5,'48/5','59/4','1216/59','4185/152','10052/279','16762/359','503016/8381','1597341/20959'])),'exact even boundary')
            for n in range(9,101):
                a,b=ab(par,n);f=n*n+2*n-4;prev=n*n-5
                require((a-f)*prev-b==(3*n+5)*(n-9),'even induction identity')
        mu=moments(par,20);alpha=[Fraction(0),ell[0]]
        for n in range(1,21):alpha.extend([Fraction(ab(par,n)[1])/ell[n-1],ell[n]])
        # Literal weighted Dyck walks: up weight1, down from height h weight alpha_h.
        dyck=[Fraction(1)]+[Fraction(0)]*40;ds=[1]
        for step in range(1,41):
            nxt=[Fraction(0)]*41
            for h,v in enumerate(dyck):
                if h<40:nxt[h+1]+=v
                if h:nxt[h-1]+=v*alpha[h]
            dyck=nxt
            if step%2==0:ds.append(dyck[0])
        require(ds==mu,'Stieltjes Dyck moment identity')
        shifted=[]
        for N in range(1,8):
            product0=1
            for k in range(1,N):product0*=ab(par,k)[1]**(N-k)
            prodell=Fraction(1)
            for x in ell[:N]:prodell*=x
            value1=det([[mu[i+j+1] for j in range(N)] for i in range(N)])
            require(value1==product0*prodell,'shifted determinant product')
            for s in range(5):require(det([[mu[i+j+s] for j in range(N)] for i in range(N)])>0,'shifted positive minor')
            shifted.append(value1)
        nonprincipal=0
        for rows in combinations(range(6),3):
            for cols in combinations(range(6),3):
                require(det([[mu[i+j] for j in cols] for i in rows])>0,'strict nonprincipal minor');nonprincipal+=1
        a=(0,) if par else transitive(2);physical=[]
        for n in range(5):
            actual=h_subsets(a)[-1];require(actual==mu[n],'physical vacuum count');physical.append(actual);a=wrapper(a)
        reports.append({'parity':'odd' if par else 'even','first_pivots':[str(x) for x in ell[:9]],'first_alphas':[str(x) for x in alpha[1:10]],'moments':mu,'shifted_leading_determinants':shifted,'nonprincipal_3x3_minors':nonprincipal,'literal_physical_moments':physical})
    hostiles=[]
    for n in (4,5):
        a=transitive(n);hs=[]
        for k in range(3):hs.append(h_subsets(a)[-1]);a=wrapper(a)
        hdet=hs[0]*hs[2]-hs[1]**2;require(hdet==(-48 if n==4 else -444),'nonvacuum hostile');hostiles.append({'start_order':n,'counts':hs,'determinant':hdet})
    print(json.dumps({'vacua':reports,'nonvacuum_hostiles':hostiles},sort_keys=True,indent=2))
if __name__=='__main__':run()
