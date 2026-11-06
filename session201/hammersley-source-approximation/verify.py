"""Exact rational controls for approximate physical sources, not target occupancy."""
from fractions import Fraction as F
from functools import lru_cache
import json

def need(c,m):
    if not c:raise RuntimeError(m)

@lru_cache(None)
def g(n):
    need(n>=1,'positive index')
    if n==1:return 1
    if n%2==0:return g(n//2)
    h=(n+1)//2;return g(h)+h-h%2

def rev(k,A):
    r=0
    for _ in range(k):r=2*r+(A&1);A>>=1
    need(A==0,'k-bit input');return r

def source(k,A):
    B=1<<k;need(0<=A<B,'word range')
    return 2*B-rev(k,A)

def trajectory(k,A):
    u=2;G=1;E=0;t=F(1,2);error=F(0)
    runs=[];i=0
    while i<k:
        bit=(A>>i)&1
        if bit==0:
            u*=2;t/=2;E=0;i+=1
            need(G-u*t==error,'O preserves error')
        else:
            start=t;r=0;old_error=error
            need(E==0 and 0<start<=F(1,2),'Z-run initial state')
            while i<k and (A>>i)&1:
                expected=(t+1)/2-E
                G+=u-E;u=2*u-1;t=(t+1)/2;E=1
                new_error=G-u*t
                need(new_error-error==expected,'single Z error identity')
                error=new_error;r+=1;i+=1
            increment=start+(1-start)/(1<<r)
            need(error-old_error==increment,'whole Z-run contribution')
            need(0<increment<=F(3,4),'Z-run positive upper bound')
            runs.append((start,r,increment))
    need(u==source(k,A) and G==g(u),'actual source trajectory')
    need(t==F(2*A+1,2*(1<<k)),'midpoint coordinate')
    need(error==G-u*t,'final exact error')
    return u,t,error,runs

def slope_discrepancy(slopes):
    n=len(slopes);ordered=sorted(slopes)
    return max([F(0)]+[F(i+1,n)-x for i,x in enumerate(ordered)]+[x-F(i,n) for i,x in enumerate(ordered)])

def approximate(k,theta):
    need(0<=theta<=1,'target slope')
    B=1<<k;A=min(B-1,(B*theta.numerator)//theta.denominator)
    return A,source(k,A)

def joint_discrepancies(k):
    B=1<<k
    P=[[0]*(B+1) for _ in range(B+1)]
    for A in range(B):P[A+1][rev(k,A)+1]=1
    for a in range(1,B+1):
        for b in range(1,B+1):P[a][b]+=P[a-1][b]+P[a][b-1]-P[a-1][b-1]
    minus=max(a*b-B*P[a][b] for a in range(B+1) for b in range(B+1))
    plus=max(B*P[a+1][b+1]-a*b for a in range(B) for b in range(B))
    aligned_plus=max(B*P[a][b+1]-a*b for a in range(B+1) for b in range(B))
    both=max(abs(B*P[a][b]-a*b) for a in range(B+1) for b in range(B+1))
    return F(max(minus,plus),B*B),F(max(minus,aligned_plus),B*B),F(both,B*B)

def rational(x):return f'{x.numerator}/{x.denominator}'

def main():
    words=0;runs=0;constructors=0;slopes_report=[]
    for k in range(13):
        B=1<<k;slopes=[];eps=F(3*(k+1),8*(B+1));sharp=F(3*((k+1)//2),4*(B+1))
        seen=set();max_error=F(0);max_slope_error=F(0)
        for A in range(B):
            u,t,error,zruns=trajectory(k,A);seen.add(u)
            need(u>=2 and B<u<=2*B,'actual dyadic source interval')
            need(0<=error<=F(3*len(zruns),4)<=F(3*(k+1),8),'one-sided run bound')
            need(len(zruns)<=(k+1)//2,'run count')
            slope=F(g(u),u);slopes.append(slope)
            need(t<=slope<=t+sharp<=t+eps,'slope error')
            max_error=max(max_error,error);max_slope_error=max(max_slope_error,slope-t)
            words+=1;runs+=len(zruns)
        need(seen==set(range(B+1,2*B+1)),'full actual-index permutation')
        discrepancy=slope_discrepancy(slopes)
        need(discrepancy<=F(1,2*B)+sharp<=F(1,2*B)+eps,'slope star discrepancy')
        # Cell endpoints and midpoints, all exact; the proof covers all real theta.
        for j in range(2*B+1):
            theta=F(j,2*B);A,u=approximate(k,theta);slope=F(g(u),u)
            need(theta-F(1,2*B)<=slope<=theta+F(1,2*B)+eps,'approximate-slope constructor')
            need(abs(slope-theta)<=F(1,2*B)+eps,'absolute constructor error')
            constructors+=1
        need(trajectory(k,0)[2]==0,'all-O boundary')
        if k==0:
            need(source(0,0)==2 and slopes==[F(1,2)],'empty-word boundary')
        else:
            need(source(k,B-1)==B+1 and g(B+1)==B+1,'all-Z upper slope boundary')
            need(trajectory(k,B-1)[2]==F(B+1,2*B),'all-Z exact run error')
        slopes_report.append({'k':k,'B':B,'max_error':rational(max_error),
                              'max_slope_error':rational(max_slope_error),
                              'star_discrepancy':rational(discrepancy)})
    joint=[]
    for k in range(9):
        B=1<<k;real,aligned,both=joint_discrepancies(k)
        if k==0:need(real==aligned==1 and both==0,'single origin boundary')
        else:
            need(aligned<=F(k,B),'grid-aligned first cutoff discrepancy')
            need(real<=F(k+1,B),'arbitrary real cutoff discrepancy')
        if k==1:need(real==F(3,4)>F(k,B),'hostile to unrounded joint claim')
        joint.append({'k':k,'real_star':rational(real),'first_grid_aligned':rational(aligned),'both_grid_aligned':rational(both)})
    ap_checks=0
    for k in range(7):
        B=1<<k
        for ell in range(k+1):
            M=1<<ell;step=B//M
            for q in range(B//M):
                vals=sorted(rev(k,A) for A in range(q*M,(q+1)*M))
                offset=rev(k-ell,q)
                need(vals==[offset+j*step for j in range(M)],'dyadic reverse block is AP')
                for h in range(2*B+1):
                    # y=h/(2B); compare count to M*y with integer arithmetic.
                    count=sum(2*r<h for r in vals)
                    need(abs(2*B*count-M*h)<2*B,'AP cutoff count error strictly below one')
                    ap_checks+=1
    bilinear=0
    for k in range(5):
        for ell in range(5):
            Dk=F(3*((k+1)//2),4);Dl=F(3*((ell+1)//2),4)
            for A in range(1<<k):
                u=source(k,A);tu=F(2*A+1,2*(1<<k));eu=g(u)-u*tu
                for C in range(1<<ell):
                    v=source(ell,C);tv=F(2*C+1,2*(1<<ell));ev=g(v)-v*tv
                    J=g(u)*g(v)+u*v-(u%2)*(v%2)
                    residual=J+(u%2)*(v%2)-u*v*(1+tu*tv)
                    need(residual==u*tu*ev+v*tv*eu+eu*ev,'exact bilinear residual')
                    need(0<=residual<=u*Dl+v*Dk+Dk*Dl,'bilinear approximate bound')
                    bilinear+=1
    return {'status':'PASS','bilinear_approximation_controls':bilinear,'source_words_checked':words,'Z_runs_checked':runs,
            'approximate_slope_constructor_checks':constructors,'AP_cutoff_checks':ap_checks,
            'slope_details':slopes_report,'joint_details':joint,
            'real_joint_bound_correction':'k/B requires grid-aligned first cutoff; use (k+1)/B for real cutoffs, k=0 separately',
            'exact_target_occupancy':'OPEN; this packet constructs approximate source slopes only'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
