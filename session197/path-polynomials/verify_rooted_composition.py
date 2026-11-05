#!/usr/bin/env python3
from verify_marked_contexts import substitute
from verify_matching_amplitude import *

def compose(outer,inner):
    out={}
    for (a,i),c in outer.items():
        for (b,j),d in inner.items():
            for k,e in enumerate(shift(mon(b),i)):
                if e:out[(a+k,i+j)]=out.get((a+k,i+j),0)+c*d*e
    return {key:v for key,v in out.items() if v}
def forward(op,nu):
    degree=max(a for a,j in op if j==nu);out=[0]*(degree+1)
    for (a,j),c in op.items():
        if j==nu:out[a]+=c
    return trim(out)
def scalar_shifts(op,x):
    out={}
    for (power,s),c in op.items():out[x+s]=out.get(x+s,0)+c*x**power
    return {b:c for b,c in out.items() if c}

def run():
    roots=[(tournament(n,m),r) for n in range(1,4) for m in range(1<<(n*(n-1)//2)) for r in range(n)]
    data={(a,r):(skeleton_operator(a,r),cross_matching(a,r)) for a,r in roots};checks=0
    for P,p in roots:
        AP,(nup,mup)=data[(P,p)]
        for Q,q in roots:
            AQ,(nuq,muq)=data[(Q,q)];R=substitute(P,p,Q);r=len(P)-1+q;AR=skeleton_operator(R,r)
            nur,mur=cross_matching(R,r)
            require(nur==nup+nuq and mur==mup*muq,('cross matching union',P,p,Q,q))
            require(AR==compose(AP,AQ),('full operator composition',P,p,Q,q))
            require(forward(AR,nur)==mul(forward(AP,nup),shift(forward(AQ,nuq),nup)),'top shift product');checks+=1
    Q=(2,12,145,149,99,143,175,19);q=7;AQ=skeleton_operator(Q,q);nu,mu=cross_matching(Q,q)
    require((nu,mu)==(2,4),'clock cross graph')
    require(forward(AQ,2)==[0,3,14,23,16,4],'clock leading profile')
    s1=scalar_shifts(AQ,1);s2=scalar_shifts(AQ,2)
    require(s1=={1:-1,3:60},'clock scalar1')
    require(s2=={0:6,2:-176,3:480,4:630},'clock scalar2')
    witnesses=[]
    for m in range(1,5):
        T=clone(Q,q,m);h=h_subsets(T)[-1];hd=h_subsets(wrapper(T))[-1];he=h_subsets(edge_wrapper(T))[-1]
        require(h==60*3**m-1,'physical clock')
        require(hd==630*4**m+540*3**m-176*2**m-1,'physical outer D')
        require(he==1890*4**m+1380*3**m-528*2**m+1,'physical outer E')
        witnesses.append({'m':m,'H':h,'D_H':hd,'E_H':he})
    print(json.dumps({'rooted_pairs':checks,'clock_rows':Q,'clock_root':q,'clock_matching_number':nu,'maximum_matchings':mu,'top_forward_polynomial':forward(AQ,nu),'scalar_shifts_at_1':s1,'scalar_shifts_at_2':s2,'physical_clone_controls':witnesses},sort_keys=True,indent=2))
if __name__=='__main__':run()
