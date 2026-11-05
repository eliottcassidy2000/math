#!/usr/bin/env python3
from verify import *
Q=(14,32,2,22,103,13,47)

def graph(T,P,U,R,V):
    blocks=[T,P,U,(0,),R,(0,),V];starts=[];n=0
    for a in blocks:starts.append(n);n+=len(a)
    out=[0]*n
    for i,a in enumerate(blocks):
        for k,row in enumerate(a):
            u=starts[i]+k;out[u]=row<<starts[i]
            for j,b in enumerate(blocks):
                if Q[i]>>j&1:out[u]|=((1<<len(b))-1)<<starts[j]
    validate(out);return tuple(out)

def profile(a):
    p=bin_polynomial(a);return value(p,1),value(p,2),value(p,3)
def coeff(P,U,R,V):
    hp,bp,_=profile(P);hu,bu,_=profile(U);hr,br,_=profile(R);hv,bv,_=profile(V)
    return (2*bu*bv*(bp-hp)*(br-hr),
            bv*(br-hr)*(4*hp-bp)*hu+bu*(bp-hp)*(4*hr-br)*hv,
            hu*hv*(bp-3*hp)*(br-3*hr))

def run():
    physical=0
    for T,P,U,R,V in product(((0,),(2,4,1)),repeat=5):
        A,B,C=coeff(P,U,R,V);ht,bt,bt3=profile(T)
        actual=h_subsets(graph(T,P,U,R,V))[-1]
        require(actual==A*bt3+B*bt+C*ht,('arbitrary module formula',T,P,U,R,V,actual,A,B,C));physical+=1
    unit_controls=[]
    for u,v,m in product((1,2),(1,2),(1,2,3)):
        P=R=transitive(2);U=transitive(u);V=transitive(v);T=transitive(m)
        A,B,C=coeff(P,U,R,V);actual=h_subsets(graph(T,P,U,R,V))[-1]
        require((A,B,C)==(18*2**(u+v),0,1),'unit coefficients')
        require(actual==A*3**m+1,'unit physical control')
        unit_controls.append({'u':u,'v':v,'m':m,'A':A,'H':actual})
    cuts=0;plus=[];minus=[]
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);h,b,_=profile(a)
            require(b>=(n+1)*h,('cut inequality',a,h,b))
            if b-3*h==1:plus.append((n,mask))
            if b-3*h==-1:minus.append((n,mask))
            cuts+=1
    require(plus==[(2,0),(2,1)] and minus==[(1,0)],'boundary inventory')
    choices=[(0,),transitive(2),transitive(3),(2,4,1),wrapper(transitive(2))];classification=0
    for P,U,R,V in product(choices,repeat=4):
        A,B,C=coeff(P,U,R,V)
        if C in (-1,1):
            hp,bp,_=profile(P);hr,br,_=profile(R);hu,_,_=profile(U);hv,_,_=profile(V)
            require(hu==hv==1 and len(P)<=2 and len(R)<=2,'unit boundary necessity')
            expected=(len(P)==len(R)==2)
            require((B==0)==expected,'complete unit cancellation')
            if expected:require(C==1 and A==18*2**(len(U)+len(V)),'unit amplitude')
        classification+=1
    print(json.dumps({'arbitrary_module_graphs':physical,'positive_unit_controls':unit_controls,'cut_inequality_graphs':cuts,'plus_one_boundary':plus,'minus_one_boundary':minus,'physical_profile_classification_cases':classification},sort_keys=True,indent=2))
if __name__=='__main__':run()
