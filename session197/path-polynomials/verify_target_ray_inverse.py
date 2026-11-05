#!/usr/bin/env python3
from verify_ratio_lattice import *

def build_counts(bank,counts):
    out=()
    for a,n in zip(bank,counts):
        for _ in range(n):out=join(out,a)
    require(out,'target modules must already be nonempty');return out

def prescribed(P,R,bank,A,C):
    require(bank and all(bank),'nonempty bank')
    hp,bp,_=profile(P);hr,br,_=profile(R);A0=2*(bp-hp)*(br-hr);C0=(bp-3*hp)*(br-3*hr)
    require(A0>0 and C0!=0,'physical target factors')
    kind,kappa=leaf_case(P,R)
    if kind=='none':return {'solvable':False,'reason':'leaf signs'}
    if not isinstance(A,int) or not isinstance(C,int) or A<=0 or A%A0 or C%C0 or C//C0<=0:
        return {'solvable':False,'reason':'target arithmetic'}
    X,Y=A//A0,C//C0;limit=X.bit_length()-1;data=[profile(a) for a in bank]
    if kind=='ratio':
        unrestricted=solve_ratios(bank,kappa)
        if not unrestricted['solvable']:return {'solvable':False,'reason':'lattice obstruction','certificate':unrestricted}
        M,t=unrestricted['matrix'],unrestricted['target']
    else:M=t=None
    def factorizations(i,x,y):
        if i==len(bank):
            if x==y==1:yield []
            return
        h,b,_=data[i];n=0
        while True:
            for rest in factorizations(i+1,x,y):yield [n]+rest
            if x%b or y%h:break
            x//=b;y//=h;n+=1
    tested_w=tested_z=0
    for w in factorizations(0,X,Y):
        require(sum(w)<=limit,'amplitude copy bound');tested_w+=1
        for z in product(*(range(-n,n+1,2) for n in w)):
            tested_z+=1
            if M is not None and matvec(M,z)!=t:continue
            u=[(n+x)//2 for n,x in zip(w,z)];v=[(n-x)//2 for n,x in zip(w,z)]
            if not sum(u) or not sum(v):continue
            U=build_counts(bank,u);V=build_counts(bank,v)
            require(coeff(P,U,R,V)==(A,0,C),'reconstructed exact target coefficients')
            return {'solvable':True,'A':A,'C':C,'w':w,'z':list(z),'u':u,'v':v,'U':U,'V':V,'copy_bound':limit,'factorizations_tested':tested_w,'allocations_tested':tested_z}
    return {'solvable':False,'reason':'bounded integer allocations exhausted','copy_bound':limit,'factorizations_tested':tested_w,'allocations_tested':tested_z}

def run():
    K=(0,);W=(8,1,3,6);P=K;R=transitive(4);bank=[K,W];cases=[]
    for power in range(13):
        A=15360*2**power;ans=prescribed(P,R,bank,A,-65)
        require(ans['solvable']==(power%2==0),'complete fixed-offset parity classification')
        if ans['solvable']:
            t=power//2;require(ans['u']==[4+t,0] and ans['v']==[t,1],'exact physical copy formula')
        cases.append({'power_of_two':power,'A':A,'solvable':ans['solvable'],'w':ans.get('w'),'u':ans.get('u'),'v':ans.get('v'),'reason':ans.get('reason')})
    physical=[]
    for t in (0,1):
        ans=prescribed(P,R,bank,15360*4**t,-65);a=graph((0,),P,ans['U'],R,ans['V']);H=h_subsets(a)[-1]
        require(H==15360*4**t*3-65 and len(a)==16+2*t,'physical target source')
        physical.append({'A':ans['A'],'C':ans['C'],'n':len(a),'H':H,'rows':a})
    good=prescribed(transitive(2),transitive(2),[K],72,1)
    require(good['solvable'] and good['u']==good['v']==[1],'nonempty positive boundary')
    a=graph((0,),transitive(2),good['U'],transitive(2),good['V']);require(h_subsets(a)[-1]==217,'physical217 boundary')
    physical.append({'A':72,'C':1,'n':len(a),'H':217,'rows':a})
    bad=prescribed(transitive(2),transitive(2),[K],36,1);require(not bad['solvable'],'no silent empty module')
    for A,C in ((15361,-65),(15360,65),(15360,0),(0,-65)):
        require(not prescribed(P,R,bank,A,C)['solvable'],'invalid target arithmetic')
    no=prescribed(P,R,[transitive(3),W],15360,-65);require(not no['solvable'] and no['reason']=='lattice obstruction','early Smith rejection')
    print(json.dumps({'fixed_offset_cases':cases,'physical_sources':physical,'nonempty_boundary_obstruction':bad,'early_lattice_obstruction':no,'status':'PASS'},sort_keys=True,indent=2))
if __name__=='__main__':run()
