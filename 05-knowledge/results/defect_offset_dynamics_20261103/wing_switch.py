"""An infinite exact observer-neutral switch, exposed by a single cloned port."""
from wing import *
from carrier import native
from return_program import literal_Q4


def check():
    rows=[];q4=0
    for high in (False,True):
        seed,_=source4(1,1,3 if high else 2,1 if high else 2)
        specs=[(compile_runs(seed,(0,)),{(1,):71,(2,):24}),
               (compile_runs(seed,(0,),starts=((0,3),)),{(1,):3,(2,):1}),
               (compile_runs(induced(seed,[i for i in range(8) if i!=3]),(0,)),{(1,):3,(2,):1}),
               (compile_runs(seed,(0,),bins=2,ends=((0,3),)),{(1,):75 if high else 87,(2,):25 if high else 29})]
        for poly,expected in specs:
            require({p:w for p,w in poly.items() if min(p)>0}==expected,'all-size selector profile certificate')
    for m in range(1,7):
        A,ad=source4(m,1,2,2);B,bd=source4(m,1,3,1)
        root=m;u=m+5;v=m+2;power=2**m
        require(native.flip(A,u,v)==B,'same labels, exactly one changed arc')
        path=(m+3,u,v,m+4,m+6)+tuple(range(m))+(m,m+1)
        require(native.legal(A,path),'uniform traced-edge path')
        X,Y,Z,p1,p2,case=native.traced_gadget(A,u,v,path)
        require(native.restrict(Z,len(A))==B,'two actual complements and explicit helper restriction')
        require(len(Z)==m+8 and case=='append_single_out','implementation clock')
        observed=[]
        for T in (A,B):
            P,E,S=bin_values(T,root,2)
            obs=(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2)
            require(obs==(24*power+23,12*power,12*power+23),'full integer observer equality')
            observed.append(obs)
            if m==1:q4+=literal_Q4(T,root,obs)
            tables.cache_clear()
        lo=birth(A,((1<<len(A))-1)^A[root])
        hi=birth(B,((1<<len(B))-1)^B[root])
        require(strong(lo) and strong(hi),'actual cloned-port query stays strong')
        require(tables(lo)[0][-1]==72*power+1 and tables(hi)[0][-1]==72*power+5,'one legal birth sees the hidden bit')
        rows.append(dict(m=m,source=A,target=B,root=root,changed_arc=(u,v),path=path,
                         native_paths=(p1,p2),helper=len(A),peak_order=len(Z),
                         observer=observed[0],phase=field(observed[0]),query_H=(72*power+1,72*power+5)))
        tables.cache_clear()
    # Closed chart: grow either left block; toggle the right-wing bit.
    states=list(product(range(1,4),range(1,4),range(2)))
    actions=list(product(range(3),range(3),range(2)))
    obs=lambda m,l,e:hwing(m,l,3 if e else 2,1 if e else 2)%7
    sig={tuple(obs((m+i-1)%3+1,(l+j-1)%3+1,e^b) for i,j,b in actions) for m,l,e in states}
    require(len(sig)==18,'minimal closed modulo-seven switch/growth carrier')
    (HERE/'wing_switch.json').write_text(json.dumps(dict(examples=rows,literal_Q4_checks=q4,future_classes=18),indent=2)+'\n')
    print('observer-neutral wing switch, six literal sizes, Q4 checks',q4,flush=True)
    return dict(examples=len(rows),literal_Q4_checks=q4,future_classes=18,scope='one-edge switch with a two-native traced implementation; helper restriction is explicit')


if __name__=='__main__':check()
