"""Literal two-entry programs, joint-cover transport, and an infinite source switch."""
from research import *
from verify import rooted_data


def load(name,relative):
    spec=importlib.util.spec_from_file_location(name,HERE.parent/relative)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


graded=load('inherited_graded','source_selection_closure_20261103/graded_carrier.py')
golden=load('inherited_golden','golden_inverse_defects_20261103/inverse.py')
native=load('inherited_native','adaptive_metagraph_reformulation_20261103/metagraph_probe.py')


def falling(a,k):
    result=1
    for i in range(k):result*=a-i
    return result


def extract(vector,m,slot=0):
    P=E=S=C=0
    for (skel,anonymous),weight in vector.items():
        c=len(skel)+anonymous
        P+=weight*falling(m,c)
        E+=weight*sum(b==slot for a,b in skel)*falling(m-1,c-1)
        S+=weight*sum(a==slot for a,b in skel)*falling(m-1,c-1)
        pairs=sum(i!=j and b>=0 and a>=0 for i,(_,b) in enumerate(skel) for j,(a,_) in enumerate(skel))
        if pairs:C+=weight*pairs*falling(m-1,c-2)
    return P,E,S,C  # C is the two-port C_J(m+1), not a marginal product.


def observer(vector,slot):
    h,e,_,_=extract(vector,1,slot);_,_,s2,_=extract(vector,2,slot)
    require((s2-e)%2==0,'golden integral parity')
    return h,(s2-e)//2,(s2+e)//2


def field(v):
    h,x,y=v
    return (4*(h+2*x-y))%7,y%7,(h-y)%7


def step(T,ports,vector,command):
    U,p=graded.birth(T,ports,*command)
    w=graded.apply(vector,2,graded.operator(2,command))
    return U,p,w


def compile_funnel(T,ports,J,word,ready=True):
    """Retain new root in slot 0 and the old second physical port in slot 1."""
    require(set(J)<=set(ports),'named entry ports')
    vector=graded.literal(T,ports)
    mask=sum(1<<j for j,v in enumerate(ports) if v not in J)
    T,ports,vector=step(T,ports,vector,(1,mask,0))
    for bit in word:T,ports,vector=step(T,ports,vector,(1,2+(1-bit),2))
    if ready:T,ports,vector=step(T,ports,vector,(0,1,1))
    return T,ports,vector


def compile_reset(T,ports,vector,slot):
    other=1-slot
    T,ports,vector=step(T,ports,vector,(0,1<<slot,other))
    T,ports,vector=step(T,ports,vector,(0,1<<other,2))
    return T,ports,vector,other


def common_source(t,lane):
    require(t>=1 and lane in (1,2),'family range')
    word=(1,)*t+((1,0,1) if lane==1 else (0,0,1))
    return compile_funnel((2,4,1),(0,1),{0} if lane==1 else {0,1},word)


def switched_by_native(t):
    T,ports,v=common_source(t,1);target,target_ports,w=common_source(t,2)
    # Construction order: C3=0,1,2; old funnel root=3; tail=4,...; ready root last.
    r=3;tail=tuple(range(4,t+7));x=t+7;u=tail[t]
    require(ports==target_ports==(r,x),'physical alignment')
    paths=((1,r,x,2,0)+tail,tail[:t+1]+(r,x,0,1,2)+tail[t+1:])
    intermediate=[]
    for a,p in zip((1,u),paths):
        require(native.legal(T,p) and (a,r) in tuple(zip(p,p[1:])),'uniform traced guard')
        A,B,C,p1,q,case=native.traced_gadget(T,a,r,p)
        require(case=='append_single_out','one out-port temporary birth')
        T=native.restrict(C,len(T))
        intermediate.append(dict(arc=[a,r],birth_out=[a],first_path=p1,second_path=q))
    require(T==target,'exact two-edge switch target')
    return dict(t=t,n=len(T),H=90*2**t+5,commands=intermediate,
                source_phase=field(observer(v,1)),target_phase=field(observer(w,1)))


def matmul(A,B):
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*B)) for row in A)


def tail_matrix(bit,m):
    states=list(product(range(3),range(2)));M=[[0]*6 for _ in states]
    for i,(last,started) in enumerate(states):
        M[i][i]+=m-1
        M[i][states.index((1 if bit else 2,started))]+=1
        if started or not bit:M[i][states.index((last,1))]+=1
    return tuple(map(tuple,M))


def main():
    extraction=literal_funnels=scalar_controls=0
    for T in ((2,4,1),(2,4,9,3)):
        for ports in combinations(range(len(T)),2):
            vector=graded.literal(T,ports)
            for m in (1,2,3):
                for slot in (0,1):
                    actual=bin_values(T,ports[slot],3)
                    p,e,s,c=extract(vector,m,slot)
                    require((p,e,s)==tuple(v[m] for v in actual),'bin extraction')
                    require(c==joint_values(T,set(ports),m),'joint extraction')
                    extraction+=1
    for J in ({0},{0,1}):
        for word in ((1,),(0,1),(1,0,1)):
            T,p,v=compile_funnel((2,4,1),(0,1),J,word,ready=False)
            U,r=funnel((2,4,1),J,word)
            # Match the old implementation's core, tail, root order.
            reorder=list(range(3))+list(range(4,len(T)))+[3]
            require(induced(T,reorder)==U,'literal funnel alignment')
            if len(T)<=6:
                require(graded.literal(T,p)==v,'all joint coefficients literal')
                literal_funnels+=1
            T,p,v=step(T,p,v,(0,1,1))
            require(rooted_data(T,p[1])[3]-rooted_data(T,p[1])[0]==len(J),'ready gap')
            for slot in (0,1):
                actual=bin_values(T,p[slot],2)
                for m in (1,2):require(extract(v,m,slot)[:3]==tuple(z[m] for z in actual),'actual ready bins')
            R,rp,rv,rs=compile_reset(T,p,v,1)
            require(graded.scalar(rv)==tables(R)[0][-1]==3*tables(T)[0][-1]+2*len(J),'native birth reset scalar')
            scalar_controls+=1;tables.cache_clear()
    # This is a symbolic operator identity checked at several integer m;
    # the all-m proof uses two commuting idempotents in TRANSPORT.md.
    annihilations=0
    for m in range(1,6):
        for bit in (0,1):
            M=tail_matrix(bit,m);I=tuple(tuple(int(i==j) for j in range(6)) for i in range(6))
            R=I
            for lam in (m-1,m,m+1):
                R=matmul(R,tuple(tuple(M[i][j]-lam*I[i][j] for j in range(6)) for i in range(6)))
            require(not any(any(row) for row in R),'three roots')
            annihilations+=1
    # Exact common family, with retained source-specific phases and all Q4 maps.
    phases={1:[],2:[]};family=[];q4literal=0
    for t in range(1,13):
        row=dict(t=t,H=90*2**t+5,lanes=[])
        for lane in (1,2):
            T,p,v=common_source(t,lane);obs=observer(v,1)
            ax,bx,cx=(147,132,-15) if lane==1 else (238,136,-54)
            xx=ax*3**t+bx*2**t+cx
            require(obs==(row['H'],xx,xx+3) and extract(v,2,1)[1]-obs[0]==lane,'uniform common source and full observer')
            R,rp,rv,rs=compile_reset(T,p,v,1)
            require(graded.scalar(rv)==3*row['H']+2*lane,'uniform branches')
            phases[lane].append(field(obs))
            q4={name:field(golden.step(obs,golden.BANK[name])) for name in ('q0','q1','q2','q4','q5')}
            a,b,g=field(obs)
            require(q4['q2']==(2*b%7,(a+4*b)%7,g),'F49 transport')
            if t<=3:
                require(rooted_data(T,p[1])[0]==row['H'],'independent family H')
                require(tables(R)[0][-1]==3*row['H']+2*lane,'independent common-family branch H')
            row['lanes'].append(dict(lane=lane,ports=p,observer=obs,phase=field(obs),Q4_phase=q4,reset_H=graded.scalar(rv)))
            tables.cache_clear()
        family.append(row)
    require(all(phases[b][:6]==phases[b][6:] for b in (1,2)),'six-clock phase')
    require(all((phases[1][t-1]==phases[2][t-1])==(t%3==0) for t in range(1,13)),'exact phase-neutral clock')
    # Literal Q4 attachments on modest new sources; all five presentations kept.
    for lane in (1,2):
        T,p,v=compile_funnel((2,4,1),(0,1),{0} if lane==1 else {0,1},(1,))
        obs=observer(v,1);order=[p[1]]+[i for i in range(len(T)) if i!=p[1]];T=induced(T,order)
        for name in ('q0','q1','q2','q4','q5'):
            block=golden.literal((name,));d=len(block)
            rows=list(block)+[(row<<d)|((1<<d)-1) for row in T]
            rows[d-1]|=1<<d;rows[d]&=~(1<<(d-1))
            P,E,S=bin_values(tuple(rows),0,2)
            actual=(P[1],(S[2]-E[1])//2,(S[2]+E[1])//2)
            require(actual==golden.step(obs,golden.BANK[name]),('literal Q4',lane,name))
            tables.cache_clear();q4literal+=1
    switches=[switched_by_native(t) for t in (1,2,3,6,20,100)]
    from math import isqrt
    prime_controls={p:all(p%d for d in range(2,isqrt(p)+1)) for p in (311,2179)}
    require(all(prime_controls.values()) and 3*(90*8+5)+2==7*311,'first phase-fixed prime lanes')
    # Long checks use the proved polynomial program, not subset enumeration.
    long=[]
    for lane in (1,2):
        T,p,v=common_source(100,lane)
        require(graded.scalar(v)==90*2**100+5,'long exact family')
        long.append(dict(lane=lane,n=len(T),H=graded.scalar(v),phase=field(observer(v,1)),joint_nonzero_coefficients=len(v)))
    print(json.dumps(dict(status='PASS',endpoint_skeletons=len(graded.skeletons(2)),
                         extraction_controls=extraction,literal_joint_funnels=literal_funnels,
                         independent_reset_controls=scalar_controls,tail_matrix_checks=annihilations,
                         literal_Q4_controls=q4literal,common_family=family,
                         source_switches=switches,long_joint_programs=long,
                         trial_division_prime_controls=prime_controls),indent=2,sort_keys=True))


if __name__=='__main__':main()
