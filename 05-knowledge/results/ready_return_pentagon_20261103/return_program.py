"""Two exact return programs, joint-port transport, and clone-run certificates."""
from explore_returns import *
from carrier import graded,golden,native,extract,observer,field,step
from math import comb,factorial


def converse(vector):
    out=Counter()
    for (skel,anon),weight in vector.items():
        out[tuple(sorted((b,a) for a,b in skel)),anon]+=weight
    return out


def source(k,lane):
    T=(2,0);ports=(0,1);v=graded.literal(T,ports)
    for _ in range(k):T,ports,v=step(T,ports,v,(1,0,2))
    T,ports,v=step(T,ports,v,(0,3,0 if lane==1 else 1))
    if lane==1:T,v=opposite(T),converse(v)
    return T,ports,v


def returned(k,lane):
    T,p,v=source(k,lane);input_observer=observer(v,1)
    # Keep the spare in slot zero throughout; overwrite the selected slot by x.
    T,p,v=step(T,p,v,(0,2,1))
    T,p,v=step(T,p,v,(0,2,2))
    sealed_observer=observer(v,1)
    T,v=opposite(T),converse(v)
    return T,p,v,input_observer,sealed_observer


def project(vector,keep):
    rename={slot:i for i,slot in enumerate(keep)};out=Counter()
    for (skel,anon),weight in vector.items():
        parts=tuple((rename.get(a,-1),rename.get(b,-1)) for a,b in skel)
        removed=sum(part==(-1,-1) for part in parts)
        out[tuple(sorted(part for part in parts if part!=(-1,-1))),anon+removed]+=weight
    return out


def common_three(k):
    T,r=incoming_pair_source(((1,1),));ports=(0,1,r);v=graded.literal(T,ports)
    for _ in range(1,k):
        T,ports=graded.birth(T,ports,0,4,3)
        v=graded.apply(v,3,graded.operator(3,(0,4,3)))
    order=[0,1]+list(reversed(range(4,len(T))))+[2,3]
    expected,_=incoming_pair_source(((1,1),)*k)
    require(induced(T,order)==expected,'common three-port physical source')
    for lane,keep in ((1,(2,1)),(2,(0,2))):
        w=project(v,keep)
        if lane==1:w=converse(w)
        _,_,target=source(k,lane)
        require(w==target,'joint projection supplies each lane without recovering forgotten ports')
    return dict(k=k,n=len(T),ports=ports,joint_terms=len(v),skeleton_types=len(graded.skeletons(3)))


def compositions(n,m):
    if m==1:yield (n,);return
    for a in range(n+1):
        for rest in compositions(n-a,m-1):yield (a,)+rest


def clone_coefficients(T,clone,root,m):
    """Exact unlabeled-run skeletons; coefficient j multiplies j! S(k,j)."""
    external=[a for a in range(len(T)) if a!=clone]
    answer=[[0]*(len(external)+m+1) for _ in range(3)] # P,E,S
    for perm in permutations(external):
        for sizes in compositions(len(external),m):
            pieces=[];offset=0
            for size in sizes:pieces.append(perm[offset:offset+size]);offset+=size
            for mask in range(1<<(len(external)+m)):
                paths=[];gap=0
                for part in pieces:
                    path=[]
                    for i in range(len(part)+1):
                        if mask>>gap&1:path.append(clone)
                        gap+=1
                        if i<len(part):path.append(part[i])
                    paths.append(path)
                if not all(T[a]>>b&1 for path in paths for a,b in zip(path,path[1:])):continue
                j=mask.bit_count();answer[0][j]+=1
                if paths[0] and paths[0][-1]==root:answer[1][j]+=1
                if paths[0] and paths[0][0]==root:answer[2][j]+=1
    return answer


def onto(k,j):
    return sum((-1)**(j-i)*comb(j,i)*i**k for i in range(j+1))


def evaluate(coeff,k):return sum(c*onto(k,j) for j,c in enumerate(coeff))


def phase_closure(lane):
    """Synthetic phase language only: affine growth and field Frobenius."""
    def sigma(v):a,b,g=v;return ((a+4*b)%7,-b%7,g)
    shift=(-2,0,-6) if lane==1 else (-1,2,-3)
    def grow(v):return tuple((2*x+t)%7 for x,t in zip(v,shift))
    def ungrow(v):return tuple((4*(x-t))%7 for x,t in zip(v,shift))
    T,p,v,_,_=returned(1,lane);start=field(observer(v,0))
    states=[start];seen={start}
    for state in states:
        for nxt in (grow(state),sigma(state)):
            if nxt not in seen:seen.add(nxt);states.append(nxt)
    require(len(states)==(6 if lane==1 else 21),'synthetic phase closure cardinality')
    delta=(0,0,0) if lane==1 else (1,3,0)
    for a,b,g in states:
        state=(a,b,g)
        require(sigma(grow(sigma(ungrow(state))))==tuple((x+y)%7 for x,y in zip(state,delta)),
                'Frobenius-growth commutator is the specified translation')
        normal=((b+2*a)%7,2*a%7)
        sa,sb,_=sigma(state)
        require(((sb+2*sa)%7,2*sa%7)==normal[::-1],'literal normal-basis swap')
    # Exact permutation group on this scope, with independent composition closure.
    states=sorted(states);ix={s:i for i,s in enumerate(states)}
    gens=[tuple(ix[f(s)] for s in states) for f in (grow,sigma)]
    identity=tuple(range(len(states)));group=[identity];known={identity}
    for perm in group:
        for gen in gens:
            nxt=tuple(gen[i] for i in perm)
            if nxt not in known:known.add(nxt);group.append(nxt)
    require(len(group)==(6 if lane==1 else 42),'synthetic operation group')
    quotient=None
    if lane==2:
        phi=lambda s:(5*s[1]+3)%7
        fibres=Counter(phi(s) for s in states)
        require(fibres==Counter({i:3 for i in range(7)}),'literal three-to-one Paley quotient')
        for s in states:
            require(phi(grow(s))==2*phi(s)%7 and phi(sigma(s))==(-phi(s)-1)%7,'equivariant quotient')
        paley=tuple(sum(1<<j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7))
        automorphisms=[p for p in permutations(range(7)) if induced(paley,p)==paley]
        affine={tuple((a*i+b)%7 for i in range(7)) for a in (1,2,4) for b in range(7)}
        require(set(automorphisms)==affine,'exact Paley automorphism group')
        quotient=dict(map='x=5b+3 modulo7',fibre_size=3,Paley=paley,automorphism_order=21,
                      growth='x -> 2x',Frobenius='x -> -x-1')
    return dict(lane=lane,scope='synthetic affine growth/Frobenius, no general source lift claimed',
                states=states,state_count=len(states),group_order=len(group),commutator=delta,Paley_quotient=quotient)


def literal_Q4(T,root,v):
    checks=0;T=induced(T,[root]+[i for i in range(len(T)) if i!=root])
    for name in ('q0','q1','q2','q4','q5'):
        block=golden.literal((name,));d=len(block)
        rows=list(block)+[(row<<d)|((1<<d)-1) for row in T]
        rows[d-1]|=1<<d;rows[d]&=~(1<<(d-1))
        p,e,s=bin_values(tuple(rows),0,2)
        actual=p[1],(s[2]-e[1])//2,(s[2]+e[1])//2
        require(actual==golden.step(v,golden.BANK[name]),'all five literal Q4 presentations')
        checks+=1;tables.cache_clear()
    return checks


def pentagon(k):
    tails=((),((1,0),),((1,0),(0,1)),((0,1),(1,0)),((0,1),))
    graphs=[incoming_pair_source(((1,1),)*k+tail)[0] for tail in tails]
    AB=graphs[2];u=k+2;v=k+3;r=k+4
    path=(0,r)+tuple(range(2,k+2))+(u,v,1)
    require(native.legal(AB,path),'uniform pentagon swap guard')
    _,_,C,p1,q,case=native.traced_gadget(AB,u,v,path)
    flip=native.restrict(C,len(AB))
    order=[0,1]+list(range(2,k+2))+[v,u,r]
    require(induced(flip,order)==graphs[3],'literal AB to BA target')
    require(case=='append_single_out','two native moves, one temporary vertex')
    # Prefix growth commutes with every displayed boundary word operation.
    for tail,T in zip(tails,graphs):
        r=len(T)-1;U,x=birth(T,r)
        expected,_=incoming_pair_source(((1,1),)*(k+1)+tail)
        require(induced(U,[0,1,x]+list(range(2,r))+[r])==expected,'pentagon growth transport')
    return dict(k=k,swap_path=path,first_native_path=p1,second_native_path=q,
                persistent_ports_in_AB_frame=[0,1,k+4],peak_order=k+6)


def main():
    family=[];coefficients={};literal=direct=phase=0
    for lane in (1,2):
        T,p,v,_,_=returned(1,lane)
        one=clone_coefficients(T,2,p[0],1);two=clone_coefficients(T,2,p[0],2)
        coefficients[lane]=dict(H=one[0],E1=one[1],S2=two[2])
    for k in range(1,13):
        for lane in (1,2):
            T,p,v,inp,sealed=returned(k,lane);obs=observer(v,0)
            a,b=0,1;r=k+2;x=k+3;y=k+4
            if lane==1:
                word=((1,1),)*k+((0,0),(0,1))
                order=[a,b]+list(range(2,k+2))+[y,x,r]
            else:
                word=((1,0),)*k+((0,0),(1,0))
                order=[r,b]+list(reversed(range(2,k+2)))+[y,x,a]
            expected,root=incoming_pair_source(word)
            require(induced(T,order)==expected,'all old arcs and return-root alignment')
            require(order[-1]==p[0] and p[1]==x,'spare port retained through both births')
            h=9*2**k+(-1 if lane==1 else 1)
            require(obs[0]==h==graded.scalar(v)==3*inp[0]+2*lane,'both exact ternary outputs')
            co=coefficients[lane]
            hv,e,s=(evaluate(co[key],k) for key in ('H','E1','S2'))
            require(obs==(hv,(s-e)//2,(s+e)//2),'independent clone skeleton formula')
            xx=12*2**k-8 if lane==1 else 7*3**k+4*2**k-4
            require(obs==(h,xx,xx+lane),'closed golden observer')
            if k<=3:
                profile=rooted_data(T,p[0]);actual=bin_values(T,p[0],2)
                require(extract(v,1,0)[:3]==tuple(z[1] for z in actual),'independent endpoint bins')
                require(extract(v,2,0)[:3]==tuple(z[2] for z in actual),'independent two-bin marginals')
                expected_profile=(h,3,3,h+4) if lane==1 else (h,3*2**k-1,3*2**k+3,12*2**k+5)
                require(profile==expected_profile,'complete returned root profile')
                if k==1:
                    require(graded.literal(T,p)==v,'every joint coefficient by literal covers')
                    literal+=1
                    direct+=literal_Q4(T,p[0],obs)
                tables.cache_clear()
            q4={name:field(golden.step(obs,golden.BANK[name])) for name in ('q0','q1','q2','q4','q5')}
            aa,bb,g=field(obs)
            require(q4['q2']==(2*bb%7,(aa+4*bb)%7,g),'literal theta multiplication')
            frob=((aa+4*bb)%7,-bb%7,g)
            require(((frob[0]+4*frob[1])%7,-frob[1]%7,g)==field(obs),'Frobenius involution')
            qa,qb,_=q4['q2']
            require(((qa+4*qb)%7,-qb%7,g)==((4*frob[0]-2*frob[1])%7,-frob[0]%7,g),'Frobenius conjugates theta transport')
            family.append(dict(k=k,lane=lane,H=h,ports=p,return_word=word,
                source_phase=field(inp),sealed_phase=field(sealed),returned_observer=obs,
                returned_phase=field(obs),Frobenius_phase=frob,Q4_phase=q4))
    for lane in (1,2):
        rows=[r for r in family if r['lane']==lane]
        require([r['returned_phase'] for r in rows[:6]]==[r['returned_phase'] for r in rows[6:]],'six-clock phase')
        require([r['returned_phase'] for r in rows[:3]]==[r['returned_phase'] for r in rows[3:6]],'actual shorter three-clock phase')
        phase+=1
    # A literal F49 phase is a projection; the joint integer vectors are retained.
    long=[]
    for lane in (1,2):
        T,p,v,inp,_=returned(100,lane)
        require(graded.scalar(v)==9*2**100+(-1 if lane==1 else 1),'long actual return program')
        co=coefficients[lane];h,e,s=(evaluate(co[key],100) for key in ('H','E1','S2'))
        require(observer(v,0)==(h,(s-e)//2,(s+e)//2),'long independent run skeletons')
        long.append(dict(lane=lane,n=len(T),phase=field(observer(v,0)),joint_terms=len(v)))
    precision=[]
    for a in range(1,5):
        modulus=7**a
        for lane in (1,2):
            length=3*7**(a-1)*(2 if lane==2 and a>=2 else 1)
            def residue(k):
                t=pow(2,k,modulus)
                h=9*t+(-1 if lane==1 else 1)
                xx=12*t-8 if lane==1 else 7*pow(3,k,modulus)+4*t-4
                return (h%modulus,xx%modulus,(xx+lane)%modulus)
            require(all(residue(k)==residue(k+length) for k in range(1,length+1)),'precision period')
            for prime in (2,3,7):
                if length%prime==0:require(any(residue(k)!=residue(k+length//prime) for k in range(1,8)),'exact precision period')
            precision.append(dict(lane=lane,exponent=a,modulus=modulus,exact_period=length))
    print(json.dumps(dict(status='PASS',return_families=family,clone_coefficients=coefficients,
        literal_full_joint_checks=literal,literal_Q4_checks=direct,six_clock_controls=phase,
        pentagons=[pentagon(k) for k in (1,2,3,20,100)],long_returns=long,
        precision_clocks=precision,synthetic_phase_closures=[phase_closure(lane) for lane in (1,2)],
        common_three_port_controls=[common_three(k) for k in (1,2,3,6,12)]),indent=2,sort_keys=True))


if __name__=='__main__':main()
