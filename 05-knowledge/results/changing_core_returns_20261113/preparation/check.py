"""Generalized resource wrapper: independent literal subset counts and marked return checks."""
from pathlib import Path
from itertools import permutations
from math import gcd
from functools import lru_cache
import json
HERE=Path(__file__).resolve().parent
checks=0

def need(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def tournament(T):
    return all(not(row>>i&1) for i,row in enumerate(T)) and all(((T[i]>>j&1)+(T[j]>>i&1)==1) for i in range(len(T)) for j in range(i+1,len(T)))

def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def opposite(T):
    mask=(1<<len(T))-1
    return tuple(mask^row^(1<<i) for i,row in enumerate(T))

@lru_cache(None)
def tables(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n);h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for u in range(n):
            value=dp[mask][u]
            if not value:continue
            rem=T[u]&~mask
            while rem:
                bit=rem&-rem;v=bit.bit_length()-1;rem-=bit
                dp[mask|bit][v]+=value
        h[mask]=sum(dp[mask])
    return h,dp

def convolution(a,b):
    out=[]
    for mask in range(len(a)):
        sub=mask;value=0
        while True:
            value+=a[sub]*b[mask^sub]
            if sub==0:break
            sub=(sub-1)&mask
        out.append(value)
    return out

def profile(B):
    h,_=tables(B);p=convolution(h,h);t=convolution(p,h)
    return h[-1],p[-1],t[-1]

def state(T,r):
    h,e=tables(T);_,s=tables(opposite(T));full=(1<<len(T))-1
    F=sum(e[mask][r]*h[full^mask] for mask in range(1,full+1))
    return h[-1],s[-1][r],h[full^(1<<r)],F

def strong(T):
    n=len(T);full=(1<<n)-1
    for u in range(n):
        seen=1<<u
        while True:
            more=seen
            for v in range(n):
                if seen>>v&1:more|=T[v]
            if more==seen:break
            seen=more
        if seen!=full:return False
    return True

def carrier(B,m,l):
    """Root0 keeps its name while its position in C changes to l."""
    order=list(range(1,l+1))+[0]+list(range(l+1,m));x=m;offset=m+1
    T=[0]*(m+1+len(B))
    for i,u in enumerate(order):
        T[u]|=1<<x
        for v in order[i+1:]:T[u]|=1<<v
    T[x]=sum(1<<(offset+i) for i in range(len(B)))
    for i,row in enumerate(B):T[offset+i]=(row<<offset)|((1<<m)-1)
    return tuple(T)

def formula(B,m,l):
    h,p,t=profile(B);r=m-1-l
    H=((1<<m)-1)*(p-h)+2*h
    q=((1<<(m-1))-1)*(p-h)+2*h
    s=(1<<r)*h
    F=2*((1<<r)*3**l-(1<<l))*(t-p)+3*(1<<l)*p-2*(1<<l)*h-(p-h if r==0 else 0)
    return H,s,q,F

def birth(T,r,move=False):
    n=len(T);rows=[]
    for i,row in enumerate(T):rows.append(row|((1<<n) if i==r else 0))
    rows.append(((1<<n)-1)^(1<<r))
    return tuple(rows),n if move else r

def returned(T,r,k):
    for _ in range(k):T,r=birth(T,r)
    T,r=birth(T,r,True);T,r=birth(T,r)
    return T,r

def state_return(z,k):
    h,s,q,F=z;a=h-s;d=F-q;delta=d-a
    t=(1<<k)*d+q-delta;b=3*t+2*delta
    return b,t,t,b

def labelled(n):
    for mask in range(1<<(n*(n-1)//2)):
        T=[0]*n;i=0
        for a in range(n):
            for b in range(a+1,n):
                u,v=(a,b) if mask>>i&1 else (b,a);T[u]|=1<<v;i+=1
        yield tuple(T)

def permstate(T,r):
    n=len(T);H=S=F=0
    for word in permutations(range(n)):
        arcs=[bool(T[a]>>b&1) for a,b in zip(word,word[1:])]
        if all(arcs):H+=1;S+=word[0]==r
        for cut in range(1,n+1):
            if word[cut-1]==r and all(arcs[:cut-1]) and all(arcs[cut:]):F+=1
    other=[v for v in range(n) if v!=r]
    q=sum(all(T[a]>>b&1 for a,b in zip(word,word[1:])) for word in permutations(other))
    return H,S,q,F

def ord2(q):
    if q==1:return 1
    x=2%q;e=1
    while x!=1:x=2*x%q;e+=1
    return e

def oddpart(n):
    while n%2==0:n//=2
    return n

def main():
    need(profile((6,12,8,1))==(5,32,117),'strong4 source profile')
    contexts=0;permutation_controls=0
    for n in range(1,5):
        for B in labelled(n):
            need(tournament(B),'source type')
            for m in range(2,5):
                for l in range(m):
                    T=carrier(B,m,l);need(tournament(T) and strong(T),'carrier type and strongness')
                    need(state(T,0)==formula(B,m,l),('marked formula',B,m,l,state(T,0),formula(B,m,l)))
                    if n+m+1<=6:
                        need(permstate(T,0)==state(T,0),'literal permutation/cut independent check');permutation_controls+=1
                    contexts+=1;tables.cache_clear()
    B=(6,12,8,1);rows=[]
    for m in (3,5,7,9,11):
        r=m-2;z=1<<r;delta=353*z-134;d=456*z-151;q=54*z-17
        prepared=formula(B,m,1);old=formula(B,m,m-1)
        need(prepared==(108*z-17,5*z,54*z-17,510*z-168),'prepared state all-parameter formula')
        need((old[3]-old[0]-old[2]+old[1])%3==0,'last-root obstruction')
        op=oddpart(delta);period=ord2(op);step=2*period//gcd(2,period)
        need(gcd(2*prepared[3]-q,op)==1,'uniform odd-branch coprimality')
        outputs=[]
        for k in (1,1+step,1+2*step):
            out=state_return(prepared,k);b,t=out[:2]
            need(gcd(b,t)==1 and t%3==0 and b%3!=0,'native balanced entry')
            outputs.append(dict(k=k,H_bits=b.bit_length(),t_bits=t.bit_length(),gcd=gcd(b,t),H_mod3=b%3,t_mod3=t%3))
        rows.append(dict(m=m,H=prepared[0],prepared=prepared,unprepared=old,delta=delta,odd_delta=op,period=period,odd_progression_step=step,returns=outputs))
        if m<=5:
            T=carrier(B,m,1)
            for k in (1,2):
                U,root=returned(T,0,k)
                need(strong(U) and tournament(U),'literal return graph type')
                need(state(U,root)==state_return(prepared,k),'literal return state')
                tables.cache_clear()
    need(ord2(113)==28 and 58 not in {pow(2,r,113) for r in range(28)},'113 forbidden phase absent')
    need(pow(2,28,113)==1 and pow(58,28,113)==98 and pow(2,74,149)==148 and pow(81,74,149)==1,'short prime-power certificates')
    need(ord2(149)==148 and [r for r in range(148) if pow(2,r,149)==81]==[52],'149 forbidden phase even')
    need(353*966-966*353==0 and -353*319+966*134==113*149,'resultant')
    # A uniform C3 source keeps content3 even after neutral root preparation.
    C3=(2,4,1);traps=[]
    for l in range(4):
        a=formula(C3,4,l);need(a[0]==141 and gcd(*a)%3==0,'141 uniform-module content trap')
        for k in range(1,7):need(gcd(*state_return(a,k))%3==0,'all tested141 exits trapped')
        traps.append(a)
    old,oldroot=returned(carrier((0,),4,2),0,3)
    new=carrier(B,5,1);newroot=0
    need(state(old,oldroot)==(847,269,269,847) and state(new,newroot)==(847,40,415,3912),'H847 cross-template preparation')
    repaired,rr=returned(new,0,1)
    need(state(repaired,rr)==(19537,4719,4719,19537),'H847 prepared controller entry')
    oldnames=[oldroot]+[v for v in range(len(old)) if v!=oldroot][:len(new)-1]
    mapping=dict(zip(oldnames,range(len(new))))
    flips=[[u,v] for i,u in enumerate(oldnames) for v in oldnames[i+1:] if (old[u]>>v&1)!=(new[mapping[u]]>>mapping[v]&1)]
    bridge=dict(source=old,source_root=oldroot,target=new,target_root=0,retained_name_to_target_index=mapping,deleted=[v for v in range(len(old)) if v not in mapping],flipped_pairs=flips,source_state=[847,269,269,847],target_state=[847,40,415,3912],exit_state=[19537,4719,4719,19537],scope='Endpoint-neutral relabel/replace macro; same physical root; root-star may change; no claim for intermediate single edits')
    result=dict(status='PASS',H847_bridge=bridge,explicit_checks=checks,source_universe='Every labelled tournament B of orders1..4, all resource widths m2..4 and every root position; no quotient filters',carrier_contexts=contexts,independent_permutation_controls=permutation_controls,prepared_family=rows,H141_trap_states=traps,source_profile_B=[5,32,117],B_rows=B,formula_scope='All nonempty tournament B and m>=2 proved inPROOF.md; finite computations are independent controls',prime_phase_certificates={'113':{'order2':28,'forbidden':58,'found':False},'149':{'order2':148,'forbidden':81,'log2':52}})
    name='receipt.json' if __debug__ else 'receipt.optimized.json'
    (HERE/name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':'PASS','checks':checks,'contexts':contexts,'permutation_controls':permutation_controls,'targets':[a['H'] for a in rows]}))

if __name__=='__main__':main()
