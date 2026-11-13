"""An actual five-base guarded grammar, with literal full-graph controls."""
from context import *
from synchrony_probe import expand
from numba import njit
from collections import Counter

@njit(cache=True)
def count_capped(rows,cap):
    n=len(rows);N=1<<n;dp=np.zeros((N,n),dtype=np.int64)
    for v in range(n):dp[1<<v,v]=1
    for mask in range(1,N):
        for v in range(n):
            c=dp[mask,v]
            if not c:continue
            for w in range(n):
                if not mask>>w&1 and rows[v]>>w&1:
                    target=mask|1<<w;old=dp[target,w]
                    dp[target,w]=cap if old>=cap-c else old+c
    total=0
    for v in range(n):
        c=dp[N-1,v]
        total=cap if total>=cap-c else total+c
    return total

def state(r,m,k,s=None):
    need(m>=1 and k>=1,'positive synchronized and parity-phase parameters')
    z=2*m+r['z_offset']+k;j=2*m+r['A_valuation_offset']
    c=r['anchor']+2**k*sum(d*b**m for b,d in r['B_odd'])
    seed_order=len(r['rows'])+2*(m-1)+z-1
    out={'family':r['family'],'m':m,'k':k,'z':z,'j':j,'c':c,'seed_order':seed_order}
    if s is not None:
        need(s>=j+1,'actual ray head guard');out.update(s=s,H=1+c*2**s,order=seed_order+s-j+1)
    return out

def construct(r,m,k,s):
    st=state(r,m,k,s);rows,O,_=expand(r['rows'],r['selected'],dict(zip(r['roots']+[r['clock']],[m,m,st['z']])))
    n=len(rows);t=s-st['j'];apex=n+t
    full=[row|sum(1<<v for v in range(n,n+t))|((1<<apex) if not O>>i&1 else 0) for i,row in enumerate(rows)]
    full.extend(sum(1<<w for w in range(v+1,n+t))|(1<<apex) for v in range(n,n+t));full.append(O)
    need(all(not row>>i&1 for i,row in enumerate(full)) and all(bool(full[i]>>j&1)!=bool(full[j]>>i&1) for i in range(len(full)) for j in range(i+1,len(full))),'actual complete loop-free orientation')
    return st,full

def recognize(c,bank):
    """Exact scalar recognition for these nine charts, not all unit seeds."""
    hits=[]
    for name,r in bank.items():
        t=c-r['anchor']
        if t<=0 or t%2:continue
        k=v2(t)
        if k<1:continue
        target=t>>k
        def B(m):return sum(d*b**m for b,d in r['B_odd'])
        lo=1;hi=1
        while B(hi)<target:hi*=2
        while lo<hi:
            mid=(lo+hi)//2
            if B(mid)<target:lo=mid+1
            else:hi=mid
        if B(lo)==target:hits.append((name,lo,k))
    return hits

def run():
    data=load('synchronized_descent.json.gz',HERE);bank={r['family']:r for r in data['records']};parents={r['child_family']:r['parent_family'] for r in data['cross_core_links']}
    terminals=sorted(set(bank)-set(parents));need(len(terminals)==5,'five terminal charts')
    checks=0;byorder=Counter();diamonds=0;recognition=0
    def step(st,op):
        family=parents[st['family']] if op=='core' else st['family']
        return state(bank[family],st['m']-int(op=='m'),st['k']-int(op=='k'))
    for name,r in bank.items():
        for m in range(1,11):
            for k in range(1,11):
                a=state(r,m,k);byorder[a['seed_order']]+=1
                need((name,m,k) in recognize(a['c'],bank),'scalar recognition recovers actual family and coordinates');recognition+=1
                moves=[]
                if k>1:moves.append(state(r,m,k-1))
                if m>1:moves.append(state(r,m-1,k))
                if name in parents:moves.append(state(bank[parents[name]],m,k))
                for b in moves:
                    need(b['c']<a['c'] and b['seed_order']<a['seed_order'] and b['j']<=a['j'],'decreasing legal coefficient, size and threshold')
                    checks+=1
                ops=(['m'] if m>1 else [])+(['k'] if k>1 else [])+(['core'] if name in parents else [])
                for i,op in enumerate(ops):
                    for other in ops[i+1:]:
                        need(step(step(a,op),other)==step(step(a,other),op),'entire state diamond including j and z');diamonds+=1
    # Literal verification of an actual even parent and its one-exit-small child.
    # The count is saturated at target+1, so no fixed-width overflow can hide
    # a mismatch even for the order21 child.
    literal=[]
    for name in ['even:0','odd:0']:
        st,rr=construct(bank[name],1,1,8);target=st['H']
        need(target<2**62,'safe cap');actual=int(count_capped(np.array(rr,dtype=np.int64),target+1))
        need(actual==target,'independent literal complete Hamilton count')
        literal.append({**st,'rows':rr,'actual_H':actual})
        print('literal',name,st['order'],actual,flush=True)
    weights=Counter(len(r['rows'])+2+r['z_offset'] for r in bank.values())
    need(weights=={17:3,18:2,19:4},'formal-state GF numerator')
    seq=[sum(mult*(1+(n-b)//4) for b,mult in weights.items() if n>=b) for n in range(17,37)]
    need(not recognize(29,bank) and (5*29+24)%13==0,'congruence alone does not certify chart membership')
    out={'scope':'Nine certified infinite charts and their four literal core links; inverse k, synchronized m, and parent moves. Formal grammar states are distinguished, not unmarked tournament isomorphism classes.','terminal_families':terminals,'terminal_coefficients':{n:state(bank[n],1,1)['c'] for n in terminals},'formal_state_GF':'(3*z^17+2*z^18+4*z^19)/((1-z)*(1-z^4))','formal_counts_orders17_through36':seq,'finite_grid_states':900,'decreasing_edge_controls':checks,'commuting_diamond_controls':diamonds,'scalar_recognition_controls':recognition,'congruence_only_hostile':29,'literal_complete_counts':literal}
    save('descent_grammar.json.gz',out);print({k:v for k,v in out.items() if k!='literal_complete_counts'},flush=True)
if __name__=='__main__':run()
