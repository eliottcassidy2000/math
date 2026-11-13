"""Exact quotient-run compiler and literal controls for the three-arc mirror."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok: raise RuntimeError(why)

def opposite(T):
    return tuple(((1<<len(T))-1)^(1<<i)^row for i,row in enumerate(T))

def inventory(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            rem=T[v]&~mask
            while value and rem:
                bit=rem&-rem;rem-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    H=[1]+[sum(row) for row in dp[1:]]
    return H,dp

def construct(M,prepared=False):
    blocks=[(0,),(0,),M,(0,),opposite(M),(0,),(0,)]
    sizes=list(map(len,blocks));starts=[];n=0
    for k in sizes:starts.append(n);n+=k
    rows=[]
    for block,start in zip(blocks,starts):
        tail=((1<<n)-1)^((1<<(start+len(block)))-1)
        rows.extend((row<<start)|tail for row in block)
    selected={0,1,4,5} if prepared else {0,3,4}
    apex=sum(((1<<sizes[i])-1)<<starts[i] for i in selected)
    T=tuple(row|((1<<n) if not(apex>>i&1) else 0) for i,row in enumerate(rows))+(apex,)
    return T,starts

def run_count(Q,r,s,observer):
    counts=[1]*8;counts[2]=r;counts[4]=s
    if observer=='delete':counts[0]=0
    @lru_cache(None)
    def visit(rem,last):
        if not any(rem):return 1
        ans=0
        for v,c in enumerate(rem):
            if c and (last<0 or Q[last]>>v&1):
                nxt=list(rem);nxt[v]-=1;ans+=visit(tuple(nxt),v)
        return ans
    if observer=='start':
        counts[0]-=1;return visit(tuple(counts),0)
    return visit(tuple(counts),-1)

def compiler(prepared):
    Q,_=construct((0,),prepared)
    return {obs:{f'{r},{s}':run_count(Q,r,s,obs) for r,s in product((1,2),repeat=2)} for obs in ('H','start','delete')}

def boundary_compiler(prepared):
    """Independent enumeration of tail assignments before/after the apex.

    1=prefix only, 2=suffix only, 3=both nonempty. Every HP has exactly
    one such assignment; within each side the ordinal blocks are forced.
    """
    wins={0,1,4,5} if prepared else {0,3,4}
    result={}
    for obs in ('H','start','delete'):
        coeff={f'{r},{s}':0 for r,s in product((1,2),repeat=2)}
        choices=[(1,2,3) if i in (2,4) else (1,2) for i in range(7)]
        if obs=='start':choices[0]=(1,)
        if obs=='delete':choices[0]=(0,)
        for assignment in product(*choices):
            prefix=[i for i,x in enumerate(assignment) if x&1]
            suffix=[i for i,x in enumerate(assignment) if x&2]
            if prefix and prefix[-1] in wins:continue
            if suffix and suffix[0] not in wins:continue
            key=','.join(str(2 if assignment[i]==3 else 1) for i in (2,4))
            coeff[key]+=1
        result[obs]=coeff
    return result

def polynomial_coefficients(coeff):
    """Coefficients in the exact basis (h^2,h*p,p^2), k=p-2*h."""
    a,b,c,d=(coeff[x] for x in ('1,1','1,2','2,1','2,2'))
    return (a-2*(b+c)+4*d,b+c-4*d,d)

def main():
    matrices=[compiler(False),compiler(True)]
    expected_polynomials=[{'H':(-1,8,12),'start':(-1,6,0),'delete':(-1,6,0)},
              {'H':(-1,8,12),'start':(0,4,4),'delete':(1,4,4)}]
    for flag,C in enumerate(matrices):
        need(C==boundary_compiler(flag),'independent apex-boundary compiler')
        need({obs:polynomial_coefficients(c) for obs,c in C.items()}==expected_polynomials[flag],
             'formal polynomial coefficient identity')
    cores=[(0,),(2,0),(6,4,0),(2,4,1),(2,4,9,3),(14,12,8,0)]
    controls=[]
    for M in cores:
        H,_=inventory(M);h=H[-1];full=len(H)-1;p=sum(H[S]*H[full^S] for S in range(full+1))
        pair=[]
        for flag in (False,True):
            T,starts=construct(M,flag);A,E=inventory(T);_,S=inventory(opposite(T));mask=(1<<len(T))-1
            obs=(A[-1],S[-1][0],A[mask^1],sum(E[m][0]*A[mask^m] for m in range(mask+1)))
            z=12*p*p+8*h*p-h*h;q=(2*p+h)**2 if flag else 6*h*p-h*h
            expected=(z,q-h*h if flag else q,q,z+h*h if flag else z)
            need(obs==expected,'literal full rooted state')
            need(sum(row&1 for row in T)==1,'root indegree one')
            pair.append(dict(rows=T,state=obs))
        T,U=pair[0]['rows'],pair[1]['rows'];a=len(T)-1
        changed=[(i,j) for i in range(len(T)) for j in range(i+1,len(T)) if ((T[i]^U[i])>>j)&1]
        need(changed==[(starts[i],a) for i in (1,3,5)],'exact three apex arcs')
        need(T[0]==U[0],'retained root star')
        need(3*pair[1]['state'][0]+2*h*h==(6*p+2*h)**2-5*h*h,'exact norm after NR')
        controls.append(dict(core=M,h=h,p=p,changed_arcs=changed,pair=pair,norm_target=(6*p+2*h)**2-5*h*h))
    receipt=dict(status='PASS',scope='Two independent universal quotient-run and apex-boundary enumerations; exact polynomial coefficients; six named literal modules, both mirrors, full root states',run_coefficients=matrices,polynomial_coefficients=expected_polynomials,controls=controls)
    (HERE/'mirror_module_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS: two symbolic run compilers and',len(controls),'literal module pairs')

if __name__=='__main__':main()
