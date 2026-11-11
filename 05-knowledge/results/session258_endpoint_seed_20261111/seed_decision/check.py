"""Complete fixed-multiplier seed search via nine positive SCC-cover states.
No inherited code imports; classical strong-order floor is an explicit input.
"""
from pathlib import Path
from collections import deque
from itertools import permutations
import json
HERE=Path(__file__).resolve().parent

def need(x,m):
    if not x:raise RuntimeError(m)

def tournament(n,code):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if (code>>k)&1 else (j,i)
            rows[a]|=1<<b;k+=1
    return tuple(rows)

def strong(T):
    n=len(T)
    if n==1:return True
    for rev in (False,True):
        seen=1;todo=[0]
        while todo:
            v=todo.pop()
            nxt=([u for u in range(n) if T[u]>>v&1] if rev else [u for u in range(n) if T[v]>>u&1])
            for u in nxt:
                if not (seen>>u&1):seen|=1<<u;todo.append(u)
        if seen!=(1<<n)-1:return False
    return True

def endpoint_tables(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            a=dp[S][v]
            if not a:continue
            choices=T[v]&~S
            while choices:
                bit=choices&-choices;choices-=bit;w=bit.bit_length()-1
                dp[S|bit][w]+=a
    return dp

def block_data(T,O):
    n=len(T);full=(1<<n)-1;ends=endpoint_tables(T)
    rev=tuple(sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n))
    starts=endpoint_tables(rev);L=[0]*9
    for S in range(1<<n):
        Q=full^S
        first=[1,0,0] if not S else [0,sum(ends[S][v] for v in range(n) if not(O>>v&1)),sum(ends[S][v] for v in range(n) if O>>v&1)]
        second=[1,0,0] if not Q else [0,sum(starts[Q][v] for v in range(n) if not(O>>v&1)),sum(starts[Q][v] for v in range(n) if O>>v&1)]
        for i in range(3):
            for j in range(3):L[3*i+j]+=first[i]*second[j]
    return tuple(L)

# concatenate ordinal blocks: end of first is overwritten when local first is nonempty;
# start of second is overwritten only when the old second is empty.
DEST=tuple(tuple(3*(a or i)+(j or b) for a in range(3) for b in range(3)) for i in range(3) for j in range(3))
def advance(v,L,K=None):
    z=[0]*9
    for old,x in enumerate(v):
        if not x:continue
        for loc,y in enumerate(L):
            if y:
                d=DEST[old][loc];z[d]+=x*y
    return tuple(min(x,K) for x in z) if K is not None else tuple(z)

def observables(v):return v[1]+v[2],v[8],v[2]+v[5]+v[8]
def join(A,B):
    n=len(A);m=len(B);high=((1<<m)-1)<<n
    return tuple(x|high for x in A)+tuple(x<<n for x in B)

def floor(n):
    if n==1:return 1
    need(n>=3,'no strong order two')
    q,r=divmod(n,3)
    return 3*5**(q-1) if r==0 else 5**q if r==1 else 9*5**(q-1)

def bound(h):
    n=3
    while floor(n)<=h:n+=1
    return n-1 if n>3 else 1

def make_bank(qmax):
    bank={};raw=0;classes={}
    nmax=bound(qmax+1)
    for n in [1]+list(range(3,nmax+1)):
        counts=0
        for code in range(1<<(n*(n-1)//2)):
            T=tournament(n,code)
            if not strong(T):continue
            h=sum(endpoint_tables(T)[-1])
            if h>qmax+1:continue
            counts+=1
            for O in range(1<<n):
                L=block_data(T,O);raw+=1
                if L not in bank:bank[L]={'rows':T,'selected':O,'h':h}
        classes[n]=counts
    return sorted(bank.items(),key=lambda z:(len(z[1]['rows']),z[1]['selected'],z[1]['rows'])),dict(raw_marked_presentations=raw,distinct_transfers=len(bank),labelled_strong_counts=classes,maximum_order=nmax)

def solve(q,bank):
    K=q+2;start=(1,0,0,0,0,0,0,0,0);queue=deque([start]);back={start:None};winner=None
    alphabet=[(L,b) for L,b in bank if b['h']<=q+1]
    while queue:
        v=queue.popleft();h,w,u=observables(v)
        if u==q and h==w+1 and h<=q+1:
            winner=v;break
        for idx,(L,bl) in enumerate(alphabet):
            z=advance(v,L,K);hz,wz,uz=observables(z)
            # H and U cannot decrease under adjoining an ordinal block.
            if hz>q+1 or uz>q:continue
            if z not in back:back[z]=(v,idx);queue.append(z)
    out={'q':q,'status':'SEED' if winner else 'NO_SEED','states_seen':len(back),'states_expanded':len(back)-len(queue),'alphabet_size':len(alphabet),'strong_order_bound':bound(q+1)}
    if winner:
        word=[];v=winner
        while back[v] is not None:
            prev,idx=back[v];word.append(idx);v=prev
        word.reverse();B=();O=0;literal=(1,0,0,0,0,0,0,0,0)
        for idx in word:
            L,bl=alphabet[idx];O|=bl['selected']<<len(B);B=join(B,bl['rows']);literal=advance(literal,L)
        direct=block_data(B,O)
        need(direct==literal,'independent subset DP matches block transfer')
        h,w,u=observables(direct);need(h-w==1 and u==q,'exact literal unit seed')
        tails=[]
        for t in (1,2):
            N=len(B);T=list(B)
            for _ in range(t):T=list(join(tuple(T),(0,)))
            v=len(T);newrow=O
            for a in range(v):
                if not(O>>a&1):T[a]|=1<<v
            T.append(newrow);H=sum(endpoint_tables(tuple(T))[-1])
            need(H==1+(1<<t)*q,'independent actual decorated-tail count')
            tails.append(H)
        out.update(rows=list(B),selected_mask=O,H=h,W=w,U=u,block_word=word,block_count=len(word),tail_counts=tails)
    return out

def main():
    # Exact positive composition controls, independent direct subset counts.
    composition_checks=0
    for n in range(1,4):
        for code in range(1<<(n*(n-1)//2)):
            T=tournament(n,code)
            for O in range(1<<n):
                for mark in (0,1):
                    need(advance(block_data(T,O),block_data((0,),mark))==block_data(join(T,(0,)),O|(mark<<n)),'ordinal composition control')
                    composition_checks+=1
    # Saturation must retain zeros and commute with nonnegative transition.
    saturation_checks=0
    for K in range(1,8):
        for code in range(8):
            T=tournament(3,code)
            for O in range(8):
                v=block_data(T,O)
                for mark in (0,1):
                    L=block_data((0,),mark)
                    need(advance(tuple(min(x,K) for x in v),L,K)==tuple(min(x,K) for x in advance(v,L)),'saturation homomorphism')
                    saturation_checks+=1
    bank,meta=make_bank(13)
    print('BANK',meta,flush=True)
    controls=[]
    for q in range(1,14):
        row=solve(q,bank);controls.append(row);print(q,row['status'],row['states_seen'],flush=True)
    need(controls[2]['status']=='NO_SEED','q3 would force forbidden total7')
    need(controls[9]['status']=='NO_SEED','q10 would force forbidden total21')
    out={'status':'PASS','scope':'Complete all-tournament fixed-q searches q1..13 using cited Busch bound; same code supplies general finite algorithm, not executed at unbounded q','strong_floor_input':'Busch EJC13(2006)N3 via repository ORDER_OBSTRUCTION section5; inherited not reproved','composition_checks':composition_checks,'saturation_checks':saturation_checks,'bank':meta,'results':controls}
    (HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: nine-state global seed compiler and exact literal witnesses')
if __name__=='__main__':main()
