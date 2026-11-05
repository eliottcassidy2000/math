"""Synchronized joint inverse verifier. Python stdlib; explicit checks only."""
import argparse, hashlib, heapq, json, math
from pathlib import Path
from itertools import product

T=(1,0,2,1); B=(3,1,3,1); Z0=(1,0,2,0)
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(w,s=T):
    for c in w:s=step(s,c)
    return s
def K(v,w):return sum(v[i]*Q[i][j]*w[j] for i in range(4)for j in range(4))
def fastK(v,w):
    D,E,U,d=v; A,B,C,e=w
    return D*e+d*A+E*B+U*C-E*e-d*B-U*e-d*C+d*e

def F(s):return s[0]-s[1]+s[2]
def G(s):return s[0]+2*s[2]
def H(s):return 2*s[2]-s[3]
def dagger(w):return ''.join({'O':'Z','Z':'O','J':'J'}[x]for x in w[::-1])
def cutoff(m):return math.isqrt(4*(m-1))
def positive(s):
    D,E,U,d=s; return (D-E-d,U-D+E-d,E,d-E)
def diffcone(s):
    D,E,U,d=s;return D>=d>=E>=0 and U>=d

def bank(cap,decorated=False):
    seed=(B,Z0) if decorated else T
    height=lambda x:x[0][2] if decorated else x[2]
    if height(seed)>cap:return {}
    out={seed:(0,'')};heap=[(height(seed),seed)]
    while heap:
        _,s=heapq.heappop(heap);cost,word=out[s]
        for c in 'OZJ':
            q=(step(s[0],c),step(s[1],c)) if decorated else step(s,c)
            need(height(q)>height(s),'strict inventory height')
            if height(q)>cap:continue
            candidate=(cost+(c=='J'),word+c)
            if q not in out:
                out[q]=candidate;heapq.heappush(heap,(height(q),q))
            elif candidate<out[q]:out[q]=candidate
    return out

def graph(word):
    rows=[2,0];root=0
    for c in word:
        n=len(rows);rows.append(0)
        for v in range(n):
            if c!='O' and v==root:rows[n]|=1<<v
            else:rows[v]|=1<<n
        if c=='J':root=n
    return rows

def path_cycle_counts(rows):
    n=len(rows);full=(1<<n)-1
    def count(starts,cycles):
        dp=[{}for _ in range(1<<n)]
        for v in starts:dp[1<<v][v]=1
        for mask in range(1,full+1):
            for last,c in dp[mask].items():
                todo=rows[last]&~mask&full
                while todo:
                    bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                    d=dp[mask|bit];d[j]=d.get(j,0)+c
        return sum(c for last,c in dp[full].items()if not cycles or rows[last]&1)
    return count(range(n),False),count([0],True)

def main(N):
    need(N>=1001,'cap includes both canonical controls')
    basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
    algebra=0
    for v in basis:
        need(K(B,v)==G(v)and K(Z0,v)==H(v),'boundary readouts')
        for w in basis:
            need(K(v,w)==fastK(v,w),'independent K formulas')
            for c in 'OZJ':
                need(K(step(v,c),w)==K(v,step(w,dagger(c))),'adjoint')
                algebra+=1
        for k in range(1,12):
            ck=state('O'*k+'Z');L=2*(2**k-1)
            need(ck==tuple(B[i]+L*Z0[i]for i in range(4)),'physical Mersenne boundary')
            need(K(ck,v)==G(v)+L*H(v)==F(state('O'+'Z'*k,v)),'Mersenne readout')
            algebra+=1
    rays_positive=[(1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1)]
    rays_difference=[(1,0,0,0),(1,0,1,1),(1,1,1,1),(0,0,1,0)]
    for s in rays_positive:
        for c in 'OZJ':need(min(positive(step(s,c)))>=0,'extended positive cone')
    for s in rays_difference:
        for c in 'OZJ':need(diffcone(step(s,c)),'difference cone')
    word_checks=0
    for length in range(8):
        for letters in product('OZJ',repeat=length):
            word=''.join(letters);s=state(word);m,h=G(s),H(s)
            need(3<=h<m<2*h,'joint cone')
            ok=False
            for cut in range(length+1):
                alpha,gamma=word[:cut],word[cut:]
                v=state(alpha);w=state(dagger(gamma),B);q=state(dagger(gamma),Z0)
                need(K(v,w)==m and K(v,q)==h,'joint common-word transfer')
                need(all(0<=q[i]<=w[i]for i in range(4)),'bounded auxiliary')
                if max(v[2],w[2])<=cutoff(m):ok=True
            need(ok,'two-sided cutoff cut exists');word_checks+=1
    # Complete forward inventory: F(s)>=U+1 bounds all targets <=N.
    S=bank(N-1)
    direct_joint={};direct_F={3:0};direct_G={3:0}
    for s,(cost,word)in S.items():
        m,h=G(s),H(s)
        if m<=N:direct_joint[m,h]=min(cost,direct_joint.get((m,h),10**9))
        if m<=N:direct_G[m]=min(cost,direct_G.get(m,10**9))
        n=F(s)
        if n<=N:direct_F[n]=min(cost,direct_F.get(n,10**9))
    C=cutoff(N);small=bank(C);decorated=bank(C,True)
    joint={}; balanced={3:0};families={};witnesses={}
    pair_count=0
    for v,(cv,alpha)in small.items():
        for (w,q),(cw,beta)in decorated.items():
            m=fastK(v,w)
            if m>N or max(v[2],w[2])>cutoff(m):continue
            h=fastK(v,q);cost=cv+cw;pair_count+=1
            need(3<=h<m<2*h,'paired joint cone')
            joint[m,h]=min(cost,joint.get((m,h),10**9))
            balanced[m]=min(cost,balanced.get(m,10**9))
            k=1
            while m+2*(2**k-1)*h<=N:
                n=m+2*(2**k-1)*h
                Ck=math.isqrt(4*(n//2**k-1))
                need(6*2**k<=n+1 and max(v[2],w[2])<=Ck,'target-dependent boundary bounds')
                families[n]=min(cost,families.get(n,10**9))
                if n in [11,23,45,51,227,443,617]and n not in witnesses:
                    word=alpha+dagger(beta)+'O'+'Z'*k
                    need(F(state(word))==n,'compiled boundary recipe')
                    witnesses[n]=dict(m=m,h=h,k=k,v=v,w=w,q=q,alpha=alpha,beta=beta,source_word=word,J_cost=cost)
                k+=1
    need(joint==direct_joint,'full joint inverse with exact minimum J costs')
    need(balanced==direct_G,'balanced inverse with exact minimum J costs')
    combined=dict(balanced)
    for n,cost in families.items():combined[n]=min(cost,combined.get(n,10**9))
    need(combined==direct_F,'full Mersenne criterion with exact minimum J costs')
    # Hostiles use complete source inventory and exact arithmetic.
    v=state('O');w=state('J',B);q=state('O',Z0)
    need((K(v,w),K(v,q))==(15,13),'unsynchronized false positive')
    need((15,13)not in direct_joint,'hostile rejected physically')
    need(sorted(h for m,h in direct_joint if m==15)==[9,11],'complete m15 joint fiber')
    a=(state('OJJZ',B),state('OJJZ',Z0));b=(state('ZZOJ',B),state('ZZOJ',Z0))
    need(a[0]==b[0]==(27,9,27,9)and a[1]!=b[1],'same w distinct q')
    need((F(a[1]),F(b[1]))==(25,39),'distinct continuation values')
    for target in [443,617]:
        need(target in direct_F and target not in direct_G,'nonbalanced hostile target')
        need(not any(m+2*h==target for m,h in direct_joint),'no one-step parent')
    literal=[]
    for target in [11,23,45,51,227,443,617]:
        ww=witnesses[target]['source_word']+'J';rows=graph(ww)
        need(len(rows)<=16,'literal graph count cap')
        hp,hc=path_cycle_counts(rows)
        need((hp,hc)==(target,1),'literal C1 reconstructed graph')
        literal.append(dict(target=target,word=ww,vertices=len(rows),H=hp,C=hc))
    return dict(status='PASS',target_cap=N,complete_forward_U_cap=N-1,complete_source_profiles=len(S),
        paired_U_cap=C,physical_small_profiles=len(small),synchronized_states=len(decorated),
        synchronized_distinct_physical_boundaries=len({w for w,q in decorated}),
        accepted_pair_tests=pair_count,distinct_joint_values=len(joint),represented_C1_targets=len(direct_F),
        balanced_targets=len(direct_G),formal_basis_checks=algebra,common_word_checks=word_checks,
        joint_minimum_J_costs_match=True,C1_minimum_J_costs_match=True,
        budget_counts={r:sum(c<=r for c in direct_F.values())for r in range(4)},
        hostiles=dict(unsynchronized=[15,13],m15_actual_h=[9,11],same_w=a[0],different_q=[a[1],b[1]],no_k1=[443,617]),
        sample_compiled_witnesses=witnesses,literal_graph_checks=literal,
        dependencies={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((Path(__file__).parent/'dependencies').glob('*.md'))},
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        scope='All-order proof candidate; exact finite implementation checks; no universal target coverage.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cap',type=int,default=10001);a=p.parse_args()
    print(json.dumps(main(a.cap),indent=2,sort_keys=True))
