"""Adjacent-cut and diagonal-repair verifier. Stdlib; explicit checks only."""
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
def cutoff(m):return math.isqrt(3*(m-1))
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

def new_checks():
    # Exact matrix checks use raw-coordinate generator images, independent of
    # the displayed positive-coordinate matrices in the proof.
    rays=[(1,0,1,0),(0,0,1,0),(2,1,2,1),(1,0,2,1)]
    matrices={
      'O':[[2,2,3,3],[2,2,0,0],[4,1,3,3],[4,1,0,0]],
      'Z':[[2,2,4,4],[2,2,1,1],[3,0,3,0],[3,0,3,0]],
      'J':[[2,1,4,3],[1,2,2,3],[4,2,2,3],[3,3,3,3]]}
    for c in 'OZJ':
        actual=[[3*K(step(v,c),w)-step(v,c)[2]*step(w,dagger(c))[2]-3*v[3]*w[3]
                 for w in rays]for v in rays]
        need(actual==matrices[c], 'raw-basis adjacent-cut matrix')
        need(all(x>=0 for row in actual for x in row),'nonnegative residual matrix')
    checked=0; boundary_counts={}
    for seed,label in [(T,'full'),(B,'balanced')]:
        boundary_counts[label]=0
        for length in range(9):
            for letters in product('OZJ',repeat=length):
                word=''.join(letters);left=[T];right=[seed]
                for c in word:left.append(step(left[-1],c))
                for c in word[::-1]:right.append(step(right[-1],dagger(c)))
                right=right[::-1];n=K(T,right[0]);heights=[]
                for v,w in zip(left,right):
                    need(K(v,w)==n,'same scalar at every cut')
                    heights.append(max(v[2],w[2]))
                need(min(heights)<=cutoff(n),'sharper cutoff at some cut')
                crossing=next((i for i,(v,w)in enumerate(zip(left,right))if v[2]>=w[2]),None)
                if crossing is None:
                    need(seed==B and word=='' and n==5,'only empty balanced no-crossing case')
                    chosen=0;boundary_counts[label]+=1
                elif crossing==0:
                    need(seed==T and max(left[0][2],right[0][2])==2,'initial equality boundary')
                    chosen=0;boundary_counts[label]+=1
                else:
                    i=crossing;c=word[i-1];v=left[i-1];w=right[i]
                    A=step(v,c)[2];BB=step(w,dagger(c))[2]
                    need(A*BB<=3*(n-v[3]*w[3]),'one-letter physical inequality')
                    need(heights[i-1]==BB and heights[i]==A,'crossed cut orientations')
                    chosen=i-1 if BB<A else i
                need(heights[chosen]<=cutoff(n),'explicit chosen cut meets bound')
                checked+=1
    sharp=[]
    for k in range(21):
        d=2**(k+1)+1;v=state('Z'+'O'*k);w=state('Z'*(k+1)+'J')
        need(w==(2*d-1,d-1,2*d,d),'sharp right source')
        need(v==(3,1 if k==0 else 0,3*2**k,1),'sharp left source boundary')
        n=K(step(v,'O'),w);A=step(v,'O')[2];BB=step(w,'Z')[2]
        need((n,A,BB)==(3*d*d,3*(d-1),3*d),'sharp crossing formulas')
        word='Z'+'O'*(k+1)+'J'+'O'*(k+1)
        hh=[]
        for cut in range(len(word)+1):
            hh.append(max(state(word[:cut])[2],state(dagger(word[cut:]))[2]))
        need(min(hh)==3*(d-1),'exact best cut height sharp family')
        short=word[:-1];hhb=[]
        for cut in range(len(short)+1):
            hhb.append(max(state(short[:cut])[2],state(dagger(short[cut:]),B)[2]))
        need(G(state(short))==n and min(hhb)==3*(d-1),'balanced sharp cut family')
        sharp.append(dict(k=k,target=n,best_cut_height=min(hh),word=word))
    diagonal=[]
    for p in range(11):
        for q in range(1,12):
            x=state('O'*p+'Z'*q+'O'*p)
            need(G(x)==3+2**q*(2**(2*p+2)-1),'palindromic balanced family')
    for a in range(2,32,2):
        p=(a-2)//2;q=a+2;r=2**p;n=(2**(a+1)-1)**2+2
        words=['O'*a+'Z'*a,'O'*p+'Z'*q+'O'*(p+1),'O'*p+'Z'*q+'O'*p+'J','Z'*(p+1)+'O'*q+'Z'*p]
        profiles=list(map(state,words));x=state('O'*p+'Z'*q+'O'*p)
        need(len(set(profiles))==4 and all(F(s)==n for s in profiles),'four distinct diagonal profiles')
        need(G(x)==n,'even diagonal balanced repair')
        need([profiles[i][2]for i in [0,1,3]]==[32*r**4-4*r*r+1,64*r**4-32*r**3+2*r,32*r**4+16*r**3-r+1],'three distinct cover formulas')
        if a<=4:
            graphword='O'*p+'Z'*q+'O'*p+'OJ';rows=graph(graphword)
            need(len(rows)==2*a+4,'one-extra-vertex balanced repair')
            need(path_cycle_counts(rows)==(n,1),'literal even diagonal C1 graph')
            diagonal.append(dict(a=a,target=n,graph_word=graphword,vertices=len(rows),profiles=profiles))
    return dict(status='PASS',raw_matrix_checks=48,word_cut_checks=checked,word_length_cap=8,
                boundary_counts=boundary_counts,sharp_family=sharp,diagonal_literal_graphs=diagonal,
                diagonal_even_indices_checked=list(range(2,32,2)),palindromic_family_checks=121)

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
    full_pair={}
    for v,(cv,alpha)in small.items():
        for w,(cw,beta)in small.items():
            n=fastK(v,w)
            if n<=N and max(v[2],w[2])<=cutoff(n):
                full_pair[n]=min(cv+cw,full_pair.get(n,10**9))
    need(full_pair==direct_F,'smaller full C1 pair inverse preserves minimum J costs')
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
                Ck=math.isqrt(3*(n//2**k-1))
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
    return dict(status='PASS',new_checks=new_checks(),target_cap=N,complete_forward_U_cap=N-1,complete_source_profiles=len(S),
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
        scope='New adjacent-cut cutoff and even diagonal repair, pending independent audit; exact finite checks; no universal coverage or globally optimal cutoff claim.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cap',type=int,default=5001);a=p.parse_args()
    print(json.dumps(main(a.cap),indent=2,sort_keys=True))
