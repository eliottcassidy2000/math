"""Exact verifier for canonical physical cores. Standard library, no asserts."""
import argparse
import hashlib
import heapq
import itertools
import json
from pathlib import Path

T=(1,0,2,1)
def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def step(s,c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(word,s=T):
    for c in word:s=step(s,c)
    return s
def F(s):return s[0]-s[1]+s[2]
def pair(s):return s[0]+2*s[2],2*s[2]-s[3]
def triple(s):return s[0],s[2],s[3]
def bank(cap):
    """All physical profiles with U<=cap, minimizing (J cost,word)."""
    if cap<2:return {}
    out={T:(0,'')};todo=[(2,T)]
    while todo:
        _,s=heapq.heappop(todo);cost,w=out[s]
        for c in 'OZJ':
            q=step(s,c)
            need(q[2]>s[2],'strict U increase')
            if q[2]>cap:continue
            v=(cost+(c=='J'),w+c)
            if q not in out:
                out[q]=v;heapq.heappush(todo,(q[2],q))
            elif v<out[q]:out[q]=v
    return out

def cores(R,predecessors):
    """Independent core construction via bounded last-letter predecessors."""
    out={(1,2,1):(0,'')}
    for p,(cost,w) in predecessors.items():
        if p[2]>R-2:continue
        for c in 'ZJ':
            q=step(p,c);D,E,u,d=q
            if D+d>R:continue
            key=(D,u,d);rec=(cost+(c=='J'),w+c)
            if key not in out or rec<out[key]:out[key]=rec
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

def literal_counts(rows):
    n=len(rows);full=(1<<n)-1
    def dp(starts,closed):
        table=[{}for _ in range(full+1)]
        for s in starts:table[1<<s][s]=1
        for mask in range(1,full+1):
            for v,c in table[mask].items():
                dest=rows[v]&~mask&full
                while dest:
                    bit=dest&-dest;dest-=bit;u=bit.bit_length()-1
                    nxt=table[mask|bit];nxt[u]=nxt.get(u,0)+c
        return sum(c for v,c in table[full].items()if not closed or rows[v]&1)
    return dp(range(n),False),dp([0],True)

def minput(d,key,cost):d[key]=min(cost,d.get(key,10**9))

def verify(N,L):
    need(N>=1001,'declared cap includes hostile 443 and small slices')
    word_checks=0
    for length in range(L+1):
        for letters in itertools.product('OZJ',repeat=length):
            w=''.join(letters);s=state(w);D,E,U,d=s
            need((E==0)==(U>=2*D),'endpoint-height dichotomy')
            if E:need(D<=U<=2*D-3,'sharp positive endpoint height')
            a=len(w)-len(w.rstrip('O'))
            q=state(w[:-a]if a else w)
            need(q==T or q[1]>0,'stripped word is a core')
            if D==1:
                need(w=='O'*a and U==2**(a+1)and d==1,'transitive exception')
            else:
                b=(U//D).bit_length()-1
                need(a==b and triple(q)==(D,U//2**a,d),'canonical core projection')
            m,h=pair(s);Dc,u,dc=triple(q)
            need(h==2**(a+1)*u-dc and m-h==Dc+dc,'joint ray')
            word_checks+=1
    # Infinite-family identities: tested parameters are finite controls only.
    for j in range(32):
        q=state('O'*j+'Z');r=q[0]+q[3]
        need(state('O'*j)[2]==r-2 and q[2]==2*r-5,'sharp core cutoff family')
    # Complete physical inventory covers all F<=N and all G<=N.
    S=bank(N-1);direct_joint={};direct_F={3:0}
    for s,(cost,w)in S.items():
        D,E,U,d=s
        need((E==0)==(U>=2*D),'full-bank endpoint dichotomy')
        m,h=pair(s)
        if m<=N:minput(direct_joint,(m,h),cost)
        if F(s)<=N:minput(direct_F,F(s),cost)
    R=N//2
    C=cores(R,S)
    need(len(C)<=(R**3-R)//6,'coarse integer core count')
    for (D,u,d),(cost,w)in C.items():
        q=state(w)
        need(triple(q)==(D,u,d),'stored core physical recipe')
        need(q==T or(q[1]>0 and D<=u<=2*D-3),'stored core chamber')
        need(w.count('J')==cost,'stored core cost')
    # Closed integer ray formulas generate pairs without forward continuation.
    ray_joint={};ray_F={3:0};witnesses={};choices=0
    for (D,u,d),(cost,w)in C.items():
        a=0
        while True:
            h=2**(a+1)*u-d;r=D+d;m=h+r
            if m>N:break
            minput(ray_joint,(m,h),cost)
            k=0
            while True:
                n=(2**(k+1)-1)*h+r
                if n>N:break
                minput(ray_F,n,cost);choices+=1
                if n in (5,11,15,23,45,51,227,443,617):
                    word=w+'O'*(a+1)+'Z'*k
                    need(F(state(word))==n,'target recipe reconstructs')
                    rec=dict(core=[D,u,d],a=a,k=k,deficit=r,source_word=word,J_cost=cost)
                    if n not in witnesses or(len(word),word)<(len(witnesses[n]['source_word']),witnesses[n]['source_word']):witnesses[n]=rec
                k+=1
            a+=1
    need(ray_joint==direct_joint,'exact complete joint image and minimum J costs')
    need(ray_F==direct_F,'exact complete C1 image and minimum J costs')
    # Core completeness cross-check and fixed-deficit spectra in a small universe.
    Rsmall=64;Csmall=cores(Rsmall,S)
    direct_cores={triple(s)for s in S if s[0]+s[3]<=Rsmall and(s==T or s[1]>0)}
    need(set(Csmall)==direct_cores,'independent complete small core inventory')
    expected={2:{(1,2,1)},4:{(3,3,1)},6:{(5,5,1),(5,7,1)}}
    for r,expected_set in expected.items():
        need({c for c in C if c[0]+c[2]==r}==expected_set,'small deficit slice')
    need((15,13)not in ray_joint and(11,7)not in ray_joint,'nonphysical hostiles')
    joint_twins=[state('ZZZJJ'),state('ZZJJJ')]
    need(joint_twins==[(27,9,35,17),(29,10,34,15)],'joint-pair collision recipes')
    need(all(pair(s)==(97,53)for s in joint_twins)and
         triple(joint_twins[0])!=triple(joint_twins[1]),
         'joint pair does not determine core projection')
    need(all(s[1]>0 and(s[2]//s[0]).bit_length()-1==0 for s in joint_twins),
         'joint-pair collision consists of actual cores')
    hostile=(3,0,4,1);D,E,U,d=hostile
    need(min(D-E-d,U-D+E-d,E,d-E)>=0,'hostile satisfies relaxed cone')
    need(D%2==d%2==1 and E%2==U%2,'hostile satisfies parity')
    need(hostile not in S and E==0 and U<2*D,'hostile fails physical dichotomy')
    # Verify consequence objects in literal tournaments, separate from profiles.
    literal=[]
    for coreword in ('','Z','ZZ','OZ'):
        q=state(coreword);D,E,u,d=q
        for a in range(2):
            for k in range(2):
                n=(2**(k+1)-1)*(2**(a+1)*u-d)+D+d
                word=coreword+'O'*(a+1)+'Z'*k+'J'
                paths,cycles=literal_counts(graph(word))
                need(paths==n and cycles==1,'literal C1 consequence')
                literal.append(dict(word=word,paths=paths,cycles=cycles))
    # Scalar nonmembership, not D!=U, is the hypothesis for two-ended descent.
    balanced={m for m,h in direct_joint}|{3};descent_checks=0
    scope_hostile=state('ZJ')
    need(scope_hostile==(5,2,6,3) and scope_hostile[0]!=scope_hostile[2],
         'scope hostile has unequal D and U')
    need(F(scope_hostile)==9==pair(state('O'))[0] and 9 in balanced,
         'scope hostile belongs to balanced scalar image')
    v=scope_hostile[2]-scope_hostile[3]
    need((v&-v).bit_length()-1==0,'scope hostile lacks positive outer exponent')
    for s,(cost,w)in S.items():
        n=F(s)
        if n>N or n in balanced:continue
        D,E,U,d=s;v=U-d;k=(v&-v).bit_length()-1;h=v//2**k
        m=n-2*(2**k-1)*h
        parent=(m-h-d,(h+d)//2,d)
        need(k>=1 and(m,h)in ray_joint,'physical canonical outer descent')
        Dp,Up,dp=parent
        if Dp==1:a=Up.bit_length()-2
        else:a=(Up//Dp).bit_length()-1
        core=(Dp,Up//2**a,dp)
        need(core in C and Up%2**a==0,'canonical double descent core')
        need(n==(2**(k+1)-1)*(2**(a+1)*core[1]-dp)+Dp+dp,'canonical descent reconstruction')
        descent_checks+=1
    # Check the actual count bound for every small budget, using target outputs.
    budgets={}
    for R0 in (2,4,6,8,16,32,64):
        vals={3};cnt=0
        for(D,u,d),(cost,w)in C.items():
            if D+d>R0:continue
            cnt+=1;a=0
            while D+2**(a+1)*u<=N:
                h=2**(a+1)*u-d;k=0
                while (2**(k+1)-1)*h+D+d<=N:
                    vals.add((2**(k+1)-1)*h+D+d);k+=1
                a+=1
        need(cnt<=(R0**3-R0)//6,'individual core-count bound')
        need(len(vals)<=1+(R0**3-R0)//6*N.bit_length()**2,'target-count bound')
        budgets[R0]=dict(cores=cnt,represented_targets=len(vals))
    return dict(status='PASS',scope='finite exact implementation checks; proof candidate awaits independent audit',cap=N,word_length=L,word_checks=word_checks,source_profiles=len(S),core_deficit_bound=R,core_triples=len(C),joint_pairs=len(direct_joint),C1_targets=len(direct_F),ray_parameter_choices=choices,canonical_nonbalanced_profiles=descent_checks,budget_counts=budgets,hostile_pairs=[[15,13],[11,7]],witnesses=witnesses,literal_graphs=literal,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cap',type=int,default=5001);p.add_argument('--word-length',type=int,default=8);args=p.parse_args()
    print(json.dumps(verify(args.cap,args.word_length),indent=2,sort_keys=True))
