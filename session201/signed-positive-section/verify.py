"""Physical signed-section checks. Python stdlib; explicit exceptions."""
from collections import deque
from pathlib import Path
import argparse,hashlib,json
T=(1,0,2,1);ONE=(1,1,1,1);B=(3,1,3,1);Z=(-1,-7,-1,-7)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(word,s=T):
    for c in word:s=step(s,c)
    return s
def W(s):return state('OJJJJ',s)
def signed(s):return 6*s[2]-7*s[0]
def in_slice(s):
    D,E,U,d=s;m=U-D
    return E==2*m and d==3*m

def bank(cap):
    out={T:''};todo=[T]
    for s in todo:
        for c in 'OZJ':
            q=step(s,c)
            need(q[2]>s[2],'strict source growth')
            if q[2]<=cap and q not in out:out[q]=out[s]+c;todo.append(q)
    return out

def graph(word,singleton=False):
    rows=[0]if singleton else[2,0];root=0
    for c in word:
        n=len(rows);rows.append(0)
        for v in range(n):
            if c!='O'and v==root:rows[n]|=1<<v
            else:rows[v]|=1<<n
        if c=='J':root=n
    return rows,root

def literal_profile(rows,root):
    n=len(rows);full=(1<<n)-1
    def path_dp(starts):
        dp=[{}for _ in range(1<<n)]
        for v in starts:dp[1<<v][v]=1
        for mask in range(1,full+1):
            for v,c in dp[mask].items():
                avail=rows[v]&~mask&full
                while avail:
                    bit=avail&-avail;avail-=bit;j=bit.bit_length()-1
                    row=dp[mask|bit];row[j]=row.get(j,0)+c
        return dp
    allpaths=path_dp(range(n));rootpaths=path_dp([root])
    H=[sum(row.values())for row in allpaths];H[0]=1
    P=[sum(row.values())for row in rootpaths]
    U=sum(P[m]*H[full^m]for m in range(1,full+1)if m>>root&1)
    response=(H[full],allpaths[full].get(root,0),U,H[full^(1<<root)])
    cycles=sum(c for v,c in rootpaths[full].items()if rows[v]>>root&1)
    return response,P[full],cycles

def main(cap):
    need(cap>=1000,'test cap includes nontrivial slice fibers')
    need(state('JJJOJ',Z)==ONE,'formal singleton boundary')
    need(state('JJJOJO',Z)==T,'formal F boundary')
    need(state('JJJOJOZ',Z)==B,'formal G boundary')
    formal=0
    for i in range(4):
        s=tuple(int(i==j)for j in range(4));D,E,U,d=s;m=D+2*U
        need(W(s)==(6*m-D,2*m,7*m-D,3*m),'rank two transfer')
        need(state('JZJJJ',s)==W(s),'global suffix rewrite')
        need(signed(state('JZJJJ',s))==D,'unsigned D transfer')
        need(signed(state('ZJZJJJ',s))==D-E+U,'F transfer')
        need(signed(state('OZJZJJJ',s))==D+2*U,'G transfer')
        formal+=5
    S=bank(cap);all_sources=dict(S);all_sources[ONE]='singleton'
    projected={}
    for s,w in all_sources.items():projected.setdefault((s[0],s[2]),(s,w))
    actual={s for s in S if in_slice(s)}
    canonical={W(s)for s in all_sources if W(s)[2]<=cap}
    exceptional={state('OJJ',s)for s in all_sources if s[0]==s[2]and state('OJJ',s)[2]<=cap}
    need(canonical.isdisjoint(exceptional),'two branches disjoint')
    need(actual==canonical|exceptional,'complete intrinsic physical slice')
    samples=[];branch_counts={'canonical':0,'exceptional':0}
    targets=set()
    for q in sorted(actual):
        D,E,U,d=q;m=U-D;n=signed(q);targets.add(n)
        need(n>0 and n%2==m%2==1,'positivity and parity forced')
        if n==m:
            need((n,n)in projected,'balanced exceptional parent')
            s,w=projected[n,n]
            need(state('OJJ',s)==q,'exceptional direct construction')
            need(state('JJ',q)==W(s)and signed(state('JJ',q))==n,'target-preserving JJ rebase')
            need(6*s[2]==U,'sharp exceptional inverse bound')
            branch_counts['exceptional']+=1
        else:
            need(m>=3*n,'no intermediate chamber')
            need((n,(m-n)//2)in projected,'canonical physical parent')
            s,w=projected[n,(m-n)//2]
            need(W(s)==q,'canonical section reconstructs exactly')
            need(14*s[2]<U,'strict canonical inverse bound')
            branch_counts['canonical']+=1
        if len(samples)<8:samples.append(dict(q=q,n=n,m=m,parent=s,parent_word=w))
    corebank=bank(11)
    cores=sorted(s for s in corebank if s[1]>0)
    need(not any(s[0]==7 for s in cores),'all-height D7 core control')
    hostile=state('OZ');need(hostile==(5,1,7,1)and signed(hostile)==7 and not in_slice(hostile),'signed occupancy hostile')
    literal=[]
    for singleton,word in [(True,''),(False,''),(False,'Z'),(False,'OZ'),(False,'ZJ')]:
        seed=ONE if singleton else T;s=state(word,seed);fullword=word+'OJJJJ'
        response,start,cycles=literal_profile(*graph(fullword,singleton))
        need(response==W(s)and start==1 and cycles==1,'literal signed section profile')
        need(signed(response)==s[0],'literal signed target')
        literal.append(dict(singleton_seed=singleton,word=fullword,vertices=len(fullword)+(1 if singleton else 2),profile=response,root_start=start,cycles=cycles,target=s[0]))
    return dict(status='PASS',complete_cover_cap=cap,physical_profiles=len(S),slice_profiles=len(actual),branches=branch_counts,
        distinct_slice_targets=len(targets),formal_basis_checks=formal,core_U11_profiles=len(corebank),core_U11_positive_endpoint_profiles=cores,
        signed7_hostile=hostile,literal_graphs=literal,sample_inverse_witnesses=sorted(samples,key=lambda x:x['q']),
        verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dependency_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((Path(__file__).parent/'dependencies').glob('*.md'))},
        scope='All-order physical-section proof candidate; exact declared finite controls; no all-target occupancy claim.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cap',type=int,default=20000);args=parser.parse_args()
    print(json.dumps(main(args.cap),indent=2,sort_keys=True))
