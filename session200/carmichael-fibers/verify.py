"""Exact source-fiber and canonical-descent verifier; Python stdlib only."""
from collections import defaultdict
from itertools import product, permutations
from pathlib import Path
import argparse, hashlib, json

T=(1,0,2,1)
def need(x,msg):
    if not x: raise RuntimeError(msg)
def step(s,c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)
def state(w,s=T):
    for c in w:s=step(s,c)
    return s
def F(s):return s[0]-s[1]+s[2]
def G(s):return s[0]+2*s[2]
def dagger(w):return ''.join({'O':'Z','Z':'O','J':'J'}[c] for c in w[::-1])
def v2(x):
    need(x>0,'positive valuation input')
    return (x&-x).bit_length()-1

def bank(cap):
    out={T:''};todo=[T]
    for s in todo:
        for c in 'OZJ':
            q=step(s,c)
            need(q[2]>s[2],'strict physical height')
            if q[2]<=cap and q not in out:
                out[q]=out[s]+c;todo.append(q)
    return out

def graph(w):
    rows=[2,0];root=0
    for c in w:
        n=len(rows);rows.append(0)
        for v in range(n):
            if c!='O' and v==root:rows[n]|=1<<v
            else:rows[v]|=1<<n
        if c=='J':root=n
    return rows,root

def counts(w):
    rows,_=graph(w);n=len(rows);H=C=0
    for p in permutations(range(n)):
        if all(rows[p[i]]>>p[i+1]&1 for i in range(n-1)):
            H+=1
            if p[0]==0 and rows[p[-1]]&1:C+=1
    return H,C

def isomorphic(A,B):
    n=len(A)
    return any(all(((A[i]>>j)&1)==((B[p[i]]>>p[j])&1)
                   for i in range(n) for j in range(n)) for p in permutations(range(n)))

def main(N):
    algebra=0
    for s in [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]:
        D,E,U,d=s
        for k in range(9):
            need(F(state('O'+'Z'*k,s))==G(s)+(2**(k+1)-2)*(2*U-d),'O split identity')
            need(F(state('J'+'Z'*k,s))==D+2**(k+1)*U,'J split identity')
            algebra+=2
    words=0
    for length in range(7):
        for letters in product('OZJ',repeat=length):
            w=''.join(letters);s=state(w);D,E,U,d=s
            need(D%2==d%2==1 and (E-U)%2==0,'physical parity')
            need(D-E>=d>=1 and U>=D and U>=2*d,'physical cone')
            need(F(step(s,'O'))==F(step(s,'J'))==G(s),'physical twin')
            need(step(s,'O')!=step(s,'J'),'distinct profiles')
            need(F(s)==F(state(dagger(w))),'physical dagger')
            h=2*U-d;m=G(s)
            need(3<=h<m<2*h,'strict parent cone')
            for k in range(1,5):
                q=state('O'+'Z'*k,s)
                need(v2(q[2]-q[3])==k,'canonical valuation')
                need(F(q)==m+2*(2**k-1)*h,'parent target')
                need(m<F(q)/2,'strict scalar descent')
            words+=1
    for a in range(1,16):
        for b in range(1,16):
            s=state('O'*a+'Z'*b)
            expected=(2+(2**b-1)*(2**(a+1)-1),1,1+2**b*(2**(a+1)-1),1)
            need(s==expected,'two-run profile')
            need(F(s)==2+(2**(a+1)-1)*(2**(b+1)-1),'Mersenne product')
            need((s==state('O'*b+'Z'*a))==(a==b),'off-diagonal twins')
    literal=[]
    for w in ['', 'O', 'Z', 'J', 'OZ']:
        for c in 'OJ':
            for k in range(2):
                q=w+c+'Z'*k;actual=counts(q+'J')
                need(actual==(F(state(q)),1),'literal C1 graph')
                literal.append(dict(word=q+'J',H=actual[0],C=actual[1]))
    need(isomorphic(graph('OJ')[0],graph('JJ')[0]),'H5 graph collision')
    S=bank(N-1);fibers=defaultdict(list);balanced={3}
    for s,w in S.items():
        if F(s)<=N:fibers[F(s)].append((s,w))
        if G(s)<=N:balanced.add(G(s))
    singleton={n:dict(profile=ss[0][0],word=ss[0][1])for n,ss in sorted(fibers.items())if len(ss)==1}
    need(set(singleton)=={3,11},'bounded singleton census')
    nonbalanced=[]
    for n in sorted(set(fibers)-balanced):
        descents=[]
        for q,w in fibers[n]:
            k=v2(q[2]-q[3]);h=(q[2]-q[3])//2**k;m=n-2*(2**k-1)*h
            need(k>=1 and w.endswith('O'+'Z'*k),'canonical recipe suffix')
            s=state(w[:-(k+1)])
            need((G(s),2*s[2]-s[3])==(m,h),'joint physical parent')
            need(m in balanced and 3<=h<m<2*h and 2*m<n,'balanced strict descent')
            descents.append(dict(profile=q,word=w,k=k,h=h,m=m,parent_word=w[:-(k+1)]))
        nonbalanced.append(dict(target=n,minimum_k=min(x['k']for x in descents),descents=descents))
    # An independent direct search for the one-step parent uses n>=5U+1.
    one_step=[]
    for s,w in bank((443-1)//5).items():
        if F(state('OZ',s))==443:one_step.append(w)
    need(not one_step,'443 no physical one-step parent')
    need({s for s,w in fibers[443]}=={(219,1,225,1),(191,1,253,1)},'443 complete fiber')
    missing_one=[x['target']for x in nonbalanced if x['minimum_k']>1]
    need(missing_one==[443,617],'bounded one-step exceptions')
    # Independent all-k family generation from the parent bank, including all-Z boundary.
    family=set(balanced)
    for s in S:
        m=G(s);h=2*s[2]-s[3];k=1
        while m+2*(2**k-1)*h<=N:
            family.add(m+2*(2**k-1)*h);k+=1
    need(family==set(fibers),'complete target compiler')
    return dict(status='PASS',target_cap=N,complete_U_cap=N-1,physical_profiles=len(S),
                represented_targets=len(fibers),balanced_targets=len(balanced),singleton_fibers=singleton,
                formal_basis_checks=algebra,source_word_checks=words,two_run_checks=225,literal_graphs=literal,
                H5_profiles_close_to_isomorphic_graphs=True,nonbalanced=nonbalanced,
                one_step_exceptions=missing_one,all_k_compiler_matches_complete_C1=True,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Exact stated finite universe; all-order twins/splitting/descent proved in PROOF.md. No universal coverage or singleton theorem.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cap',type=int,default=10001);a=p.parse_args()
    need(a.cap>=1001,'verification cap must include 617')
    print(json.dumps(main(a.cap),indent=2,sort_keys=True))
