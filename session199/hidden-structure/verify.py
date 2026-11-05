"""Exact physical bounded-J census and independent small controls (Python 3)."""
from collections import deque, Counter
from itertools import product
from pathlib import Path
import argparse, hashlib, json

HERE=Path(__file__).resolve().parent
T=(1,0,2,1)
BASIS=[tuple(int(i==j) for i in range(4)) for j in range(4)]

def require(test,label):
    if not test:raise RuntimeError(label)

def raw(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)

def read(s):return s[0]+2*s[2]

def replay(s,w):
    for c in w:s=raw(s,c)
    return s

def block(eps):
    k=len(eps)
    if not k:return ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1))
    B=2**k;A=sum(e*2**i for i,e in enumerate(eps));R=sum(e*2**(k-1-i) for i,e in enumerate(eps))
    C=sum(eps[j]*eps[i]*2**(i-j-1) for i in range(k) for j in range(i))+sum(eps[i-1]*eps[i] for i in range(1,k))
    require(C%2==0,'C parity')
    return ((1,-eps[0],A,-C),(0,0,0,eps[-1]),(0,0,B,-R),(0,0,0,1))

def mv(M,s):return tuple(sum(a*b for a,b in zip(row,s)) for row in M)

def census(cap,maxj):
    best={T:0};q=deque([T])
    while q:
        s=q.popleft();cost=best[s]
        for i,c in enumerate('OZJ'):
            cost2=cost+(i==2)
            if cost2>maxj:continue
            z=raw(s,c);g=read(z)
            require(g>read(s),'strict physical readout growth')
            if g>cap:continue
            if cost2<best.get(z,999):
                best[z]=cost2
                (q.append if i==2 else q.appendleft)(z)
    scalar={3:0}
    for s,c in best.items():
        g=read(s)
        if c<scalar.get(g,999):scalar[g]=c
    out={'cap':cap,'max_j':maxj,'source_count':len(best),
         'source_min_j_histogram':dict(sorted(Counter(best.values()).items())),
         'scalar_min_j_histogram':dict(sorted(Counter(scalar.values()).items())),
         'gaps':[n for n in range(3,cap+1,2) if n not in scalar],
         'gaps_at_most_two':[n for n in range(3,cap+1,2) if scalar.get(n,999)>2]}
    return out,scalar

def positive_census(cap,maxj):
    def grow(x,c):
        a,f,e,p=x
        if c=='O':return a+e,a+2*f+e+2*p,0,e+p
        if c=='Z':return 2*a+f+e+p,f,e+p,0
        return f+p,0,a+e+p,e
    def G(x):a,f,e,p=x;return 3*a+2*f+6*e+5*p
    known=set();last=[];values={3:0}
    for cost in range(maxj+1):
        stack=[(0,0,0,1)] if cost==0 else [grow(x,'J') for x in last]
        last=[]
        while stack:
            x=stack.pop()
            if x in known or G(x)>cap:continue
            known.add(x);last.append(x);values.setdefault(G(x),cost)
            stack.extend((grow(x,'O'),grow(x,'Z')))
    return len(known),values

def append(T,c,label):
    n=len(T);out=list(T)+[0]
    for v in range(n):
        if label in ('Z','J') and v==c:out[n]|=1<<v
        else:out[v]|=1<<n
    return out,n if label=='J' else c

def graph(word):
    out=[2,0];root=0
    for c in word:out,root=append(out,root,c)
    return out,root

def counts(out,root):
    n=len(out);full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)];start=[[0]*n for _ in range(full+1)]
    for v in range(n):dp[1<<v][v]=1
    start[1<<root][root]=1
    for m in range(1,full+1):
        for v in range(n):
            if not(m>>v&1):continue
            nxt=out[v]&~m
            while nxt:
                b=nxt&-nxt;w=b.bit_length()-1;nxt-=b
                dp[m|b][w]+=dp[m][v];start[m|b][w]+=start[m][v]
    H=[sum(row) for row in dp];H[0]=1
    U=sum(sum(start[m])*H[full^m] for m in range(full+1))
    profile=H[full],dp[full][root],U,H[full^(1<<root)]
    cyc=sum(start[full][v] for v in range(n) if out[v]>>root&1)
    return profile,sum(start[full]),cyc

def literal_controls():
    total=0
    for k in range(5):
        for ww in product('OZJ',repeat=k):
            w=''.join(ww);out,c=graph(w);s,starts,_=counts(out,c)
            require(s==replay(T,w) and starts==1,'literal source counts')
            # Every reversal edge is exactly the current-root edge at its birth.
            root=0;spine=[0];expected=set()
            for v,letter in enumerate(w,2):
                if letter!='O':expected.add((root,v))
                if letter=='J':spine.append(v);root=v
            actual={(a,b) for b in range(len(out)) for a in range(b) if out[b]>>a&1}
            require(actual==expected,'literal ordered caterpillar reversal edges')
            require(len(spine)==w.count('J')+1,'spine count')
            x,r=append(*append(out,c,'O'),'J');p,starts,cyc=counts(x,r)
            require(p==(read(s),s[0],read(s),s[0]),'balanced literal profile')
            require(starts==cyc==1,'balanced literal source and Hamilton cycle')
            total+=1
    return total

def controls():
    words=0
    for k in range(11):
        for eps in product((0,1),repeat=k):
            M=block(eps);w=''.join('Z' if e else 'O' for e in eps)
            for b in BASIS:require(mv(M,b)==replay(b,w),'binary block on coordinate basis')
            words+=1
    g=[1]*20001
    for n in range(2,len(g)):
        if n%2==0:g[n]=g[n//2]
        else:
            m=(n+1)//2;g[n]=g[m]+m-m%2
        # Independent inverse from u: O at even u; Z at odd u.
        u=n;back=[]
        while u>2:
            if u%2:back.append('Z');u=(u+1)//2
            else:back.append('O');u//=2
        s=replay(T,reversed(back))
        require(s==(g[n],n%2,n,1),'unit source integer inverse')
        require(s==mv(block(tuple(int(c=='Z') for c in reversed(back))),T),'block energy formula')
    for k in range(1,41):
        require(replay(T,'O'*(k-1)+'J')==replay(T,'Z'*k),'unit-deletion prefix rebase')
        if k<=12:
            a,ar=graph('O'*(k-1)+'J');b,br=graph('Z'*k)
            mapping=list(range(1,k+2))+[0]
            require(mapping[ar]==br,'unit-deletion rebase root')
            for i in range(len(a)):
                for j in range(len(a)):
                    require(((a[i]>>j)&1)==((b[mapping[i]]>>mapping[j])&1),'unit-deletion rebase all edges')
    require(replay(T,'OZ')==(5,1,7,1),'ordered-word positive control')
    require(replay(T,'ZO')==(3,0,6,1),'digit-count hostile')
    n,scalar=positive_census(50001,3)
    normal,other=census(50001,3)
    require(n==normal['source_count'] and scalar==other,'independent positive layered census')
    return {'binary_words_through_length_10':words,'coordinate_basis_checks':4*words,
            'unit_sources_through_cover_20000':19999,'unit_deletion_rebases':40,
            'literal_unit_deletion_isomorphisms':12,'literal_words':literal_controls(),
            'independent_census_cap':50001,'independent_source_count':n}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cap',type=int,default=1000001)
    args=parser.parse_args()
    checks=controls()
    data,scalar=census(args.cap,3)
    expected=json.loads((HERE/'GAPS_THREE.json').read_text())
    require(data['gaps']==[n for n in expected if n<=args.cap],'frozen three-J gaps')
    if args.cap==1000001:
        require(data['source_count']==19852699,'frozen source count')
        require(data['scalar_min_j_histogram']=={0:213873,1:266697,2:19137,3:220},'frozen scalar census')
        require(len(data['gaps_at_most_two'])==293 and max(data['gaps_at_most_two'])==879001,'two-J hostile frontier')
    receipt={'status':'PASS','controls':checks,'census':data,
             'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'verify.py',HERE/'PROOF.md',HERE/'GAPS_THREE.json')}}
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':main()
