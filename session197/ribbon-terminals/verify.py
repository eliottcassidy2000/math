"""Exact graph-level verifier for distinct-terminal one-seam transfers."""
from itertools import combinations
import hashlib,json

def require(v,msg):
    if not v:raise RuntimeError(msg)

def all_graphs(n):
    pairs=list(combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        r=[0]*n
        for k,(i,j) in enumerate(pairs):
            if mask>>k&1:r[j]|=1<<i
            else:r[i]|=1<<j
        yield r

def dp_all(r,start=None):
    n=len(r);full=(1<<n)-1;dp=[[0]*n for _ in range(full+1)]
    for i in range(n):
        if start is None or i==start:dp[1<<i][i]=1
    for S in range(1,full+1):
        for i in range(n):
            if not dp[S][i]:continue
            todo=r[i]&(full^S)
            while todo:
                bit=todo&-todo;todo-=bit
                dp[S|bit][bit.bit_length()-1]+=dp[S][i]
    return dp

def profile(r,root):
    dp=dp_all(r);sp=dp_all(r,root);full=len(dp)-1
    H=[sum(x) for x in dp];H[0]=1
    R=sum(H[full^S]*sum(sp[S]) for S in range(1,full+1))
    return (H[full],dp[full][root],R)

def constants(r,b,c):
    dp=dp_all(r);sp=dp_all(r,b);full=len(dp)-1
    H=[sum(x) for x in dp];H[0]=1
    L=sum(dp[S][b]*H[full^S] for S in range(1,full+1))
    D=sum(dp[S][b]*dp[full^S][c] for S in range(1,full))
    return H[full],sum(sp[full]),L,dp[full][c],sp[full][c],D

def augment(r,c):
    n=len(r);return [x|(1<<n if i!=c else 0) for i,x in enumerate(r)]+[1<<c]

def matrix(r,b,c):
    g,t,L,E,F,D=constants(r,b,c)
    gp,tp,Lp,_,_,_=constants(augment(r,c),b,c)
    return [[g,-t,L],[E,-F,D],[gp-g+E,t-tp-F,Lp-L+D]]

def glue(A,a,G,b):
    m=len(A);n=len(G);r=[x|(((1<<n)-1)<<m) for x in A]+[x<<m for x in G]
    r[a]^=1<<(m+b);r[m+b]|=1<<a
    return r

def mv(M,v):return tuple(sum(a*b for a,b in zip(row,v)) for row in M)
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def rm(r,M):return [sum(r[i]*M[i][j] for i in range(3)) for j in range(3)]
def span(r):return max([0]+[j-i for j in range(len(r)) for i in range(j) if r[j]>>i&1])

def wordgraph(word,w=3):
    r=[0]*len(word)
    for j,a in enumerate(word):
        for i in range(j):
            back=i>=j-w and (int(a)>>(i-j+w))&1
            r[j if back else i]|=1<<(i if back else j)
    return r

def main():
    sources=[(r,a,profile(r,a)) for n in range(1,4) for r in all_graphs(n) for a in range(n)]
    cases=0;same=0;distinct=0;digest=hashlib.sha256()
    for n in range(1,5):
        for G in all_graphs(n):
            for b in range(n):
                for c in range(n):
                    M=matrix(G,b,c)
                    for A,a,v in sources:
                        T=glue(A,a,G,b);actual=profile(T,len(A)+c);expected=mv(M,v)
                        require(actual==expected,('transfer',A,a,G,b,c,M,actual,expected))
                        cases+=1;same+=b==c;distinct+=b!=c
                        digest.update(str((A,a,G,b,c,actual)).encode())
    G=[6,4,8,3];b=2;c=3;M=matrix(G,b,c)
    require(constants(G,b,c)==(5,1,7,1,0,1),'base constants')
    require(constants(augment(G,c),b,c)[:3]==(11,2,20),'augmentation constants')
    require(M==[[5,-1,7],[1,0,1],[7,-1,14]],'example matrix')
    ell=[-2,1,1];require(rm(ell,M)==ell,'invariant')
    M2=mm(M,M);M3=mm(M2,M)
    require([M2[0][j]-18*M[0][j]+5*(j==0) for j in range(3)]==[6*x for x in ell],'affine recurrence')
    require(all(M3[i][j]-19*M2[i][j]+23*M[i][j]-5*(i==j)==0 for i in range(3) for j in range(3)),'Cayley Hamilton')
    state=profile(G,c);states=[state];T=G
    physical=[]
    for k in range(1,5):
        if k>1:T=glue(T,len(T)-1,G,b);state=mv(M,state);states.append(state)
        require(T==wordgraph('0013'*k),'physical word')
        actual=profile(T,len(T)-1)
        require(actual==state,('literal clock',k,actual,state))
        require(span(T)<=3,'physical span')
        physical.append(dict(k=k,order=len(T),profile=state,span=span(T)))
    while len(states)<30:states.append(mv(M,states[-1]))
    require(all(-2*h+e+R==-2 for h,e,R in states),'invariant states')
    require(all(states[k+2][0]==18*states[k+1][0]-5*states[k][0]-12 for k in range(len(states)-2)),'recurrence states')
    print(json.dumps(dict(source_root_pairs=len(sources),transfer_checks=cases,same_terminal_checks=same,distinct_terminal_checks=distinct,graph_data_sha256=digest.hexdigest(),matrix=M,physical_clock_checks=physical,first_counts=[v[0] for v in states[:7]],thirtieth_count=states[-1][0]),sort_keys=True,indent=2))

if __name__=='__main__':main()
