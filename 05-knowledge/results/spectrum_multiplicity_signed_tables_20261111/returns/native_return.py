"""Literal controls and exact ordered-run certificate for the native +2 return."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
Q=(162,1020,1017,1009,993,960,897,768,513,1)

def ends(rows):
    n=len(rows);M=(1<<n)-1;dp=[[0]*n for _ in range(M+1)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,M+1):
        for i,x in enumerate(dp[mask]):
            if not x:continue
            todo=rows[i]&~mask&M
            while todo:
                bit=todo&-todo;todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=x
    return dp

def converse(rows):
    return [sum(1<<j for j in range(len(rows)) if rows[j]>>i&1) for i in range(len(rows))]

def state(rows,root):
    n=len(rows);M=(1<<n)-1;e=ends(rows);s=ends(converse(rows));h=[sum(x) for x in e];h[0]=1
    U=sum(s[t][root]*h[M^t] for t in range(1,M+1))
    V=sum(e[t][root]*h[M^t] for t in range(1,M+1))
    D=M^(1<<root);K=0;t=D
    while True:
        K+=h[t]*h[D^t]
        if not t:break
        t=(t-1)&D
    return [h[M],e[M][root],s[M][root],U,V,h[D],K]

def birth(rows,pred):
    ans=list(rows);n=len(ans);ans[pred]|=1<<n;ans.append(((1<<n)-1)^(1<<pred));return ans

def inflate(rows,star,k):
    labels=list(range(star))+[star]*k+list(range(star+1,len(rows)));out=[0]*len(labels)
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            yes=labels[i]==labels[j] or rows[labels[i]]>>labels[j]&1
            a,b=(i,j) if yes else (j,i);out[a]|=1<<b
    return out

def strong(rows):
    n=len(rows);reach=[x|(1<<i) for i,x in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(x==(1<<n)-1 for x in reach)

def construct(k):
    parent=inflate(Q,7,k);X=birth(parent,0);x=len(X)-1
    A=birth(X,x);a=len(A)-1;Y=birth(A,a)
    P=birth(Y,a);c=len(P)-1
    zero=birth(P,c);plus=birth(P,a)
    return X,Y,zero,plus,a,c

def ordered_two_cover_runs(rows,star):
    vertices=[i for i in range(len(rows)) if i!=star]
    n=len(vertices);N=1<<n;D=n+3
    dp=[[[0]*D for _ in range(n)] for _ in range(N)]
    for i,v in enumerate(vertices):
        dp[1<<i][i][0]=1;dp[1<<i][i][1]=bool(rows[star]>>v&1)
    for mask in range(1,N):
        for i,u in enumerate(vertices):
            if not(mask>>i&1):continue
            for j,v in enumerate(vertices):
                if mask>>j&1:continue
                direct=rows[u]>>v&1;via=(rows[u]>>star&1) and (rows[star]>>v&1)
                if not(direct or via):continue
                target=dp[mask|1<<j][j]
                for k,x in enumerate(dp[mask][i]):
                    if x:
                        if direct:target[k]+=x
                        if via and k+1<D:target[k+1]+=x
    paths=[[0]*D for _ in range(N)];paths[0][0]=paths[0][1]=1
    for mask in range(1,N):
        for i,u in enumerate(vertices):
            for k,x in enumerate(dp[mask][i]):
                paths[mask][k]+=x
                if rows[u]>>star&1 and k+1<D:paths[mask][k+1]+=x
    coeff=[0]*D
    for mask in range(N):
        for i,x in enumerate(paths[mask]):
            if x:
                for j,y in enumerate(paths[N-1-mask]):
                    if y and i+j<D:coeff[i+j]+=x*y
    return coeff

def main():
    controls=[]
    for k in (1,2,3):
        X,Y,zero,plus,a,c=construct(k);D=1<<k;root=len(plus)-1
        actual=state(plus,root)
        expected=[2673*D+647,1,1188*D+288,183520*3**k+104660*D-7029,
                  2970*D+718,1485*D+359,237530*3**k+134800*D-8982]
        if actual!=expected:raise RuntimeError(('state',k,actual,expected))
        if sum(ends(zero)[-1])!=expected[0]-2:raise RuntimeError(('zero sibling',k))
        flipped=zero.copy()
        for u,v in ((a,root),(c,root)):flipped[u]^=1<<v;flipped[v]^=1<<u
        if flipped!=plus:raise RuntimeError(('native relocation',k))
        if not strong(plus):raise RuntimeError(('strong',k))
        if [i for i in range(len(plus)) if plus[i]>>root&1]!=[a]:raise RuntimeError(('root indegree',k))
        if [i for i in range(len(plus)) if plus[i]>>a&1]!=[len(X)-1]:raise RuntimeError(('predecessor indegree',k))
        oldmask=(1<<len(Y))-1
        if [x&oldmask for x in plus[:len(Y)]]!=Y:raise RuntimeError(('retained old graph',k))
        controls.append({'k':k,'order':len(plus),'root':root,'predecessor':a,'HESUVdK':actual,'rows':plus})
    coefficients=ordered_two_cover_runs(construct(1)[0],7)
    expected=[21832,58612,51170,14390]+[0]*9
    if coefficients!=expected:raise RuntimeError(('run coefficients',coefficients))
    receipt={'scope':'Literal all-seven-coordinate controls at k1,2,3; exact complete ordered-run coefficient certificate for Q_X(2); no coverage scan. Arbitrary k proved in PROOF.md.', 'literal_controls':controls,'Q_X_two_cover_run_coefficients':coefficients,'root_content':1}
    (HERE/'native_return_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'pass','orders':[c['order'] for c in controls],'H':[c['HESUVdK'][0] for c in controls],'Q_X_two_cover_run_coefficients':coefficients}))

if __name__=='__main__':main()
