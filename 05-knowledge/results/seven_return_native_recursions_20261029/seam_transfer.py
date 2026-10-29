"""Universal three-coordinate seam response, with exact directed-graph controls."""
from itertools import combinations, permutations, product
from pathlib import Path
import json

def end_counts(adj):
    n=len(adj);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            if not value:continue
            for w in range(n):
                if not(mask>>w)&1 and adj[v][w]:dp[mask|1<<w][w]+=value
    return dp

def response(adj,port):
    n=len(adj);full=(1<<n)-1
    end=end_counts(adj);start=end_counts(list(zip(*adj)))
    H=[sum(row) for row in end];H[0]=1
    R=sum(end[mask][port]*H[full^mask] for mask in range(1,full+1) if mask>>port&1)
    return [H[full],start[full][port],R]

def contracted(adj,entry,exit,state):
    q=2 if state in ('DE','DB') else 1;n=len(adj)
    start=state in ('S','DB','B');end=state in ('E','DE','DB','B')
    out=[[0]*(n+q) for _ in range(n+q)]
    for i in range(n):
        for j in range(n):out[q+i][q+j]=adj[i][j]
    for i in range(n):
        for j in range(q):out[q+i][j]=int(not(i==entry and j==0 and start))
    if end:out[0][q+entry]=1
    return response(out,q+exit)

def matrix(adj,entry,exit,check=True):
    if not adj or not(0<=entry<len(adj) and 0<=exit<len(adj)):raise ValueError('nonempty graph and valid ports required')
    rows={s:contracted(adj,entry,exit,s) for s in ('I','S','E','DE','DB','B')}
    if check and (rows['DB']!=rows['DE'] or any(rows['E'][j]-rows['I'][j]!=rows['DE'][j] for j in range(3))):
        raise ValueError(('three-state reduction',adj,entry,exit,rows))
    if check and any(rows['B'][j]!=rows['E'][j]+rows['S'][j]-rows['I'][j] for j in range(3)):
        raise ValueError('singleton-prefix inclusion-exclusion')
    return [[rows['I'][j],rows['S'][j]-rows['I'][j],rows['DE'][j]] for j in range(3)]

def mv(a,v):return [sum(x*y for x,y in zip(row,v)) for row in a]
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def adjoint(a):
    signs=(1,-1,1)
    return [[signs[i]*a[j][i]*signs[j] for j in range(3)] for i in range(3)]
def tournament(n,mask,spine=False):
    adj=[[0]*n for _ in range(n)];bit=0
    for i,j in combinations(range(n),2):
        if spine and j==i+1:edge=1
        else:edge=mask>>bit&1;bit+=1
        adj[i][j]=edge;adj[j][i]=1-edge
    return adj
def seam(a,exit,b,entry):
    na,nb=len(a),len(b);out=[[0]*(na+nb) for _ in range(na+nb)]
    for i in range(na):
        for j in range(na):out[i][j]=a[i][j]
    for i in range(nb):
        for j in range(nb):out[na+i][na+j]=b[i][j]
        for j in range(na):out[na+i][j]=1
    out[exit][na+entry]=1;out[na+entry][exit]=0
    return out

def literal_response(adj,port):
    n=len(adj);H=S=R=0
    for word in permutations(range(n)):
        bad=[i for i in range(n-1) if not adj[word[i]][word[i+1]]]
        if not bad:
            H+=1;S+=word[0]==port;R+=word[-1]==port
        if len(bad)<=1:
            for cut in (bad if bad else range(n-1)):
                if word[cut]==port:R+=1
    return [H,S,R]

def verify():
    checks=0;blocks=0
    def need(condition,message):
        nonlocal checks
        checks+=1
        if not condition:raise ValueError(message)
    # All labelled tournaments through four, every ordered port pair.
    samples=[]
    for n in range(1,5):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits)
            for entry,exit in product(range(n),repeat=2):
                M=matrix(a,entry,exit);blocks+=1
                need(mv(M,[1,0,0])==response(a,exit),'homogeneous initialization')
                need(matrix(list(zip(*a)),exit,entry)==adjoint(M),'universal mirror adjoint')
                if entry==0 and exit==n-1:samples.append((a,entry,exit,M))
    # More general internal digraphs: loops absent, reciprocal arcs allowed.
    n=3;pairs=list(permutations(range(n),2))
    for bits in range(1<<len(pairs)):
        a=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):a[i][j]=bits>>k&1
        for entry,exit in product(range(n),repeat=2):
            M=matrix(a,entry,exit);blocks+=1
            need(matrix(list(zip(*a)),exit,entry)==adjoint(M),'digraph adjoint')
    # All fixed-spine presentations of sizes five and six.
    for n in (5,6):
        for bits in range(1<<((n-1)*(n-2)//2)):
            a=tournament(n,bits,True);M=matrix(a,0,n-1);blocks+=1
            need(matrix(list(zip(*a)),n-1,0)==adjoint(M),'spine adjoint')
            need(mv(M,[1,0,0])==response(a,n-1),'spine initialization')
    # Independent permutation/cut path counts, actual composite targets.
    small=samples[:10]
    composed=[]
    for (a,ea,xa,Ma),(b,eb,xb,Mb) in product(small,repeat=2):
        c=seam(a,xa,b,eb);v=mv(Mb,response(a,xa))
        need(v==literal_response(c,len(a)+xb),'actual two-block target, literal words and cuts')
        need(matrix(c,ea,len(a)+xb)==mm(Mb,Ma),'actual composite transfer')
        composed.append(v[0])
    witnesses=[tournament(2,0,True),tournament(3,1,True),tournament(3,0,True)]
    vectors=[response(a,len(a)-1) for a in witnesses]
    det=lambda a:a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    need(abs(det(vectors))==2,'native reachable vectors span three dimensions')
    singleton=matrix([[0]],0,0)
    need(singleton==[[1,-1,1],[1,-1,0],[1,0,1]],'singleton tribonacci transfer')
    v=[1,0,0];trib=[]
    for _ in range(14):trib.append(v[0]);v=mv(singleton,v)
    need(all(trib[k]==sum(trib[k-3:k]) for k in range(3,len(trib))), 'inherited tribonacci partner recovered from the universal transfer')
    # The entire seam operator still loses information tested by native complement.
    hostile=[]
    for bits in (19,35):
        a=tournament(5,bits,True)
        b=[[a[i][j] if abs(i-j)==1 else (1-a[i][j] if i!=j else 0)
            for j in range(5)] for i in range(5)]
        M=matrix(a,0,4);v=literal_response(a,4);h=literal_response(b,4)[0]
        need(v==response(a,4),'hostile independent source words/cuts')
        need(h==response(b,4)[0],'hostile independent complement words')
        hostile.append(dict(bits=bits,adjacency=a,matrix=M,response=v,complement_h=h))
    need(hostile[0]['matrix']==hostile[1]['matrix'],'same complete seam transfer')
    need([x['complement_h'] for x in hostile]==[11,13],'complement distinguishes seam-equivalent marked blocks')
    return dict(status='PASS',blocks=blocks,explicit_checks=checks,
                spanning_vectors=vectors,spanning_determinant=det(vectors),
                actual_compositions=len(composed),target_h_min=min(composed),target_h_max=max(composed),
                singleton_matrix=singleton,tribonacci_values=trib,
                native_complement_hostile=hostile)

if __name__=='__main__':print(json.dumps(verify(),indent=2,sort_keys=True))
