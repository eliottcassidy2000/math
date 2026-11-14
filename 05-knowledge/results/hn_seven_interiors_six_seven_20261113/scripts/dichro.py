# Which tournaments cannot be split into two transitive (acyclic) halves?  (dichromatic number >= 3)
import subprocess, itertools
def parse(bits,n):
    adj=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            if bits[k]=='1': adj[i][j]=1
            else: adj[j][i]=1
            k+=1
    return adj
def acyclic(adj,S):
    # a tournament on S is acyclic iff its score sequence is 0..|S|-1
    sc=sorted(sum(adj[i][j] for j in S if j!=i) for i in S)
    return sc==list(range(len(S)))
for n in range(3,9):
    out=subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split()
    bad=[]
    for bits in out:
        if bits.startswith('>'): continue
        adj=parse(bits,n); ok=False
        for mask in range(1<<(n-1)):          # vertex n-1 fixed in part B
            A=[i for i in range(n) if mask>>i&1]; B=[i for i in range(n) if not mask>>i&1]
            if acyclic(adj,A) and acyclic(adj,B): ok=True; break
        if not ok: bad.append(bits)
    print(f"n={n}: tournaments={len([b for b in out if not b.startswith('>')])}, not 2-dicolourable: {len(bad)}")
    if n==7:
        for b in bad:
            adj=parse(b,n)
            # is it a circulant with connection set {1,2,4} up to relabelling? check vertex-transitivity cheaply via scores + 3-cycles
            sc=[sum(r) for r in adj]
            c3=sum(1 for a,b2,c in itertools.combinations(range(n),3) if (adj[a][b2] and adj[b2][c] and adj[c][a]) or (adj[a][c] and adj[c][b2] and adj[b2][a]))
            print("   ",b,"scores",sc,"3-cycles",c3)
