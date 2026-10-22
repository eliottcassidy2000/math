"""Independent all-odd-cycle parity proof audit; no Hamiltonian cancellation."""
from common import core,response
from itertools import combinations,permutations

def odd_cycles(T):
    counts={}
    for r in range(3,len(T)+1,2):
        counts[r]=sum(all(T[p[i]][p[(i+1)%r]] for i in range(r)) for S in combinations(range(len(T)),r) for rest in permutations(S[1:]) for p in [(S[0],)+rest])
    return counts

def check(Q,u,v,fixed):
    k=len(Q);m=k-1;n=k+m;other=[w for w in range(k) if w not in (u,v)]
    R=[v]
    for i in range(m):
        R.append(k+i)
        if i<len(other):R.append(other[i])
    edges=list(zip(R,R[1:]));parity=0;queries=0
    for bits in range(1<<len(edges)):
        word=list(fixed)
        for j,(a,b) in enumerate(edges):
            start,end=(a,b) if bits>>j&1 else (b,a)
            ucore=start if start<k else end;base=end-k if start<k else start-k
            bit=int(start<k);word[base]=(word[base]&~(1<<ucore))|(bit<<ucore)
        word[-1]&=~(1<<u)
        T=[row+[0]*m for row in Q]+[[0]*n for _ in range(m)]
        for i,a in enumerate(word):
            for w in range(k):T[w][k+i]=(a>>w)&1;T[k+i][w]=1-T[w][k+i]
            for j in range(i+1,m):T[k+i][k+j]=1
        cycles=odd_cycles(T);bit=sum(cycles.values())%2
        if bit!=((response(Q,word)-1)//2)%2:raise RuntimeError('mod4 odd-cycle refinement failed')
        parity^=bit;queries+=1
    if parity!=Q[u][v]:raise RuntimeError(('odd-cycle Fourier failed',Q,u,v,parity))
    return queries

if __name__=='__main__':
    tests=queries=0
    for k in [2,3]:
        m=k-1
        for q in range(1<<(k*(k-1)//2)):
            for u,v in permutations(range(k),2):
                # Four deterministic frozen-attachment backgrounds at k3.
                for t in range(1 if k==2 else 4):
                    queries+=check(core(k,q),u,v,[((t+1)*3+i*5)%(1<<k) for i in range(m)]);tests+=1
    for q in [0,2,7,21,63]:
        for u,v in [(0,3),(2,1)]:queries+=check(core(4,q),u,v,[3,10,5]);tests+=1
    print('PASS:',tests,'odd-cycle XOR experiments;',queries,'extensions; all odd-cycle parities agree with Hamiltonian second bit.')
    for q in [41,12]:
        Q=core(5,q);c=odd_cycles(Q);h=response(Q,())
        print('Triangle-only hostile:',{'mask':q,'sorted_scores':sorted(map(sum,Q)),'H':h,'odd_cycle_counts':c,'second_bit':((h-1)//2)%2})
