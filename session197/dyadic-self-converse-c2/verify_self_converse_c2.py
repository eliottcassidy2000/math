import hashlib
import json
from verify import construct,evaluate,revc,seed_graph,sink,profile,need

def odd_sqrt(a,k):
    need(k>=3 and a%8==1,'odd dyadic square condition')
    x=1
    for j in range(3,k):
        if (x*x-a)%(1<<(j+1)):
            x+=1<<(j-1)
        need((x*x-a)%(1<<(j+1))==0,'square lifting bit')
    return x

def construct_N(V,c):
    n=len(V);full=(1<<n)-1
    left=[row|(1<<n) for row in V]+[1<<c]
    left[c]^=1<<n
    right=[full^(1<<v)^row for v,row in enumerate(V)]
    out=left+[row<<(n+1) for row in right]
    for i in range(n+1):
        for j in range(n):
            if j==c and i in(c,n):out[n+1+j]|=1<<i
            else:out[i]|=1<<(n+1+j)
    anti=[n+1+i for i in range(n)]+[n]+list(range(n))
    return out,anti

def check_anti(T,p):
    n=len(T)
    need(sorted(p)==list(range(n)) and all(p[p[v]]==v for v in range(n)),'anti involution permutation')
    for v in range(n):
        need(not(T[v]>>v&1),'no loop')
        for w in range(v):
            need(((T[v]>>w)&1)+((T[w]>>v)&1)==1,'tournament adjacency')
            need(((T[v]>>w)&1)==((T[p[w]]>>p[v])&1),'every anti edge')

def hp_cycles(T):
    n=len(T);size=1<<n;full=size-1
    def count(cycle):
        dp=[[0]*n for _ in range(size)]
        for v in ([0] if cycle else range(n)):dp[1<<v][v]=1
        for S in range(1,size):
            for v in range(n):
                if not dp[S][v]:continue
                unseen=T[v]&~S
                while unseen:
                    b=unseen&-unseen;unseen-=b;w=b.bit_length()-1
                    dp[S|b][w]+=dp[S][v]
        return sum(dp[full][v] for v in range(n) if not cycle or T[v]&1)
    return count(False),count(True)

def main():
    residue=large=edges=literal=0;digest=hashlib.sha256();examples=[]
    for k in range(3,11):
        M=1<<k
        for target in range(1,M,4):
            t=1 if target%8==5 else 2;U=1<<t
            root=odd_sqrt((target-2*U*(U-1))%M,k)
            D=(root-U)%M
            word=construct(k,1,D,1);n,g=evaluate(word)
            actual_D=g+n;actual_U=(1<<t)*(1+revc(word)+(1<<len(word)))
            value=(actual_D+actual_U)**2+2*actual_U*(actual_U-1)
            need(actual_D%M==D and actual_U%M==U,'actual weak source residues')
            need(value%M==target,'diagonal target residue')
            need(2*(len(word)+2+t)+1<=8*k*k+6*k+9,'physical order bound')
            digest.update(f'{k},{target},{n},{value}\n'.encode());residue+=1
            if k<=4:
                V,c=seed_graph(n+1)
                for _ in range(t):V=sink(V)
                T,anti=construct_N(V,c);check_anti(T,anti);edges+=len(T)*(len(T)-1)//2
                examples.append({'k':k,'target':target,'source_order':len(V),'final_order':len(T),'count':value,'word':word})
    for k in (16,32,64,128):
        M=1<<k
        for target in (1,5,M-3,M-7):
            t=1 if target%8==5 else 2;U=1<<t;x=odd_sqrt((target-2*U*(U-1))%M,k)
            D=(x-U)%M;word=construct(k,1,D,1);n,g=evaluate(word)
            actual_D=g+n;actual_U=(1<<t)*(1+revc(word)+(1<<len(word)))
            need(((actual_D+actual_U)**2+2*actual_U*(actual_U-1))%M==target,'large exact diagonal target')
            large+=1
    # Independent literal final counts from small one-apex sources, not just residue arithmetic.
    for L in range(2,10):
        V,c=seed_graph(L)
        for t in (1,2):
            W=V
            for _ in range(t):W=sink(W)
            D,E,s,U,_,d=profile(W,c)
            need(E==0 and s==d==1,'actual input interface')
            T,anti=construct_N(W,c);check_anti(T,anti)
            H,C=hp_cycles(T)
            need(H==(D+U)**2+2*U*(U-1) and C==2,'literal self-converse C2 count')
            literal+=1
    print(json.dumps({'status':'PASS','all_small_dyadic_targets':residue,'large_exact_targets':large,'literal_anti_edges':edges,'literal_final_H_and_C_counts':literal,'constructed_count_digest':digest.hexdigest(),'small_constructor_examples':examples,'scope':'Exact dyadic residue image conditional on named reviewed structural dependencies; not exact integer coverage.'},sort_keys=True,indent=2))

if __name__=='__main__':main()
