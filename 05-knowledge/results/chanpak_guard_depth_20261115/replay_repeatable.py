from itertools import product
from math import gcd

def need(test):
    if not test:
        raise RuntimeError("control failed")

K=(8,1,3,6)
# Independent subset Hamilton path counts for the one strong block.
def literal_counts(rows):
    n=len(rows)
    D={}
    H=[1]+[0]*((1<<n)-1)
    for m in range(1,1<<n):
        for v in range(n):
            if m>>v&1:
                a=m^(1<<v)
                D[m,v]=sum(D[a,u] for u in range(n) if a>>u&1 and rows[u]>>v&1) if a else 1
                H[m]+=D[m,v]
    return H
HC=literal_counts(K)

A=((1,2,1),(1,2,1),(1,2,0),(1,2,1),(1,2,1),(5,32,0))
B=((5,32,1),)+((1,2,0),)*5

def literal_word(word):
    # Enumerate literal paths by their two sides around the apex.
    # Reads old Hamilton path witnesses and the even rooted contribution.
    starts=[]; marked=0; k=0
    for h,p,w in word:
        n=4 if h==5 else 1
        if h==5: starts.append(k)
        if w: marked|=((1<<n)-1)<<k
        k+=n
    H=E=0
    full=(1<<k)-1
    for right in range(1<<k):
        if right and not (right&-right)&marked:
            continue
        left=full^right
        w=1
        for j in starts:
            w*=HC[(left>>j)&15]*HC[(right>>j)&15]
        if not left or not ((marked>>(left.bit_length()-1))&1):
            H+=w
        if right.bit_count()&1:
            E+=w
    return H+2*E

def run(word,D=1,C=0,N=1,z=0):
    for h,p,w in word:
        D,C,N,z=h*D,p*C+w*(p-h)*D,h*N+(2-w)*(p-h)*C+w*(p-h*z)*D,w
    return D,C,N,z

def main():
    need((HC[-1],sum(h*HC[15^i] for i,h in enumerate(HC)))==(5,32))
    for D,C,N in [(1,0,1),(1,1,9),(5,864,1623),(7,33,25)]:
        for w,a,b in [(A,1618,1903),(B,1706,2011)]:
            need(run(w,D,C,N)==(5*D,1024*C+864*D,5*N+b*C+a*D,0))
    for k in range(1,9):
        w=[5**(k-j)*(108*864*(1024**(j-1)-5**(j-1))//1019+88*5**(j-1)) for j in range(1,k+1)]
        need(all(a>204*b for a,b in zip(w[1:],w)))
        if k>=2: need(gcd(*w)==8)
        base=run(A*k)[2]
        outputs=[]
        for e in product((0,1),repeat=k):
            word=sum((B if t else A for t in e),())
            v=run(word)[2]
            need(v==base+sum(t*b for t,b in zip(e,w)))
            if k<=2: need(v==literal_word(word))
            outputs.append(v)
        need(len(set(outputs))==2**k)
        print(k,len(outputs),min(outputs),max(outputs))
    print("PASS; six literal apex partition checks at depths one and two")
if __name__=='__main__':
    main()
