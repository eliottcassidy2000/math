"""All literal two-ear extensions; masks mark old vertices pointing to newborn."""
from dynamics import *


def birth(T, mask):
    n=len(T)
    return tuple(row+((1<<n) if mask>>i&1 else 0) for i,row in enumerate(T))+(((1<<n)-1)^mask,)


def ear_values(T):
    H,E=tables(T);_,S=tables(opposite(T));n=len(T);full=(1<<n)-1
    W=[[0]*n for _ in range(n)]
    for mask in range(1,full):
        left=[(i,x) for i,x in enumerate(E[mask]) if x]
        right=[(j,x) for j,x in enumerate(S[full^mask]) if x]
        for i,x in left:
            for j,y in right:W[i][j]+=x*y
    L=[E[-1][i]-S[-1][i]+sum(W[i]) for i in range(n)]
    Q=[[W[i][j]+W[j][i] for j in range(n)] for i in range(n)]
    vals=[H[-1]]*(1<<n)
    for mask in range(1,1<<n):
        bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
        vals[mask]=vals[rest]+L[i]-sum(Q[i][j] for j in range(n) if rest>>j&1)
    return vals


def all_two(T,targets):
    out=[];n=len(T)
    # Every old Hamilton path admits an insertion of a new vertex. Deletion
    # recovers that path, so H cannot decrease under an induced extension.
    for m,hfirst in enumerate(ear_values(T)):
        if hfirst>max(targets):continue
        U=birth(T,m)
        for q,h in enumerate(ear_values(U)):
            if h in targets:out.append((m,q,h))
        tables.cache_clear()
    return out


def main():
    T,p,v=source(((1,1),))
    T,p,v=exact_vector(T,p,v,(1,5,1,7));p=p[:3]
    T,p=actual(T,p,(0,1,0,1));p=p[:3]
    h=tables(T)[0][-1];print('source',T,h,flush=True)
    hits=all_two(T,{3*h+2,3*h+4})
    print('hits',len(hits),hits[:80],flush=True)
    (HERE/'literal_53.json').write_text(json.dumps(dict(T=T,h=h,hits=hits),indent=2)+'\n')


if __name__=='__main__':main()
