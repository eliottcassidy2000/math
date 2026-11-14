from random import Random

def count(blocks, masks):
    n=sum(blocks); out=[0]*(n+1); j=0; a=n
    for w, bits in zip(blocks,masks):
        for i,b in enumerate(bits,j):
            out[i] |= ((1<<n)-1)^((1<<(j+w))-1)
            if w==3:
                out[i] |= 1<<(j+(i-j+1)%3)
            else:
                out[i] |= ((1<<(j+w))-1)^((1<<(i+1))-1)
            if b: out[a] |= 1<<i
            else: out[i] |= 1<<a
        j+=w
    N=1<<(n+1)
    dp=[{} for _ in range(N)]
    for i in range(n+1):dp[1<<i][i]=1
    for s in range(1,N-1):
        for i,v in dp[s].items():
            opts=out[i]&~s
            while opts:
                bit=opts&-opts; opts^=bit; z=bit.bit_length()-1
                D=dp[s|bit]; D[z]=D.get(z,0)+v
    return sum(dp[-1].values())

def main():
    rng=Random(361892)
    cases=0
    for b in [(3,),(3,3),(1,3,3,2),(3,1,3),(2,3,2,3,1,3),(3,3,3,3),(3,1,3,3),(3,1)]:
        for _ in range(15):
            w=[tuple(rng.randrange(2) for _ in range(k)) for k in b]
            tri=[j for j,k in enumerate(b) if k==3]
            sizes=[sum(w[j]) for j in tri]
            r=len(tri)
            n=count(b,w)
            ap=any(i+1==j for i,j in zip(tri,tri[1:]))
            if r>=2:
                t=-sum(sizes[k]*sizes[k+1] for k,(i,j) in enumerate(zip(tri,tri[1:])) if i+1==j)
                d=3**(r-2)
                if n%d or (n//d-t)%3: raise AssertionError((b,w,n,d,t))
            if not ap:
                d=3**(r-1); t=0
                for i,k in zip(tri,sizes):
                    L=w[i-1][-1] if i else 0
                    E=1-w[i+1][0] if i+1<len(w) else 0
                    t+=2*(0<k<3)-k*L+k*E
                if n%d or (n//d-t)%3: raise AssertionError((b,w,n,d,t))
            cases+=1
    print('direct rooted triangle clock controls:',cases,'PASS')

if __name__=='__main__':main()
