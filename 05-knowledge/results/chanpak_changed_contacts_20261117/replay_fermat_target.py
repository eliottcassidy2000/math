from math import isqrt
from replay_mersenne import need, ordinal, strong, table

SCAFFOLD=(8,1,3,6)
TARGET=7*65537
def converse(rows):
    return tuple(sum((row>>i&1)<<j for j,row in enumerate(rows))
                 for i in range(len(rows)))

def inflate(block):
    labels=[(0,i) for i in range(len(block))]+[(i,0) for i in (1,2,3)]
    def arc(i,j):
        x,u=labels[i]; y,v=labels[j]
        return (block[u]>>v&1) if x==y==0 else (SCAFFOLD[x]>>y&1)
    return tuple(sum(arc(i,j)<<j for j in range(len(labels))) for i in range(len(labels))),len(block)

def profile(rows,v):
    d,H=table(rows)
    b=sum(R[v]*H[(len(H)-1)^S] for S,R in enumerate(d) if S>>v&1)
    c=table(converse(rows))[0][-1][v]
    need(d[-1][v]>0)
    return H[-1],b,c

def join(left,v,right_reverse,u):
    n=len(left); mask=(1<<n)-1
    result=list(left)+[(row<<n)|mask for row in converse(right_reverse)]
    result[v]|=1<<(n+u)
    result[n+u]^=1<<v
    return tuple(result)

for block,want in [((0,),(5,9,1)),((2,0),(11,17,3))]:
    left,v=inflate(block)
    need(profile(left,v)==want)

for block,reader,u,source_pair,right_want in [
    ((2,20,9,19,37,79,31),(126,84,40,18,36,74,28),5,(49,824),(43,106,0)),
    (ordinal((40,1,3,6,15,30),(0,)),(254,160,18,166,170,4,190,36),7,(19,512),(19,208,0)),
]:
    _,counts=table(block)
    h,P=counts[-1],sum(x*counts[(len(counts)-1)^S] for S,x in enumerate(counts))
    need((h,P)==source_pair)
    left,v=inflate(block)
    a,b,c=profile(left,v)
    need((a,b,c)==(3*P-h,4*P+h,P-h))
    d,e,f=profile(reader,u)
    need((d,e,f)==right_want)
    child=join(left,v,reader,u)
    need(a*d+b*e-c*f==TARGET)
    need(strong(child))
    need(table(child)[1][-1]==TARGET)
    print(len(child), h,P,(a,b,c),(d,e,f),TARGET,child)

need(all(65537%d for d in range(2,isqrt(65537)+1)))
need(27*19==512+1 and 127*512+27*19==65537)
need([k for k in range(18) if pow(2,k,27)==26]==[9] and pow(2,18,27)==1)
