from itertools import combinations,product
from math import prod

def need(ok):
    if not ok: raise RuntimeError("profile or mixed readout mismatch")

def get(n,edges):
    dest=[0]*n
    for a,b in edges:dest[a]|=1<<b
    m=(1<<n)-1
    h=[0]*(m+1);h[0]=1
    dp=[[0]*n for _ in h]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,m+1):
        for i,value in enumerate(dp[s]):
            if value:
                z=dest[i]&~s
                while z:
                    b=z&-z;z^=b
                    dp[s|b][b.bit_length()-1]+=value
        h[s]=sum(dp[s])
    def sub(s):
        z=s
        while True:
            yield z
            if not z:return
            z=(z-1)&s
    q=[sum(h[z]*h[s^z] for z in sub(s)) for s in range(m+1)]
    return h[m],q[m],sum(q[z]*h[m^z] for z in sub(m))

def lex(n,c):
    return {e if c>>k&1 else e[::-1] for k,e in enumerate(combinations(range(n),2))}
def row(c):
    return {(i,j) for i in range(8) for j in range(8) if c>>(8*i+j)&1}
def join(*P):
    return tuple(prod(x) for x in zip(*P))

a,b,c,d=map(row,(4584474270165959682,9169081464381770758,9160127196319778056,9169072685987728394))
D=[join(get(8,a),get(8,b)),join(get(8,c),get(8,d))]
E=[join(get(5,lex(5,12)),get(6,lex(6,408))),join(get(5,lex(5,76)),get(6,lex(6,88)))]
need(D==[(1155,1943040,308859075),(1155,1946880,308859075)])
need(E==[(495,43200,1058805),(495,43200,1063125)])
need((len(a^c)+len(b^d))//2==11)
need((76^12).bit_count()+(88^408).bit_count()==4)
for s in [0,1,10,50,120]:
    box={}
    for i,j,k in product(range(2),range(2),range(3)):
        h,x,y=join(D[i],E[j],(1,2**s,3**s))
        box[i,j,k]=4*y+9*x+(-3+2*k)*h
    need(len(set(box.values()))==12)
    need({box[1,j,k]-box[0,j,k] for j,k in product(range(2),range(3))}=={9*43200*3840*2**s})
    need({box[i,1,k]-box[i,0,k] for i,k in product(range(2),range(3))}=={4*308859075*4320*3**s})
    need({box[i,j,k+1]-box[i,j,k] for i,j,k in product(range(2),repeat=3)}=={2*1155*495})
print("Eight independent DP profiles; all twelve physical-operand box positions and three exact strides.")
