from itertools import product
EPs=(((0,0),),((1,1),),((0,0),(1,1)),((0,1),(1,0)))
S=[list(product(range(5),repeat=len(ep))) for ep in EPs]
off=[0,5,10,35]
N=60
def mat(bits):
 out=[[] for _ in range(N)]
 for ep,ss,g in zip(EPs,S,off):
  ix={p:i+g for i,p in enumerate(ss)}
  for p,i in ix.items():
   out[i].append(i)
   for j,(u,v) in enumerate(ep):
    r=p[j]
    if r>=3 or (not r and not bits[u]): continue
    q=2 if r==1 else 1
    out[i].append(ix[p[:j]+(q,)+p[j+1:]])
    if not bits[v]:
     out[i].append(ix[p[:j]+(q+2,)+p[j+1:]])
 return out

def go(s,M):
 w=[0]*N
 for i,z in enumerate(s):
  if z:
   for j in M[i]: w[j]+=z
 return bytes(x%3 for x in w)
def back(s,M):
 return bytes(sum(s[j] for j in row)%3 for row in M)
M0=mat((0,1));M1=mat((1,1));MZ=mat((0,0))
init=bytes(int(i in off) for i in range(N))
def accept(arc):
 w=[0]*N
 w[4]=w[9]=2;w[34]=1
 for i,(a,b) in enumerate(S[-1],35):
  if (a>=3 or a==0 and arc) and (b>=3 or b==0 and not arc) and ((a==3)+(b==3))==1:w[i]=2
 return bytes(w)
def calc(s,v):
 return -(1+sum(x*y for x,y in zip(s,v)))%3
def advance(states,length):
    for _ in range(length):
        states={go(s,M) for s in states for M in (M0,M1)}
    return states

def main():
    cur={init}
    seen=set(cur)
    counts=[len(cur)]
    while cur:
        cur=advance(cur,1)-seen
        counts.append(len(cur))
        seen|=cur
    assert counts==[1,2,4,6,10,12,12,10,8,6,2,0]
    assert len(seen)==73

    # Pull the acceptance rows backwards through a positive all-zero suffix.
    # These pairs correspond to odd and even positive suffix lengths.
    v0=tuple(accept(i) for i in (0,1))
    v1=tuple(back(w,MZ) for w in v0)
    v2=tuple(back(w,MZ) for w in v1)
    v3=tuple(back(w,MZ) for w in v2)
    assert v3==v1 and len({v0,v1,v2})==3
    for s in seen:
        t3=advance({s},3)
        t4=advance({s},4)
        for vs in (v1,v2):
            assert {calc(z,v) for z in t3 for v in vs}=={0,1,2}
        assert {calc(z,v2[1]) for z in t4}=={0,1,2}
    print("73 prefix states, suffix period two, 3-bit general and 4-bit strong repair passed.")

if __name__=='__main__': main()
