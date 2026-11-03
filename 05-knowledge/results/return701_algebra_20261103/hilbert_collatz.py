"""Certified inverse Collatz addresses in a recursive Hilbert square chart."""
from functools import lru_cache
from collections import Counter
import hashlib,json
from common import need,save

def exponents(n):
    ks=(2,4,6) if n%3==1 else (3,5,7)
    out=tuple(k for k in ks if (2**k*n-1)%9)
    need(len(out)==2,'two exact unit children')
    return out

def child(n,bit):
    k=exponents(n)[bit];m=(2**k*n-1)//3
    need(m>n and m%2 and m%3,'increasing odd unit inverse child')
    v=m;steps=0
    while v!=n and steps<=k+1:v=3*v+1 if v%2 else v//2;steps+=1
    need(v==n and steps==k+1,'literal native source clock')
    return m,k+1

@lru_cache(None)
def point(depth,index):
    if depth==0:return (0,0)
    side=1<<(depth-1);q,rest=divmod(index,4**(depth-1));x,y=point(depth-1,rest)
    return ((y,x),(x,y+side),(x+side,y+side),(2*side-1-y,side-1-x))[q]

def main():
    sources=[(5,0)];out=[]
    for d in range(1,8):
        for _ in range(2):
            sources=[(m,t+c) for n,t in sources for m,c in (child(n,0),child(n,1))]
        side=1<<d;pts=[point(d,i) for i in range(4**d)]
        need(len(set(n for n,t in sources))==4**d and len(set(pts))==4**d,'exact source and spatial bijections')
        need(all(abs(x-a)+abs(y-b)==1 for (x,y),(a,b) in zip(pts,pts[1:])),'Hilbert scan has literal unit grid steps')
        need(all((x//2,y//2)==point(d-1,i//4) for i,(x,y) in enumerate(pts)),'two inverse steps match nested parent cells')
        # Local binary-basis translations in the final address pair.
        need(all(sum(abs(a-b) for a,b in zip(pts[i],pts[i^1]))==1 for i in range(4**d)),'last-bit XOR is a neighbour swap')
        need(all(sum(abs(a-b) for a,b in zip(pts[i],pts[i^3]))==1 for i in range(4**d)),'other final-pair XOR is a neighbour swap')
        need(all(sum((a-b)**2 for a,b in zip(pts[i],pts[i^2]))==2 for i in range(4**d)),'penultimate-bit XOR is a diagonal swap')
        values=[(i.bit_count()+i//5)%2 for i in range(4**d)]
        D=lambda v,a:[v[i]^v[i^a] for i in range(len(v))]
        need(not any(D(D(values,1),1)) and not any(D(D(values,3),3)),'two local difference automata square to zero')
        need(D(D(values,1),3)==[values[i]^values[i^1]^values[i^2]^values[i^3] for i in range(len(values))],'mixed local difference retains the four-cell joint response')
        raw=json.dumps([[n,t,*p] for (n,t),p in zip(sources,pts)],separators=(',',':')).encode()
        out.append({'depth':d,'inverse_edges_per_source':2*d,'sources':len(sources),'side':side,
          'clock_histogram':dict(sorted(Counter(t for n,t in sources).items())),
          'source_min':min(n for n,t in sources),'source_max':max(n for n,t in sources),
          'full_source_clock_xy_sha256':hashlib.sha256(raw).hexdigest(),'first_eight':[[n,t,*p] for (n,t),p in list(zip(sources,pts))[:8]]})
    save('hilbert_collatz.json',{'status':'FINITE_EXACT_CERTIFIED_SPACE_FILLING_ADDRESS_MODEL','root':5,
       'depths':out,'scope':'ancestry is cell inclusion; scan adjacency is not native Collatz adjacency; universal source membership remains separate'})
    print('PASS source/clock Hilbert charts through16384 certified sources; native clocks and scale-dependent XOR transport',flush=True)
if __name__=='__main__':main()
