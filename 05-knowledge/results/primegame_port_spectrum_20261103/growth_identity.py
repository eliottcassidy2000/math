"""Growth-only queries recover signed Boolean edge derivatives of H."""
from itertools import combinations,permutations
import json
from pathlib import Path
from port_scout import add,classes,converse,encode,h,port_counts,rows,table
HERE=Path(__file__).resolve().parent

def sink(a):
    n=len(a);return [r|(1<<n) for r in a]+[0]

def edge_set(a,x,y,bit):
    b=a[:]
    if bit:b[x]|=1<<y;b[y]&=~(1<<x)
    else:b[y]|=1<<x;b[x]&=~(1<<y)
    return b

def probe(a,x,y):
    b=sink(a)
    values=[h(add(b,y)),h(add(b,x)),h(add(add(b,x),y)),h(add(add(b,y),x))]
    return values[0]-values[1]-values[2]+values[3],values

def unpadded(a,x,y):
    hs,starts=table(a);_,ends=table(converse(a));full=(1<<len(a))-1
    cxy=sum(ends[s][x]*starts[full^s][y] for s in range(1,full) if s>>x&1 and not s>>y&1)
    cyx=sum(ends[s][y]*starts[full^s][x] for s in range(1,full) if s>>y&1 and not s>>x&1)
    b=port_counts(a)
    direct=h(add(add(a,x),y))-h(add(add(a,y),x))
    return direct,b[y]-b[x]+2*(ends[full][y]-ends[full][x])-cxy+cyx

def main():
    checks=0;examples=[];signature_counts=[]
    for n in range(2,7):
        for mask,a,orbit in classes(n):
            hs,_=table(a);bs=port_counts(a)
            if sum(bs)!=(n-2)*hs[-1]+sum(hs[s]*hs[((1<<n)-1)^s] for s in range(1<<n)):raise RuntimeError('sum rule')
            checks+=1
            for x,y in combinations(range(n),2):
                actual=h(edge_set(a,x,y,1))-h(edge_set(a,x,y,0))
                got,values=probe(a,x,y)
                if actual!=got:raise RuntimeError(('growth derivative',n,mask,x,y,actual,got))
                u,v=unpadded(a,x,y)
                if u!=v:raise RuntimeError('endpoint correction')
                checks+=2
                if (n,mask,x,y) in [(2,0,0,1),(3,0,0,1),(6,16,0,1),(6,16,0,2)]:
                    examples.append({'n':n,'mask':mask,'ports':[x,y],'counts':values,'derivative':got})
        # H and the entire first edge jet: finite injectivity, not presumed.
        if n<=5:
            seen={}
            for mask in range(1<<(n*(n-1)//2)):
                a=rows(n,mask)
                sig=(h(a),)+tuple(h(edge_set(a,x,y,1))-h(edge_set(a,x,y,0)) for x,y in combinations(range(n),2))
                seen.setdefault(sig,[]).append(mask)
            signature_counts.append({'n':n,'classes':len(seen),'cube':1<<(n*(n-1)//2),'max_fibre':max(map(len,seen.values())),'collision':next((v for v in seen.values() if len(v)>1),None)})
    result={'checks':checks,'examples':examples,'first_jet':signature_counts}
    (HERE/'growth_identity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
