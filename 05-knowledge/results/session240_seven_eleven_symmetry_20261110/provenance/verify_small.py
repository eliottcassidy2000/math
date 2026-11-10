"""Independent finite controls; standard library; explicit errors survive -O."""
import itertools,json
from collections import Counter

def require(x):
    if not x: raise RuntimeError('verification failed')
def rank(rows):
    basis={}
    for x in rows:
        while x:
            i=x.bit_length()-1
            if i in basis:x^=basis[i]
            else:basis[i]=x;break
    return len(basis)
def matmul(a,b):
    return [[sum(x*y for x,y in zip(r,c))%2 for c in zip(*b)] for r in a]
def run():
    out={}
    for p in (7,11):
        qr={i*i%p for i in range(1,p)}
        a=[[int((j-i)%p in qr) for j in range(p)] for i in range(p)]
        masks=[sum(x<<j for j,x in enumerate(r)) for r in a]
        code={0}
        for m in masks:code|={x^m for x in list(code)}
        ext={x|((x.bit_count()%2)<<p) for x in code}
        a2=matmul(a,a)
        aat=matmul(a,list(map(list,zip(*a))))
        identity=[[int(i==j) for j in range(p)] for i in range(p)]
        circle=[(x,y) for x in range(p) for y in range(p) if (x*x+y*y)%p==5]
        require(rank(masks)==(4 if p==7 else 11))
        require((a2==a)==(p==7)); require((aat==identity)==(p==11))
        require(len(circle)==(8 if p==7 else 12))
        out[str(p)]={'qr':sorted(qr),'rank':rank(masks),'projection':a2==a,'inverse_is_transpose':aat==identity,'extended_weight_distribution':dict(sorted(Counter(x.bit_count() for x in ext).items())),'norm5_circle_size':len(circle)}
    words=[w for w in itertools.product((0,1),repeat=5) if all(not(w[i] and w[(i+1)%5]) for i in range(5))]
    labels={w:sum(w[j]*4**(4-j) for j in range(5))%11 for w in words}
    require(len(words)==11 and len(set(labels.values()))==11)
    require(all(labels[w[1:]+w[:1]]==4*labels[w]%11 for w in words))
    affine=[tuple((a*x+b)%11 for x in range(11)) for a in (1,3,4,5,9) for b in range(11)]
    pairs={tuple(sorted((f[0],f[1]))) for f in affine}
    require(len(set(affine))==len(pairs)==55)
    out['golden_period5']={'words':11,'labels':11,'shift_multiplier':4,'affine_group_order':55,'unordered_pair_orbit':55}
    return out
if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
