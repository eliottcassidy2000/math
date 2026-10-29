"""Independent incidence-to-group reconstruction and lattice normalization."""
from itertools import combinations,product
from collections import Counter
from math import comb
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
def bank(g):
    S={0}
    for k in range(15):S|={v^(g<<k) for v in S}
    return S|{v^((1<<32)-1) for v in S}

def main():
    records=[]
    for g in [73235,83903,92135,102641,118733,129989]:
        octads=[set(i for i in range(32) if v>>i&1) for v in bank(g) if v.bit_count()==8]
        cnt=Counter(tuple(s) for B in octads for s in combinations(sorted(B),4))
        table={(a,31):a for a in range(32)}|{(31,a):a for a in range(32)}|{(a,a):31 for a in range(32)}
        for S,n in cnt.items():
            if n!=7 or 31 not in S:continue
            T=set(S)-{31}
            for a,b in combinations(T,2):
                c=next(iter(T-{a,b}));need((a,b) not in table,'unique sum')
                table[a,b]=table[b,a]=c
        need(len(table)==1024,'total operation')
        for a,b,c in product(range(32),repeat=3):need(table[table[a,b],c]==table[a,table[b,c]],'associative')
        for B in octads:
            origin=min(B);V={table[b,origin] for b in B}
            need(all(table[a,b] in V for a,b in product(V,repeat=2)),'affine3flat')
        rot=lambda a:(a+1)%31 if a!=31 else 31
        need(all(rot(table[a,b])==table[rot(a),rot(b)] for a,b in product(range(32),repeat=2)),'linear cycle')
        masks=[]
        for bits in range(32):
            val=31
            for i in range(5):
                if bits>>i&1:val=table[val,i]
            if val==5:masks.append(bits+32)
        need(len(masks)==1,'minimal polynomial');records.append([g,masks[0]])
    need(records==[[73235,41],[83903,55],[92135,47],[102641,37],[118733,61],[129989,59]],'six charts')
    # Lattice norm4=integer squared length8: either two +/-2 coordinates
    # or one octad with arbitrary +/-1 signs; weight4 parity is absent.
    roots=2*32;norm4=4*comb(32,2)+256*620
    need((roots,norm4)==(64,160704),'root and norm4 shells')
    limit=16
    def conv(A,B):return [sum(A[i]*B[k-i] for i in range(k+1)) for k in range(limit+1)]
    def power(A,n):
        out=[1]+[0]*limit
        for _ in range(n):out=conv(out,A)
        return out
    E=[0]*(limit+1);O=E[:]
    for a in range(-4,5):(E if a%2==0 else O)[a*a]+=1
    weights={0:1,8:620,12:13888,16:36518,20:13888,24:620,32:1}
    theta=[0]*(limit+1)
    for w,c in weights.items():
        term=conv(power(E,32-w),power(O,w))
        theta=[x+c*y for x,y in zip(theta,term)]
    shells=[theta[k] for k in [0,4,8,12,16]]
    need(shells==[1,64,160704,64543488,4845725632],'normalized theta shells')
    need(256*155==39680 and 256*35==8960 and 256*7==1792,'root port factors')
    print(json.dumps({'affine_rotation_polynomials':records,'theta_lattice_norms':[0,2,4,6,8],
        'theta_coefficients':shells,'root_port_factor':256},sort_keys=True))
if __name__=='__main__':main()
