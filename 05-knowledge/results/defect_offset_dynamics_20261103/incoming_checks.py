"""Literal port identification with incoming Session164 and canonical continuation."""
from wing import *


def guard_profile(T,q):
    incoming=[i for i in range(len(T)) if T[i]>>q&1]
    require(len(incoming)==1,'unique-in-neighbor guard')
    h,s,_,F=rooted_data(T,q)
    return s,h-s,F-h


def run():
    seed=(8,1,19,6,43,15);count=0
    for left,right in product(range(4),repeat=2):
        T=seed;L=[];R=[]
        for _ in range(left):L.append(len(T));T=birth(T,1<<4)
        for _ in range(right):R.append(len(T));T=birth(T,((1<<len(T))-1)^(1<<0))
        order=[1]+R+[0,3,2,4]+list(reversed(L))+[5]
        expected,_=source2(right+1,left+1)
        require(induced(T,order)==expected,'literal incoming protected-port transport')
        require(order.index(4)==right+4 and order.index(0)==right+1,'fixed left and right guards')
        count+=1
    data=json.loads((HERE/'native_escape_161.json').read_text())
    T=tuple(data['487']['target']);q=2
    require(guard_profile(T,q)[2]==0,'escaped zero defect')
    for command in 'CD':
        new=len(T);T=birth(T,1<<q)
        if command=='C':q=new
    require(guard_profile(T,q)==(487,974,0),'literal canonical entry')
    for command in 'DCD':
        new=len(T);T=birth(T,1<<q)
        if command=='C':q=new
    require(strong(T) and guard_profile(T,q)==(2435,4870,0),'literal canonical multiplier five')
    tables.cache_clear()
    # I2 consists of literal members of the four-module family.
    for m,l in product((1,2,3),repeat=2):
        h0=hwing(m,l,1,1);h1=hwing(m,l,2,1)
        require((h0+1)//6==wing(m,l)[0] and (h1-2*h0-1)//2==wing(m,l)[1],'two actual response probes')
        h2=hwing(m,l,1,2);K0=(h0+1)//2;K1=(h2+1)//2;K2=(h1+1)//2
        require((2*K0-K1,K2-6*(2*K0-K1))==wing(m,l),'three-probe integer left inverse')
    require(__import__('math').gcd(3,5)==1,'coprime maximal minors')
    out=dict(literal_port_identifications=count,canonical_final_h=7305,canonical_final_order=len(T),two_probe_cases=9,integral_probe_smith=(1,1))
    (HERE/'incoming.json').write_text(json.dumps(out,indent=2)+'\n')
    print('incoming literal checks',out,flush=True)
    return out


if __name__=='__main__':run()
