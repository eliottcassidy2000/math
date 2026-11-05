#!/usr/bin/env python3
from verify_unit_offset_rigidity import *
from math import gcd,isqrt
from fractions import Fraction
import hashlib

X=(54,500,152,243,64,340,261,113,157)
Y=(32,17,3,7,13,30)
A0=1998447616;C0=-5304875

def strong(a):
    for v in range(len(a)):
        seen=1<<v
        while True:
            nex=seen
            for u in range(len(a)):
                if seen>>u&1:nex|=a[u]
            if nex==seen:break
            seen=nex
        if seen!=(1<<len(a))-1:return False
    return True

def prime_certificate(q):
    bound=isqrt(q);require(bound*bound<=q<(bound+1)**2,'exact prime trial cutoff')
    digest=hashlib.sha256();minimum=q
    for d in range(2,bound+1):
        r=q%d;require(r!=0,('proper divisor found',q,d));minimum=min(minimum,r)
        digest.update(f'{d}:{r}\n'.encode())
    return {'prime':q,'floor_sqrt':bound,'integers_tested':bound-1,'minimum_nonzero_remainder':minimum,'all_remainders_sha256':digest.hexdigest()}

def cross_control(a,m):
    root=0;outside=list(range(m,len(a)))
    O=[v for v in outside if a[root]>>v&1];I=[v for v in outside if a[v]>>root&1]
    edges=[(u,v) for u in O for v in I if a[u]>>v&1]
    # P is at m, and the fixed singleton immediately after X is at m+10.
    cover={m,m+10}
    require(all(u in cover or v in cover for u,v in edges),'uniform physical root crossing cover')
    require(any(u!=w and v!=z for u,v in edges for w,z in edges),'physical rank2 matching')
    return {'O':O,'I':I,'edges':edges,'cover':sorted(cover),'matching_rank':2}

def run():
    validate(X);validate(Y)
    px=bin_polynomial(X);py=bin_polynomial(Y)
    require(px==[0,104,0,648,0,354,0,40,0,1],'entire X profile')
    require(py==[0,0,24,0,12,0,1],'entire Y profile')
    hx,bx,_=profile(X);hy,by,_=profile(Y)
    require((hx,bx,hy,by)==(1147,22352,37,352),'bank responses')
    require((hx,bx,hy,by)==(31*37,2**4*11*127,37,2**5*11),'exact response factors')
    require(coeff((0,),X,transitive(7),Y)==(A0,0,C0),'physical arbitrary-module coefficients')
    require(A0==2**10*11**2*127**2 and C0==-5**3*31*37**2,'base coefficient factors')
    allocations=0;hits=0
    for uk,ux,uy,vk,vx,vy in product(range(4),repeat=6):
        if uk+ux+uy==0 or vk+vx+vy==0:continue
        hu=1147**ux*37**uy;bu=2**uk*22352**ux*352**uy
        hv=1147**vx*37**vy;bv=2**vk*22352**vx*352**vy
        B=254*bv*hu-124*bu*hv
        expected=uk==vk and ux==vx+1 and vy==uy+1
        require((B==0)==expected,'complete nonnegative allocation classification control')
        allocations+=1;hits+=expected
    arithmetic=0
    for t,s,j in product(range(9),repeat=3):
        A=A0*4**t*(22352**2)**s*(352**2)**j
        C=C0*(1147**2)**s*(37**2)**j
        require(gcd(A,C)==1 and gcd(3*A+C,9*A+C)==1,'primitive entire ray')
        for m in range(1,19):
            H=A*3**m+C
            require(H>0 and H%2==1,'physical sign and parity')
            require((H%7==0)==((m+4*t-4)%6==0),'exact7 targeting rule')
            arithmetic+=1
    witnesses=[]
    for m in (1,4):
        a=graph(transitive(m),(0,),X,transitive(7),Y)
        require(len(a)==m+25 and strong(a),'large physical witness order/strongness')
        H=A0*3**m+C0
        q=H if m==1 else H//7
        require(H==5990037973 if m==1 else H==161868952021 and H==7*23124136003,'exact target values')
        witnesses.append({'m':m,'n':len(a),'rows':a,'H_from_proved_identity':H,'type':'prime' if m==1 else '7p','prime_certificate':prime_certificate(q),'crossing_control':cross_control(a,m)})
    print(json.dumps({'status':'PASS','bank_X':{'rows':X,'bin_polynomial':px,'H':hx,'B2':bx},'bank_Y':{'rows':Y,'bin_polynomial':py,'H':hy,'B2':by},'base_coefficients':[A0,0,C0],'nonempty_allocation_cases':allocations,'cancelling_allocations':hits,'arithmetic_controls':arithmetic,'witnesses':witnesses},sort_keys=True,indent=2))
if __name__=='__main__':run()
