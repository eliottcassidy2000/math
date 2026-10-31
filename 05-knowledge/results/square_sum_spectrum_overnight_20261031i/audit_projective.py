#!/usr/bin/env python3
"""Independent integer projective group arithmetic for the seam curves."""
from fractions import Fraction as Q
from math import isqrt,gcd
import json
from elliptic_ports import setup,frame_at,add
from one_copy_ports import need

def affine(P):
    X,Y,Z=P;return Q(X,Z*Z),Q(Y,Z*Z*Z)

def reduced(P):
    x,y=affine(P);z=isqrt(x.denominator)
    need(z*z==x.denominator,'elliptic coordinate denominator square')
    yy=y*z*z*z;need(yy.denominator==1,'elliptic coordinate denominator cube')
    return x.numerator,int(yy),z

def doubling(P,a,b):
    X,Y,Z=P;W=3*X*X+2*a*X*Z*Z+b*Z**4
    R=W*W-4*a*Y*Y*Z*Z-8*X*Y*Y
    S=W*(4*X*Y*Y-R)-8*Y**4
    return reduced((R,S,2*Y*Z))

def joining(P,Qp,a,b):
    X,Y,Z=P;x,y,z=Qp
    U=X*z*z;u=x*Z*Z;S=Y*z**3;s=y*Z**3
    H=u-U;R=s-S
    need(H!=0,'distinct odd-versus-even multiples')
    xx=R*R-a*H*H*Z*Z*z*z-H*H*(U+u)
    yy=-S*H**3+R*(H*H*U-xx)
    return reduced((xx,yy,H*Z*z))

def main():
    seeds={1:[(3,4,5),(5,12,13),(7,24,25)],2:[(1,7,5),(7,17,13),(17,31,25)],
           5:[(1,2,1),(11,2,5),(1,38,17)],10:[(1,3,1),(9,13,5),(3,19,0)]}
    # Keep only authenticated norm seeds; the last deliberate false seed is rejected.
    rejected=[];valid={}
    for d,triples in seeds.items():
        valid[d]=[]
        for t in triples:
            u,v,w=t
            if gcd(u,v)==1 and u*u+v*v==d*w*w:valid[d].append(t)
            else:rejected.append([d,list(t)])
    # Additional square classes come from primitive norms, then the exact
    # Gaussian rotation creates a second shape in the same class.
    for u in range(1,10):
        for v in range(u+1,12):
            if gcd(u,v)!=1:continue
            n=u*u+v*v;d=n;w=1;p=2
            while p*p<=d:
                while d%(p*p)==0:d//=p*p;w*=p
                p+=1
            if d in valid:continue
            t=2*d*w;x=(t*t-1)*u-2*t*v;y=(t*t-1)*v+2*t*u;z=w*(t*t+1)
            need(x>0 and gcd(x,y)==1 and x*x+y*y==d*z*z,'additional norm seed')
            valid[d]=[(u,v,w),(x,y,z)]
            if len(valid)>=20:break
        if len(valid)>=20:break
    curves=0;point_checks=0
    for d,triples in valid.items():
        for i,(u,v,w) in enumerate(triples):
            for x,y,z in triples[i+1:]:
                left=tuple(sorted((d*u*u,d*v*v)));right=tuple(sorted((d*x*x,d*y*y)))
                _,_,A,B,P,step,clock=setup(left,right,d);a=int(A+B);b=int(A*B)
                point=(int(P[0]),int(P[1]),1);even=doubling(point,a,b)
                need(affine(even)==step,'independent projective doubling')
                rational=P
                for j in range(10):
                    need(affine(point)==rational,'independent projective odd multiples')
                    f=frame_at(rational,left,right,d)
                    L,R=f['left_scale_root'],f['right_scale_root']
                    need(all(L*L*a+R*R*b==H*H for a,b,H in zip(f['left_oriented_port'],f['right_oriented_port'],f['seam_roots'])),'native port equations')
                    point=joining(point,even,a,b);rational=add(rational,step,A,B);point_checks+=1
                curves+=1
    need(curves>=7 and rejected,'positive and hostile norm controls')
    print(json.dumps({'status':'passed','square_classes':sorted(valid),'curves':curves,'independent_odd_multiples':point_checks,
          'rejected_false_norm_seeds':rejected,'scope':'Integer projective arithmetic versus the rational affine producer, with exact seam consequences. This is not a numerical rank computation.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
