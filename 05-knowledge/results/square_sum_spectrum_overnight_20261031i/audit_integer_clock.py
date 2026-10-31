#!/usr/bin/env python3
"""Independent homogeneous integer recurrence and hostile source frames."""
from fractions import Fraction
from functools import reduce
from math import gcd,isqrt
import json
from one_copy_ports import initial,double,normalize,frame,v2integer,need

def reduce_state(state):
    g=reduce(gcd,state);return tuple(x//g for x in state)

def integer_step(state,A,B):
    S,U,V,T=state
    q=S**4;tail=A*B*T**4
    out=(abs(q-tail),q+2*A*S*S*T*T+tail,q+2*B*S*S*T*T+tail,2*S*U*V*T)
    result=reduce_state(out)
    s,u,v,t=result
    need(u*u==s*s+A*t*t and v*v==s*s+B*t*t,'integer pencil invariants')
    need(all(z%2 for z in (s,u,v)) and v2integer(t)==v2integer(T)+1,'one paid binary denominator digit')
    return result

def main():
    triples=[]
    for m in range(2,10):
        for n in range(1,m):
            if gcd(m,n)!=1 or (m-n)%2==0:continue
            triples.append((m*m-n*n,2*m*n,m*m+n*n))
    stages=0;frames=0;max_stage=0;max_bound=0;examples=[]
    for i,(u,v,w) in enumerate(triples):
        for j,(x,y,z) in enumerate(triples):
            A=(v*x)**2;B=(u*y)**2
            state=reduce_state((abs(u*u*x*x-v*v*y*y),u*u*x*x+v*v*y*y+2*v*v*x*x,
                                u*u*x*x+v*v*y*y+2*u*u*y*y,2*w*z))
            rational,a,b=initial(u,v,w,x,y,z);need((a,b)==(A,B),'same roots')
            for k in range(3):
                S,U,V,T=state
                need(tuple(Fraction(z,T) for z in (S,U,V))==rational,'independent integer/fraction paths')
                need(v2integer(T)==k+1,'declared initial clock')
                state=integer_step(state,A,B);rational=double(rational,A,B);stages+=1
            # Different common factors and arbitrary endpoint order are retained.
            g=1+i%3;h=1+j%4
            left=tuple(sorted(((g*u)**2,(g*v)**2)))
            right=tuple(sorted(((h*x)**2,(h*y)**2)))
            M=10**(3+(i+j)%6)
            f=frame(left,right,M,M+17)
            l,r=f['left_scale_root'],f['right_scale_root']
            need(gcd(l,r)==1 and (r*r>M or l*l>M+17),'physical primitive separation')
            for a,b,H in zip(f['left_oriented_port'],f['right_oriented_port'],f['seam_square_roots']):
                need(l*l*a+r*r*b==H*H,'actual square connector independent of recurrence')
            need(f['stage']<=f['stage_bound'],'uniform integer stopping bound')
            frames+=1;max_stage=max(max_stage,f['stage']);max_bound=max(max_bound,f['stage_bound'])
            if len(examples)<3 and f['stage']>=2:examples.append({'left':left,'right':right,'M':M,'stage':f['stage'],'bound':f['stage_bound']})
    # The first unnormalized slope for the familiar pair need not separate large intervals.
    need(117**2*400==20**2*13689,'fixed first frame has an actual collision')
    large=frame((9,16),(196,2304),10**12)
    need(large['stage']>0,'changing the integer frame repairs the collision scale')
    print(json.dumps({'status':'passed','primitive_triangles':len(triples),'heterogeneous_frames':frames,
          'independent_integer_stages':stages,'maximum_observed_stage':max_stage,'maximum_stated_bound':max_bound,
          'examples':examples,'fixed_frame_collision':{'scale_roots':[117,20],'source_labels':[400,13689]},
          'large_frame_stage':large['stage'],'scope':'Exact finite controls of two independent arithmetic implementations; uniform termination is proved by the parity invariant, not the observed maximum stage.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
