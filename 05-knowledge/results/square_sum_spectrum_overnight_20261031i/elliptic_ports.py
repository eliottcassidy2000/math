#!/usr/bin/env python3
"""Same-class one-copy seams via odd multiples and a finite collision bound."""
from fractions import Fraction as Q
from math import gcd,isqrt
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031h'))
from pythagorean_glue import Closed,Atlas,from_cube,condition_cube,glue,scale,edge,need,digest,CounterDegrees
from port_synchronization import square_scale
from one_copy_ports import v2

def sqrtq(q):
    q=Q(q);need(q>=0,'nonnegative square class')
    a,b=isqrt(q.numerator),isqrt(q.denominator)
    need(a*a==q.numerator and b*b==q.denominator,'exact rational square root')
    return Q(a,b)

def normalized(port,d):
    a,b=[isqrt(v//d) for v in port]
    need([d*a*a,d*b*b]==list(port),'declared endpoint square class')
    g=gcd(a,b);u,v=a//g,b//g
    if d%2:
        if u%2==0:u,v=v,u
        need(u%2 and v%2==0,'odd-class primitive orientation')
    else:
        u,v=sorted((u,v));need(u%2 and v%2,'even-class primitive orientation')
    w=isqrt((u*u+v*v)//d)
    need(u*u+v*v==d*w*w and w%2,'primitive norm parity')
    return g,u,v,w

def add(P,R,A,B):
    if P is None:return R
    if R is None:return P
    x,y=P;X,Y=R
    if x==X and y==-Y:return None
    lam=(3*x*x+2*(A+B)*x+A*B)/(2*y) if P==R else (Y-y)/(X-x)
    xx=lam*lam-(A+B)-x-X;yy=-y+lam*(x-xx)
    need(yy*yy==xx*(xx+A)*(xx+B),'exact cubic group law')
    return xx,yy

def setup(left,right,d):
    g,u,v,w=normalized(left,d);h,x,y,z=normalized(right,d)
    A=Q((v*x)**2);B=Q((u*y)**2)
    need(A!=B,'different primitive port shapes')
    P=(Q((u*x)**2),Q(d*u*u*x*x*w*z))
    need(P[1]**2==P[0]*(P[0]+A)*(P[0]+B),'collision point on cubic')
    step=add(P,P,A,B);witness=step if d%2 else add(step,step,A,B)
    need(witness is not None and v2(witness[0])<0,'non-torsion binary witness')
    clock=[]
    for j in range(4):
        clock.append(v2(witness[0]));new=add(witness,witness,A,B)
        need(new is not None and v2(new[0])==v2(witness[0])-2,'infinite distinct doubling orbit')
        witness=new
    return (g,u,v,w),(h,x,y,z),A,B,P,step,clock

def frame_at(point,left,right,d):
    g,u,v,w=normalized(left,d);h,x,y,z=normalized(right,d)
    A=Q((v*x)**2);B=Q((u*y)**2);X,Y=point
    S=sqrtq(X);U=sqrtq((X+A)/d);V=sqrtq((X+B)/d)
    t=S/(u*v);r,k=t.numerator,t.denominator
    L,R=h*r,g*k;common=gcd(L,R);L//=common;R//=common
    roots=[Q(g*h*k*d,common)*U/v,Q(g*h*k*d,common)*V/u]
    need(all(z.denominator==1 for z in roots),'integer seam hypotenuses')
    return {'r':r,'k':k,'left_scale_root':L,'right_scale_root':R,'removed_common_root':common,
            'left_oriented_port':[d*(g*u)**2,d*(g*v)**2],
            'right_oriented_port':[d*(h*x)**2,d*(h*y)**2],'seam_roots':[int(z) for z in roots]}

def collision(frame,left_vertices,right_vertices):
    L,R=frame['left_scale_root'],frame['right_scale_root']
    need(gcd(L,R)==1,'primitive physical scalings')
    if R*R>max(left_vertices) or L*L>max(right_vertices):return None
    for a in sorted(left_vertices):
        if a%(R*R)==0:
            b=L*L*(a//(R*R))
            if b in right_vertices:return a,b
    return None

def find_frame(left,right,d,LV,RV):
    lp,rp,A,B,P,step,clock=setup(left,right,d)
    bound=len(LV)*len(RV)+1;point=P;failures=[]
    for m in range(bound):
        f=frame_at(point,left,right,d);bad=collision(f,LV,RV)
        if bad is None:
            f.update({'odd_multiple':2*m+1,'frames_tested':m+1,'frame_test_bound':bound,
                      'nontorsion_valuation_witness':clock,'rejected_collisions':failures,
                      'square_class':d,'left_normalized':lp,'right_normalized':rp})
            return f
        failures.append({'odd_multiple':2*m+1,'source_labels':bad})
        point=add(point,step,A,B)
        need(point is not None,'non-torsion odd multiple')
    raise RuntimeError('finite collision exhaustion theorem failed')

def join(A,B):
    A.check_port();B.check_port();need(A.square_class==B.square_class,'same square class')
    d=A.square_class;lp=normalized(A.port,d);rp=normalized(B.port,d)
    if lp[1:3]==rp[1:3]:
        # Equal shapes have the earlier direct norm-conic construction.
        left,right=square_scale(A,rp[0]),square_scale(B,lp[0])
        out,r=glue(left,right)
        return out,{'method':'matched norm conic','input_alignment_roots':[rp[0],lp[0]],'frame':r,
                   'output_vertices':len(out.vertices),'output_bits':len(out.moves)}
    f=find_frame(A.port,B.port,d,A.vertices,B.vertices)
    a=f['left_scale_root']**2;b=f['right_scale_root']**2
    LV={a*v for v in A.vertices};RV={b*v for v in B.vertices};need(not LV&RV,'actual source ownership')
    left,right=scale(A.base,a),scale(B.base,b)
    left.remove(edge(a*A.port[0],a*A.port[1]));right.remove(edge(b*B.port[0],b*B.port[1]))
    seams={edge(a*u,b*v) for u,v in zip(f['left_oriented_port'],f['right_oriented_port'])}
    need(sorted(u+v for u,v in seams)==sorted(z*z for z in f['seam_roots']),'actual native seam equations')
    moves=[(scale(r,a),scale(s,a)) for r,s in A.moves]+[(scale(r,b),scale(s,b)) for r,s in B.moves]
    out=Closed(LV|RV,left|right|seams,moves,min(seams),d,A.norm_seed);out.check_port()
    union=out.base|set().union(*(s for r,s in moves))
    need(max(CounterDegrees(union).values())<=3,'subcubic closure')
    need({e for e in union if (e[0] in LV)!=(e[1] in LV)}==seams,'exact source cut')
    f.update({'method':'elliptic odd multiples','source_scales':[a,b],'seams':sorted(seams),
              'source_vertices':[len(A.vertices),len(B.vertices)],'source_bits':[len(A.moves),len(B.moves)],
              'output_vertices':len(out.vertices),'output_bits':len(out.moves),'maximum_label_bits':max(out.vertices).bit_length()})
    return out,f

def main():
    atlas=Atlas();records=[]
    for N,port,d in [(2812,(9,16),1),(2812,(20,80),5),(2828,(578,1922),2),(2812,(160,1440),10)]:
        A,condition=condition_cube(atlas.cycle(N,'A'),port,d)
        if d==1:B=from_cube(atlas.cycle(2861,'A'),(196,2304))
        else:B,_=glue(A,A)
        C,rec=join(A,B);hashes=[]
        for j in range(4):
            bits=''.join(str(hashlib.sha256(f'{j}:{b}'.encode()).digest()[0]%2) for b in range(len(C.moves)))
            hashes.append(digest(C.walk(bits)))
        need(len(C.moves)==len(A.moves)+len(B.moves),'one copy of each branch history')
        records.append({'class':d,'source_condition':condition,'join':rec,'cycle_sha256':hashes})
    # Force eight successive frames to collide; these are finite arithmetic
    # source sets, not asserted to carry Hamiltonian cycles.
    left,right,d=(9,16),(196,2304),1
    _,_,A,B,P,step,_=setup(left,right,d);LV=set(left);RV=set(right);point=P
    for j in range(8):
        f=frame_at(point,left,right,d);LV.add(f['right_scale_root']**2);RV.add(f['left_scale_root']**2)
        point=add(point,step,A,B)
    delayed=find_frame(left,right,d,LV,RV)
    need(delayed['frames_tested']>=9,'hostile delayed escape')
    print(json.dumps({'status':'passed','records':records,
          'delayed_arithmetic_hostile':{'source_sizes':[len(LV),len(RV)],'frames_tested':delayed['frames_tested'],
          'odd_multiple':delayed['odd_multiple'],'bound':delayed['frame_test_bound']},
          'scope':'Universal one-copy compatibility for authenticated invariant ports in one square class. Actual selected-carrier joins in classes1,2,5,10; delayed finite source-set control is arithmetic only.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
