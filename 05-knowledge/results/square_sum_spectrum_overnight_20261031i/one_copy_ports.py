#!/usr/bin/env python3
"""Universal one-copy join for square-endpoint ports, with an integer clock."""
from fractions import Fraction as Q
from math import gcd,isqrt
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031h'))
from pythagorean_glue import Closed,Atlas,from_cube,scale,edge,need,digest,CounterDegrees

def v2integer(n):
    need(n!=0,'nonzero valuation argument');n=abs(n)
    return (n&-n).bit_length()-1
def v2(q):return v2integer(q.numerator)-v2integer(q.denominator)

def normalize(port):
    a,b=map(isqrt,port)
    need(a*a==port[0] and b*b==port[1] and isqrt(sum(port))**2==sum(port),'square endpoint edge')
    g=gcd(a,b);u,v=a//g,b//g
    if u%2==0:u,v=v,u
    w=isqrt(u*u+v*v)
    need(u%2==1 and v%2==0 and w*w==u*u+v*v and w%2==1,'primitive odd/even frame')
    return g,u,v,w

def initial(u,v,w,x,y,z):
    D=2*w*z
    s=Q(abs(u*u*x*x-v*v*y*y),D)
    a=Q(u*u*x*x+v*v*y*y+2*v*v*x*x,D)
    b=Q(u*u*x*x+v*v*y*y+2*u*u*y*y,D)
    A=(v*x)**2;B=(u*y)**2
    need(a*a==s*s+A and b*b==s*s+B and v2(s)==-1,'initial three-square identity')
    return (s,a,b),A,B

def double(state,A,B):
    s,a,b=state;D=2*s*a*b;x=s*s
    out=(abs((x*x-A*B)/D),(x*x+2*A*x+A*B)/D,(x*x+2*B*x+A*B)/D)
    S,U,V=out
    need(U*U==S*S+A and V*V==S*S+B,'three-square doubling identities')
    need(v2(S)==v2(s)-1,'exact binary denominator clock')
    return out

def frame(left,right,max_left,max_right=None):
    if max_right is None:max_right=max_left
    g,u,v,w=normalize(left);h,x,y,z=normalize(right)
    state,A,B=initial(u,v,w,x,y,z)
    # The first square scaling aligns the two common root factors at gh.
    aligned_max=h*h*max_left
    bound=max(0,(aligned_max.bit_length()+1)//2-1-v2integer(v))
    records=[]
    for j in range(bound+1):
        s,a,b=state;t=s/(u*v);r,k=t.numerator,t.denominator
        need(gcd(r,k)==1 and v2integer(k)==j+1+v2integer(v),'reduced slope clock')
        hyp1=Q(k)*a/v;hyp2=Q(k)*b/u
        need(hyp1.denominator==hyp2.denominator==1,'integral cross-edge hypotenuses')
        records.append({'j':j,'s_v2':v2(s),'denominator_v2':v2integer(k),
                        'numerator_bits':r.bit_length(),'denominator_bits':k.bit_length()})
        physical_left,physical_right=h*r,g*k
        common=gcd(physical_left,physical_right)
        physical_left//=common;physical_right//=common
        if physical_right**2>max_left or physical_left**2>max_right:
            roots=[Q(g*h*int(hyp1),common),Q(g*h*int(hyp2),common)]
            need(all(z.denominator==1 for z in roots),'common square-scale normalization')
            return {'left_normalized':[g,u,v,w],'right_normalized':[h,x,y,z],
                    'A':A,'B':B,'r':r,'k':k,'left_scale_root':physical_left,'right_scale_root':physical_right,
                    'common_scale_root_removed':common,'physical_scale_roots_coprime':gcd(physical_left,physical_right)==1,
                    'left_oriented_port':[(g*u)**2,(g*v)**2],
                    'right_oriented_port':[(h*x)**2,(h*y)**2],
                    'seam_square_roots':[int(z) for z in roots],
                    'aligned_left_maximum':aligned_max,'stage':j,'stage_bound':bound,'clock':records}
        state=double(state,A,B)
    raise RuntimeError('proved denominator stopping bound failed')

def join(A,B):
    A.check_port();B.check_port()
    need(A.square_class==B.square_class==1,'one-copy theorem has square endpoints')
    f=frame(A.port,B.port,max(A.vertices),max(B.vertices));a=f['left_scale_root']**2;b=f['right_scale_root']**2
    va={a*x for x in A.vertices};vb={b*x for x in B.vertices}
    need(not va&vb,'coprime aligned-source disjointness')
    left,right=scale(A.base,a),scale(B.base,b)
    left.remove(edge(a*A.port[0],a*A.port[1]));right.remove(edge(b*B.port[0],b*B.port[1]))
    seams={edge(a*u,b*v) for u,v in zip(f['left_oriented_port'],f['right_oriented_port'])}
    need(sorted(u+v for u,v in seams)==sorted(z*z for z in f['seam_square_roots']),'authenticated arithmetic seams')
    moves=[(scale(R,a),scale(S,a)) for R,S in A.moves]+[(scale(R,b),scale(S,b)) for R,S in B.moves]
    out=Closed(va|vb,left|right|seams,moves,min(seams));out.check_port()
    union=out.base|set().union(*(S for R,S in moves))
    need(max(CounterDegrees(union).values())<=3,'subcubic carrier')
    need({e for e in union if (e[0] in va)!=(e[1] in va)}==seams,'exact two-edge boundary')
    f.update({'source_scales':[a,b],'seams':sorted(seams),'new_port':out.port,
              'source_vertices':[len(A.vertices),len(B.vertices)],'source_bits':[len(A.moves),len(B.moves)],
              'output_vertices':len(out.vertices),'output_bits':len(moves),'maximum_label_bits':max(out.vertices).bit_length()})
    return out,f

def main():
    atlas=Atlas();A=from_cube(atlas.cycle(2812,'A'));B=from_cube(atlas.cycle(2861,'A'),(196,2304))
    C,rec=join(A,B);seen=set();records=[]
    for j in range(64):
        bits=format(j,'06b')[:3]+'0'*22+format(j,'06b')[3:]+'0'*22
        seen.add(digest(C.walk(bits)))
        es=C.edges(bits);a,b=rec['source_scales'];va={a*v for v in A.vertices}
        recoveredA={edge(u//a,v//a) for u,v in es if u in va and v in va}|{A.port}
        recoveredB={edge(u//b,v//b) for u,v in es if u not in va and v not in va}|{B.port}
        need(recoveredA==A.edges(bits[:25]) and recoveredB==B.edges(bits[25:]),'one-copy exact source recovery')
    need(len(seen)==64 and len(C.vertices)==5673 and len(C.moves)==50,'one-copy product clock')
    D,rec2=join(C,A)
    for j in range(4):
        bits=''.join(str(hashlib.sha256(f'{j}:{i}'.encode()).digest()[0]%2) for i in range(len(D.moves)))
        records.append(digest(D.walk(bits)))
    # Same-shape case has A=B in the auxiliary cubic; no nonsingularity is assumed.
    equal=frame((9,16),(36,64),10**12)
    need(equal['A']==equal['B'],'singular cubic control')
    # Normalize an even-first physical port and unequal common root factors.
    flipped=frame((64,225),(196,2304),10**18)
    need(flipped['left_normalized'][1:3]==[15,8],'arithmetic orientation selector')
    # Four more exact recurrence stages audit the algebra without a growth fit.
    st,a,b=initial(3,4,5,7,24,25);valuations=[]
    for j in range(5):
        need(v2(st[0])==-j-1,'all displayed binary-clock stages');valuations.append(v2(st[0]));st=double(st,a,b)
    print(json.dumps({'status':'passed','first_join':rec,'second_join':rec2,'complete_six_bit_cycles':len(seen),
          'second_join_sha256':records,'same_shape_control':equal,'even_first_control':flipped,
          'extra_clock_valuations':valuations,'scope':'All-size one-copy compatibility of authenticated square-endpoint carriers, with exact integer stopping bound. No polynomial bound on integer bit cost is claimed; the companion elliptic_ports.py now handles every common square class.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
