"""Independent rational checks of the interval-conductor audit; stdlib only."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product
import json
from math import lcm

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def floor(x): return x.numerator//x.denominator

def measure(intervals): return sum((b-a for a,b in intervals),Q(0))

def safe(intervals,v,r,phase=Q(0)):
    out=[]
    for a,b in intervals:
        for n in range(floor(v*a+phase)-1,floor(v*b+phase)+2):
            l=max(a,(n+r-phase)/v); h=min(b,(n+1-r-phase)/v)
            if l<=h: out.append((l,h))
    return out

def bad_measure(a,b,v,r,phase=Q(0)):
    total=Q(0)
    for n in range(floor(v*a+phase)-1,floor(v*b+phase)+2):
        l=max(a,(n-r-phase)/v); h=min(b,(n+r-phase)/v)
        if l<h: total+=h-l
    return total

def envelope(x,r):
    m=floor(x); d=x-m
    return (2*r*m+min(d,2*r))/x

def main():
    envelope_checks=literal_measures=conductor_checks=point_checks=0
    for k in range(2,13):
        for ratio in [Q(1),Q(1,2),Q(2,3)]:
            r=ratio/(4*k-2)
            for m in range(0,21):
                for d in [Q(0),r,2*r,3*r,Q(1,2),Q(1)-r]:
                    x=m+d
                    if x>=2*k*r:
                        need(envelope(x,r)<=Q(1,k),'analytic envelope boundary')
                        envelope_checks+=1
    for r in [Q(1,14),Q(1,10),Q(1,6),Q(1,4),Q(1,3)]:
        for a in [Q(0),Q(1,17),Q(-2,7),Q(5,11)]:
            for W in [Q(1,21),Q(2,7),Q(2,3),Q(5,4)]:
                for v in [Q(1),Q(4,3),Q(2),Q(7,2),Q(7)]:
                    for phase in [Q(0),Q(1,5),Q(3,7)]:
                        need(bad_measure(a,a+W,v,r,phase)<=W*envelope(v*W,r),'literal measure exceeds envelope')
                        literal_measures+=1
    for k in range(2,6):
        for r in [Q(1,4*k-2),Q(1,8*k-4)]:
            for u in range(1,5):
                W=2*k*r/u
                for vs in combinations_with_replacement(range(u,u+4),k):
                    for a in [Q(0),Q(1,17),Q(-2,7),Q(5,11)]:
                        I=[(a,a+W)]
                        for j,v in enumerate(vs):
                            phase=Q((2*j+1)%7,7)
                            I=safe(I,Q(v),r,phase)
                        need(bool(I),'boundary conductor failed')
                        for t in [I[0][0],I[-1][1]]:
                            for j,v in enumerate(vs):
                                z=v*t+Q((2*j+1)%7,7); z-=floor(z)
                                need(r<=z<=1-r,'invalid original-coordinate witness')
                        conductor_checks+=1
    # Sharp equality can yield only isolated safe points with shifted speeds.
    for k in range(2,9):
        r=Q(1,4*k-2); W=2*k*r; I=[(-r,(2*k-1)*r)]
        for j in range(k): I=safe(I,Q(1),r,-2*j*r)
        need(bool(I) and measure(I)==0,'point-only equality control')
        point_checks+=1
    # Outside radius condition: k=2,r=1/4,v=3 has two bad pieces over W=1/2.
    need(envelope(Q(3,2),Q(1,4))>Q(1,2),'radius hostile missing')
    # A reduced clock cutoff cannot be justified by the envelope.
    need(envelope(Q(1,7),Q(1,14))>Q(1,2),'cutoff hostile missing')
    general_envelopes=general_intersections=0
    constants=[]
    for k in range(2,10):
        for j in range(1,8):
            r=Q(j,2*k*(j+1)); a=2*r; T=a*(k-1)/(1-k*a)
            N=-floor(-T); C=k*a*N
            need(N-1+a<C<=N,'general last descending branch')
            for n in range(0,N+8):
                for d in [Q(0),a/2,a,(1+a)/2,Q(1)]:
                    x=n+d
                    if x>=C:
                        need(envelope(x,r)<=Q(1,k),'general envelope failed')
                        general_envelopes+=1
            for x in [C,C+Q(1,100),C+Q(1,3)]:
                need(envelope(x,r)<=Q(1,k),'general cutoff boundary')
                general_envelopes+=1
            for u in [1,3]:
                W=C/u
                for start in [Q(0),Q(2,13),Q(-3,7)]:
                    for seq in [tuple(u for _ in range(k)),tuple(u+i%3 for i in range(k))]:
                        I=[(start,start+W)]
                        for i,v in enumerate(seq): I=safe(I,Q(v),r,Q(i,2*k+1))
                        need(bool(I),'general actual intersection')
                        general_intersections+=1
    for k in [5,6]:
        r=Q(1,14); T=2*r*(k-1)/(1-2*k*r); N=-floor(-T)
        constants.append(str(2*k*r*N))
    need(constants==['10/7','30/7'],'five-six tail constants')
    sharp=[]
    for k in range(2,9):
        r=Q(1,4*k-2); c0=Q(99,100)*2*k*r
        centers=[Q(2*j+1,2*k)*c0 for j in range(k)]
        D=lcm(*(z.denominator for z in centers))
        margin=r-c0/(2*k); N=floor(c0/margin)+1; U=N*D
        speeds=[U+int(D*(1-z)) for z in centers]
        W=c0/U; a=Q(1,D); I=[(a,a+W)]
        for v in speeds: I=safe(I,Q(v),r)
        need(not I,'sharpness integer speeds unexpectedly safe')
        need(Q(99,100)*2*k*r<min(speeds)*W<2*k*r,'sharpness normalized clock')
        sharp.append(dict(k=k,r=str(r),speeds=speeds,interval=[str(a),str(a+W)],minimum_speed_times_width=str(min(speeds)*W)))
    return dict(status='PASS',general_envelope_checks=general_envelopes,general_shifted_intersections=general_intersections,five_six_tail_constants=constants,integer_speed_sharpness_controls=sharp,envelope_checks=envelope_checks,literal_measure_checks=literal_measures,
                shifted_boundary_conductor_checks=conductor_checks,point_only_controls=point_checks,
                hostile_controls=2,producer_imports=0)

if __name__=='__main__': print(json.dumps(main(),sort_keys=True))
