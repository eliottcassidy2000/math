#!/usr/bin/env python3
"""Exact rational identities and complete marked F2^4 deck-group control."""
from fractions import Fraction as F
from itertools import combinations
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def span(gens):
    out={0}
    for g in gens:out|={v^g for v in tuple(out)}
    return out
def dot(a,b):return (a&b).bit_count()%2
def code(labels):return {sum(dot(a,b)<<i for i,b in enumerate(labels)) for a in range(8)}
def proj(z):return (z>>3)&15
def trim(p):
    p=list(map(F,p))
    while p and p[-1]==0:p.pop()
    return p
def gcdpoly(a,b):
    a,b=trim(a),trim(b)
    while b:
        while len(a)>=len(b):
            shift=len(a)-len(b);factor=a[-1]/b[-1]
            for i,c in enumerate(b):a[i+shift]-=factor*c
            a=trim(a)
        a,b=b,a
    return [x/a[-1] for x in a]
def derivative(a):return [i*a[i] for i in range(1,len(a))]

def main():
    U=code((1,2,4,5,6,3,7));D=code((1,2,4,6,3,7,5))
    CU={proj(z) for z in U if not z&2};CD={proj(z) for z in D if not z&2}
    C=span(CU|CD)
    need(C=={z for z in range(16) if dot(z,7)==0},'unchanged quotient ports have kernel 7')
    need(span((3,5,8))==C,'u,v,y character generators')
    tests=0
    for den in range(1,18):
        for num in range(-40,41):
            s=F(num,den);B=1+4*s**4;A=1-2*s*s+4*s**4
            t=4*s*s/B;u=(1-4*s**4)/B;v=2*s*(1+2*s*s)/B
            need(u*u==1-t*t and v*v==t*(1+t),'rational parametrization')
            need(v/(1+t+u)==s,'marked inverse')
            need(2-t==2*A/B and (1+t)*B==(1+2*s*s)**2,'residual squareclass')
            tests+=1
    A=[1,0,-2,0,4];B=[1,0,0,0,4]
    need(gcdpoly(A,derivative(A))==[1] and gcdpoly(B,derivative(B))==[1]
         and gcdpoly(A,B)==[1],'squarefree degree-eight hyperelliptic polynomial')
    delta=F(5,6);t=F(4,5)
    need(1-t*t==F(9,25) and t*(1+t)==F(36,25) and delta*(2-t)==1,'rational quotient pair')
    for z in U|D:
        if not z&2:need(dot(proj(z),7)==0,'each surviving original character rational')
    ports=(1,2,4,8,15);pairs=list(combinations(range(5),2))
    label={p:ports[p[0]]^ports[p[1]] for p in pairs}
    need(len(set(label.values()))==10 and set(label.values()).isdisjoint(ports),'ten distinct free involutions')
    zero_subsets=[I for k in range(1,6) for I in combinations(ports,k) if len(span(I))<2**k]
    need(zero_subsets==[ports],'unique five-port circuit')
    graph={v:{v^a for a in ports} for v in range(16)}
    need(all(len(graph[a]&graph[b])<=2 for a,b in combinations(range(16),2)),'no K5 or K3,3 subgraph')
    petersen_edges=0;free_subgroups=set();genus_counts={1:0,3:0}
    for h in range(1,16):
        fixed=8 if h in ports else 0
        # 2g(C)-2=8=2(2g(Q)-2)+fixed.
        genus=(12-fixed)//4;genus_counts[genus]+=1
    for p,q in combinations(pairs,2):
        disjoint=not set(p)&set(q);a,b=label[p],label[q]
        need(disjoint==((a^b) in ports),'Petersen adjacency equals ramified join')
        petersen_edges+=disjoint
        H=frozenset(span((a,b)));inertia=len(set(ports)&H)
        need(inertia==int(disjoint),'complete inertia intersection')
        # 8=4(2g(Q)-2)+8*inertia.
        need((16-8*inertia)//8==(1 if disjoint else 2),'V4 quotient genus')
        if not inertia:free_subgroups.add(H)
    need(petersen_edges==15 and len(free_subgroups)==10,'Petersen and free V4 counts')
    star=[label[(0,j)] for j in (1,2,3)]
    triangle=[label[p] for p in ((0,1),(0,2),(1,2))]
    need(all(not(set(ports)&span(pair)) for pair in combinations(star,2)),'star pairwise free')
    need(set(ports)&span(star)=={ports[4]},'star fails simultaneous freeness')
    need(not set(ports)&span(triangle) and len(span(triangle))==4,'triangle free V4')
    need(label[(3,4)]==7,'marked missing involution and physical branch ports')
    need({z for z in range(16) if all(dot(z,h)==0 for h in star)}=={0,15},'elliptic star quotient character')
    need(span(C|{15})==set(range(16)) and 15 not in C,'elliptic coordinate completes original field')
    need(F(5,6)*F(4,5)*(1-F(4,5)**2)*(2-F(4,5))==F(36,125),'elliptic missing-root hostile')
    # Algebraic reconstruction on a squareclass-trivial rational full lift,
    # over formal roots represented by their exponent characters:
    # u=ab, v=ac, y=d, R=abcd, so R/(uy)=c and v/c=a.
    need((15^3^8)==4 and (5^4)==1 and (3^1)==2,'marked root reconstruction')
    elliptic_characters={7,11,13,14,15}
    need(elliptic_characters&C=={11,13,14},'three retained elliptic differential characters')
    for den in range(1,9):
        for num in range(-12,13):
            t=F(num,den);delta=F(5,6)
            need(delta*(1+t)*(1-t)*(2-t)==delta*(t+1)*(t-1)*(t-2),'first E1 twist')
            x=1-t
            need(delta*t*(1+t)*(2-t)==delta*(x+1)*(x-1)*(x-2),'second E1 twist')
            x=t-1
            need(delta*t*(1-t)*(2-t)==delta*(x**3-x),'E0 twist')
            if t:
                x=2/t
                need(4*delta*t*(1-t*t)*(2-t)/t**4==delta*(x-2)*(x+2)*(x-1),'E2 reciprocal twist')
    print(json.dumps({'status':'FINITE-EXACT PASS','rational_identity_cases':tests,
        'unit_character_spaces':[sorted(CU),sorted(CD)],'combined_character_space':sorted(C),
        'genus_by_involution_counts':genus_counts,'petersen_edges':petersen_edges,
        'free_V4_subgroups':len(free_subgroups),'star_hostile':star,'triangle_control':triangle,
        'scope':'exact chart identities and all marked deck involutions/pairs; rational-point exclusion uses the written infinite descent'},indent=2,sort_keys=True))
if __name__=='__main__':main()
