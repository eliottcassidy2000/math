#!/usr/bin/env python3
"""Independent exact controls for the honeycomb-return bridge. No dependencies."""
from collections import deque
from itertools import product
from math import gcd, isqrt
import json

COLORS='abc'
def require(x,msg):
    if not x: raise RuntimeError(msg)
def step(s,c):
    u,v,e=s
    if c=='a': return u,v,1-e
    if c=='b': return u+(1 if e else -1),v,1-e
    return u,v+(1 if e else -1),1-e

def endpoint(word):
    s=(0,0,0)
    for c in word:s=step(s,c)
    return s

def ball(r):
    d={(0,0,0):0}; q=deque(d)
    while q:
        s=q.popleft()
        if d[s]==r:continue
        for c in COLORS:
            t=step(s,c)
            if t not in d:d[t]=d[s]+1;q.append(t)
    return d

def arith(s,a,b,c,x=0):
    u,v,e=s; t=x+(b-a)*u+(c-a)*v
    return (a-t if e else t),e

def cube(s):
    u,v,e=s
    return (u+v+e)%2,u%2,v%2

def ordinary(n): return n//2 if n%2==0 else 3*n+1

def square_graph(n,targets=None):
    return {x:[y for y in range(1,n+1) if x!=y and ((x+y in targets) if targets is not None else isqrt(x+y)**2==x+y)] for x in range(1,n+1)}

def main():
    out={}
    B=ball(32)
    for r in range(33):require(sum(d<=r for d in B.values())==1+3*r*(r+1)//2,('ball',r))
    for s,d in B.items():
        if s[2]==0:
            u,v,_=s; require(d==2*max(abs(u),abs(v),abs(u+v)),('norm',s))
        for ci,c in enumerate(COLORS):
            t=step(s,c)
            require(step(t,c)==s,('involution',s,c))
            cs=list(cube(s));cs[ci]^=1
            require(tuple(cs)==cube(t),('cube-cover',s,c))
            for a,b,z in [(4,9,16),(9,16,25),(16,25,36),(0,6,15)]:
                value,e=arith(s,a,b,z); value2,e2=arith(t,a,b,z)
                require(value+value2==[a,b,z][ci] and e+e2==1,('arith-cover',s,c))
    out['honeycomb_balls']={str(r):sum(d<=r for d in B.values()) for r in (15,16,20,21,24,26,32)}
    out['cube_quotient_vertices']=len({cube((u,v,e)) for u,v,e in product(range(2),range(2),range(2))})
    out['word_controls']={}
    for word in ['abcabc','abab','ab'*9+'cb'*7,('ab'*9+'cb'*7)*2]:
        s=endpoint(word)
        n=0
        for c in word:n={'a':9,'b':16,'c':25}[c]-n
        require((n,len(word)%2)==arith(s,9,16,25),('direct-affine',word))
        out['word_controls'][word]={'length':len(word),'honeycomb_endpoint':s,'arithmetic_endpoint':arith(s,9,16,25),'cube_endpoint':cube(s)}
    # Exhaust every word through length 10, including all backtracks.
    tested=0
    for n in range(11):
        for w in product(COLORS,repeat=n):
            s=endpoint(w);p=[0,0,0];v=0
            for c in w:p[COLORS.index(c)]^=1;v={'a':9,'b':16,'c':25}[c]-v
            require(tuple(p)==cube(s) and (v,n%2)==arith(s,9,16,25),('word',w))
            tested+=1
    out['exhaustive_words_through_length_10']=tested
    # Small direct quotient BFS proves no shorter extra return at the checked targets.
    out['deck_controls']=[]
    for a,b,c in [(4,9,16),(9,16,25),(16,25,36),(0,6,15)]:
        p,q=b-a,c-a;g=gcd(p,q);k=(q//g,-p//g,0);R=2*q//g
        D=ball(R)
        extra=[(d,s) for s,d in D.items() if s!=(0,0,0) and s[2]==0 and arith(s,a,b,c)==(0,0)]
        require(min(d for d,s in extra)==R,('systole',a,b,c))
        for radius in [R//2-1,R//2]:
            states=[s for s,d in D.items() if d<=radius]
            unique=len({arith(s,a,b,c) for s in states})==len(states)
            require(unique==(radius<R//2),('injectivity',a,b,c,radius))
        out['deck_controls'].append({'targets':[a,b,c],'gcd':g,'primitive_deck':k[:2],'first_extra_even_return':R,'last_injective_ball_radius':R//2-1,'first_joint_extra_return':2*R})
    # Recheck the Collatz obstruction by actual positive inverse enumeration.
    depths={5:0};front={5}
    out['collatz_balls']={}
    first_hex=first_tri=None
    for r in range(1,27):
        new={2*x for x in front}|{(x-1)//3 for x in front if x%6==4}
        require(not (new & depths.keys()),('tree-repetition',r))
        front=new
        for x in new:depths[x]=r
        if len(depths)>1+3*r*(r+1)//2 and first_hex is None:first_hex=r
        if len(depths)>1+3*r*(r+1) and first_tri is None:first_tri=r
        if r in [20,21,24,26]:out['collatz_balls'][str(r)]={'vertices':len(depths),'shell':len(front),'honeycomb_capacity':1+3*r*(r+1)//2,'triangular_capacity':1+3*r*(r+1)}
    for x,d in depths.items():
        y=x
        for _ in range(d):y=ordinary(y)
        require(y==5,('forward-collatz-check',x,d))
    require(first_hex==21 and first_tri==26,'first growth obstruction')
    out['collatz_first_capacity_violations']={'honeycomb':first_hex,'triangular':first_tri}
    # The three-target graph at 15 is literally the inherited unique path.
    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    P=square_graph(15,{9,16,25}); full=square_graph(15)
    path_edges={tuple(sorted(e)) for e in zip(path,path[1:])}
    edges=lambda A:{tuple(sorted((x,y))) for x,ys in A.items() for y in ys}
    require(edges(P)==path_edges,'three-target path')
    require(edges(full)-edges(P)=={(1,3)},'fourth-target chord')
    pentagon=[1,3,6,10,15,1]
    require([x+y for x,y in zip(pentagon,pentagon[1:])]==[4,9,16,25,16],'anchored pentagon')
    out['square_15']={'path':path,'three_target_edges':len(edges(P)),'fourth_target_chord':[1,3],'odd_cycle':pentagon}
    # Lift this path to H using its exact target colors, retaining the physical starting value.
    s=(0,0,0);lift=[s]
    for x,y in zip(path,path[1:]):
        s=step(s,{9:'a',16:'b',25:'c'}[x+y]);lift.append(s)
        require(arith(s,9,16,25,8)[0]==y,'clipped path lift')
    out['square_15']['honeycomb_lift']=lift
    # Constructive positive cycles, independently checking all edges and winding.
    out['positive_cycle_family']=[]
    for m in range(3,101):
        a,b,c=(2*m+1)**2,(2*m+3)**2,(2*m+5)**2
        p,q=b-a,c-a;g=gcd(p,q);r=q//g
        D=a+b-c+g;L=D//2
        cyc=[];word=[]
        for j in range(r):
            x=L+g*((j*(p//g))%r)
            color='a' if x<L+q-p else 'c'
            y=(a if color=='a' else c)-x
            cyc.extend((x,y));word.extend((color,'b'))
        require(len(cyc)==4*m+6 and len(set(cyc))==len(cyc),'simple mechanical cycle')
        require(min(cyc)>0 and max(cyc)==2*m*m+14*m+13,'family height')
        sums=[cyc[j]+cyc[(j+1)%len(cyc)] for j in range(len(cyc))]
        require(all(t in {a,b,c} for t in sums),'square labels')
        e=endpoint(word)
        require(e==(q//g,-p//g,0),'primitive winding')
        require(cube(e)!=(0,0,0) and cube(endpoint(word*2))==(0,0,0),'two-lap synchronization')
        if m<=8:
            out['positive_cycle_family'].append({'m':m,'targets':[a,b,c],'cycle_length':len(cyc),'height':max(cyc),'primitive_deck':e[:2],'tournament_holonomy':cube(e),'cycle':cyc,'word':''.join(word)})
    out['positive_cycle_family_checked_m']=[3,100]
    # A concrete merged-node separation: one 18-cycle changes transitive to cyclic.
    forward={(i,j) for i in range(4) for j in range(i+1,4)}
    after=set(forward)
    for edge in [(0,2),(0,3)]:after.remove(edge);after.add(edge[::-1])
    require(all(e in after for e in [(0,1),(1,2),(2,0)]),'actual tournament triangle')
    out['order_four_tournament_witness']={'forced_path':[0,1,2,3],'free_arc_color_order':[[0,2],[0,3],[1,3]],'initial_state':'all forward','after_one_lap_triangle':[0,1,2,0],'after_two_laps':'all forward'}
    out['all_checks_passed']=True
    print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
