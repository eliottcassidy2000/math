#!/usr/bin/env python3
"""Additional portable integer checks, independent of check.py's functions."""
from collections import deque
from math import gcd
import json

def require(ok, label):
    if not ok:
        raise RuntimeError(str(label))

def honey_ball(radius):
    distances={(0,0,0):0}
    queue=deque(distances)
    while queue:
        u,v,e=queue.popleft()
        d=distances[u,v,e]
        if d==radius:
            continue
        sign=1 if e else -1
        for nxt in [(u,v,1-e),(u+sign,v,1-e),(u,v+sign,1-e)]:
            if nxt not in distances:
                distances[nxt]=d+1
                queue.append(nxt)
    return distances

def image(vertex,a,b,c):
    u,v,e=vertex
    value=(b-a)*u+(c-a)*v
    return (a-value if e else value),e

def bits(vertex):
    u,v,e=vertex
    return (u+v+e)%2,u%2,v%2

def main():
    triple_count=cycle_count=simple_count=0
    for a in range(-5,26):
        for b in range(a+1,31):
            for c in range(b+1,36):
                p,q=b-a,c-a
                g=gcd(p,q)
                r,s=q//g,p//g
                D=a+b-c+g
                if D<2:
                    continue
                triple_count+=1
                heights=[]
                for L in range(1,D):
                    xs=[L+g*((j*s)%r) for j in range(r)]
                    ys=[(a if x<L+q-p else c)-x for x in xs]
                    values=[z for pair in zip(xs,ys) for z in pair]
                    require(min(values)>0 and len(set(xs))==r and len(set(ys))==r,('vertices',a,b,c,L))
                    require(all(xs[j]+ys[j] in (a,c) and ys[j]+xs[(j+1)%r]==b for j in range(r)),('edges',a,b,c,L))
                    if (2*L-a)%g:
                        simple_count+=1
                        require(len(set(values))==2*r,('simple',a,b,c,L))
                    height=max(values)
                    heights.append(height)
                    require(height==max(L+q-g,b-L),('height',a,b,c,L))
                    cycle_count+=1
                require(min(heights)==(b+c-a-g+1)//2,('optimum',a,b,c))
    controls=[]
    for a,b,c in [(4,9,16),(9,16,25),(16,25,36),(0,6,15),(1,2,3),(3,6,9),(10,12,16)]:
        g=gcd(b-a,c-a)
        r=(c-a)//g
        B=honey_ball(4*r)
        joint=[d for vertex,d in B.items() if vertex!=(0,0,0) and image(vertex,a,b,c)==(0,0) and bits(vertex)==(0,0,0)]
        require(min(joint)==4*r,('joint',a,b,c))
        centers=[(0,0,0),(0,0,1),(5,-3,0),(5,-3,1)]
        for U,V,E in centers:
            for R in [r-1,r]:
                vertices=[(U+(-u if E else u),V+(-v if E else v),1-e if E else e) for (u,v,e),d in B.items() if d<=R]
                images=[image(vertex,a,b,c) for vertex in vertices]
                injective=len(images)==len(set(images))
                require(injective==(R<r),('radius',a,b,c,(U,V,E),R))
        controls.append({'targets':[a,b,c],'first_joint_extra_return':min(joint),'first_noninjective_ball_radius':r,'centers_checked':centers})
    print(json.dumps({'all_checks_passed':True,'general_mechanical_triples':triple_count,'general_layered_cycles':cycle_count,'simple_cycles':simple_count,'joint_and_radius_controls':controls},indent=2,sort_keys=True))

if __name__=='__main__':
    main()
