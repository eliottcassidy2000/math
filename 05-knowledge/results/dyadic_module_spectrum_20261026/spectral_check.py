#!/usr/bin/env python3
"""Integer structural certificates underlying the all-size nuclear-norm proof."""
from random import Random
import json

def check(p,m):
    if not p:raise RuntimeError(m)

def run():
    rng=Random(20261026);cases=[]
    for n in (3,4,7,11,19,31,63,127):
        for arity in (3,5,7):
            for repeat in range(4):
                vertices=list(range(n));rng.shuffle(vertices)
                arcs=set();rectangles=[];light_depth=[0]*n;bill=0;light_mass=0
                def tree(vs):
                    nonlocal bill,light_mass
                    m=len(vs)
                    if m==1:return vs
                    q=rng.randint(2,min(arity,m))
                    cuts=[0]+sorted(rng.sample(range(1,m),q-1))+[m]
                    children=[vs[a:b] for a,b in zip(cuts,cuts[1:])]
                    # An arbitrary quotient, not only a cyclic quotient.
                    win={(i,j) if rng.randrange(2) else (j,i) for i in range(q) for j in range(i+1,q)}
                    h=max(range(q),key=lambda j:len(children[j]))
                    order=[i for i in range(q) if (i,h) in win]+[h]+[i for i in range(q) if (h,i) in win]
                    pos={x:i for i,x in enumerate(order)}
                    L=sum(len(c) for i,c in enumerate(children) if i!=h)
                    light_mass+=L
                    for i,c in enumerate(children):
                        if i!=h:
                            check(2*len(c)<=m,'every light child shrinks by half')
                            for x in c:light_depth[x]+=1
                    local_bill=0
                    for i,j in win:
                        arcs.update((x,y) for x in children[i] for y in children[j])
                        if pos[i]>pos[j]:
                            check(i!=h and j!=h,'heavy incident arcs agree with observer order')
                            rectangles.append((children[j][:],children[i][:]))
                            local_bill+=2*(len(children[i])+len(children[j]))
                    check(local_bill<=2*(q-2)*L,'rectangle trace-norm invoice')
                    bill+=local_bill
                    orders={i:tree(children[i]) for i in range(q)}
                    return sum([orders[i] for i in order],[])
                order=tree(vertices);pos={v:i for i,v in enumerate(order)}
                backwards={(x,y) if pos[x]<pos[y] else (y,x) for x,y in arcs if pos[x]>pos[y]}
                rect_pairs=set();count=0
                for A,B in rectangles:
                    for a in A:
                        for b in B:
                            check(pos[a]<pos[b],'rectangle orientation')
                            rect_pairs.add((a,b));count+=1
                check(rect_pairs==backwards and len(rect_pairs)==count,'LCA rectangles own each backward pair exactly once')
                check(sum(light_depth)==light_mass,'leaf charges retain multiplicity')
                check(max(light_depth)<=n.bit_length()-1,'light depth bound')
                check(bill<=2*(arity-2)*n*(n.bit_length()-1),'all-size light invoice')
                # Balanced binary decomposition of the transitive staircase.
                levels=[]
                def balanced(m,d=0):
                    if m<=1:return
                    if len(levels)<=d:levels.append(0)
                    levels[d]+=m
                    balanced(m//2,d+1);balanced(m-m//2,d+1)
                balanced(n)
                check(len(levels)<=(n-1).bit_length() and max(levels)<=n,'staircase nuclear invoice')
                cases.append(dict(n=n,arity=arity,repeat=repeat,backwards=len(backwards),light_mass=light_mass,
                                  rectangle_invoice=bill,staircase_invoice=sum(levels)))
    return dict(status='PASS',cases=cases,scope='96 finite structural controls; all-size norm inequalities proved in README')

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
