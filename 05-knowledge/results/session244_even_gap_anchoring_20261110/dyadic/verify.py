#!/usr/bin/env python3
"""Literal permutation verification; no producer or sibling imports."""
import itertools,json,sys

def need(b,msg):
    if not b: raise RuntimeError(msg)

def profiles(a):
    n=len(a);full=(1<<n)-1;H=[0]*(full+1);H[0]=1
    starts=[[0]*n for _ in H];ends=[[0]*n for _ in H]
    for s in range(1,full+1):
        v=[x for x in range(n) if s>>x&1]
        for p in itertools.permutations(v):
            if all(a[x]>>y&1 for x,y in zip(p,p[1:])):
                H[s]+=1;starts[s][p[0]]+=1;ends[s][p[-1]]+=1
    return [(H[full],ends[full][r],sum(starts[s][r]*H[full^s] for s in range(1,full+1)),H[full^(1<<r)]) for r in range(n)]

def signs(states):return sorted({U-D-d for D,E,U,d in states if E==d and abs(U-D-d)==1})

def extend(a,r):return [x|(1<<len(a)) if i!=r else x for i,x in enumerate(a)]+[1<<r]

if '--fail' in sys.argv:need(False,'intentional failure')
a=[8,1,3,6];graphs=[a,extend(a,0),extend(a,3)];records=[]
for g,expected in zip(graphs,([-1,1],[-1],[1])):
    s=profiles(g);need(signs(s)==expected,'complete root menu')
    records.append({'rows':g,'profiles':s,'eligible_signs':signs(s)})
for D in range(1,30,2):
    for d in range(1,30,2):
        for c in (-1,1):
            E=d;U=D+d+c
            z=(D+U-E,d,2*U-d,d)
            need(z==(2*D+c,d,2*D+d+2*c,d),'formal transfer')
print(json.dumps({'status':'PASS','method':'independent literal permutations on all subsets','literal_graphs':records,'formal_ray_checks':450},indent=2))
