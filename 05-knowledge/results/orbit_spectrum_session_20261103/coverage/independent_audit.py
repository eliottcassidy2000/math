#!/usr/bin/env python3
"""Independent small literal controls for moving-mark plus-two repair.

No producer mathematical code is imported. Actual adjacency is constructed
from named arcs; Hamilton and ending-bin counts use literal permutations.
"""
from itertools import permutations
from math import gcd
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise RuntimeError(message)

def source(m):
    # Physical names x=0, y=1, c_i=i+1.
    n=m+2;a=[[False]*n for _ in range(n)]
    a[0][1]=True
    for c in range(2,n):a[1][c]=True;a[c][0]=True
    for c in range(2,n):
        for d in range(c+1,n):a[c][d]=True
    return a

def repaired(m):
    need(m>=2,'generic repair starts at m2')
    old=source(m);n=len(old);a=[row+[False,False] for row in old]+[[False]*(n+2) for _ in range(2)]
    ca=n-1;newa=n;newb=n+1
    for v in range(n):a[v][newa]=True
    for v in range(n):
        if v==ca:a[newb][v]=True
        else:a[v][newb]=True
    a[newa][newb]=True
    return a

def root_birth(a,root):
    n=len(a);b=[r+[False] for r in a]+[[False]*(n+1)]
    for v in range(n):
        if v==root:b[v][n]=True
        else:b[n][v]=True
    return b

def subset(a,vertices):
    return [[a[u][v] for v in vertices] for u in vertices]

def literal_profile(a,root):
    n=len(a);h=s=F=0
    for p in permutations(range(n)):
        bad=[i for i in range(n-1) if not a[p[i]][p[i+1]]]
        if not bad:h+=1;s+=p[0]==root
        # The unique cut after root removes exactly the next edge constraint.
        if not bad or (len(bad)==1 and bad[0]==p.index(root)):F+=1
    rest=[v for v in range(n) if v!=root]
    q=sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in permutations(rest))
    return [h,s,q,F]

def strong(a):
    n=len(a)
    for start in range(n):
        todo=[start];seen={start}
        while todo:
            v=todo.pop()
            for w in range(n):
                if a[v][w] and w not in seen:seen.add(w);todo.append(w)
        if len(seen)!=n:return False
    return True

def two_apex_word(word):
    n=len(word)+2;a=[[False]*n for _ in range(n)];a[0][1]=True
    for i,b in enumerate(word,2):
        for apex in (0,1):
            if b>>apex&1:a[apex][i]=True
            else:a[i][apex]=True
        for j in range(2,i):a[j][i]=True
    return a

def main():
    E=two_apex_word((2,3,2));exterior_profile=literal_profile(E,2)
    need(exterior_profile==[9,3,5,11],'five-vertex exterior profile')
    bcontrols=[]
    for m in (2,3,4):
        b=repaired(m);root=0;profilex=literal_profile(b,root)
        need(strong(b),'repaired strong source')
        need(profilex[0]==3*2**m+5,'all-size repair count at control')
        need(profilex[3]==6*2**m-3,'ending-two-bin count at control')
        bcontrols.append({'m':m,'profile_at_old_x':profilex,
                          'rows':[sum((1<<j) for j,v in enumerate(row) if v) for row in b]})
    balanced=[]
    for m in (2,3):
        b=repaired(m);first=len(b)
        a=root_birth(b,0);need(strong(a),'first reset intermediate strong')
        a=root_birth(a,first);need(strong(a),'completed reset strong')
        profile=literal_profile(a,first)
        t=3*2**m+5;h=15*2**m-1
        need(profile==[h,t,t,h],'literal balanced root profile')
        balanced.append({'m':m,'order':len(a),'profile':profile,'root':first,
                         'rows':[sum((1<<j) for j,v in enumerate(row) if v) for row in a]})
    # The first source must use the separately stated base suffix.
    normal1=literal_profile(E,0)[0]
    exceptional=two_apex_word((2,2,0))
    exception=literal_profile(exceptional,0)[0]
    need(normal1==9 and exception==11 and strong(exceptional),'m1 exceptional branch')
    residues=[]
    for m in range(2,98):
        h=15*2**m-1;t=3*2**m+5
        need(gcd(h,t)==(13 if m%12==11 else 1),'exact primitive guard')
        need((h%7==0)==(m%3==0),'7 divisibility guard')
        if m%3==0:
            j=m//3
            need((h%49==0)==(j%7==5),'49 divisibility guard')
            residues.append({'j':j,'mod49':h%49})
    result={'status':'PASS','scope':'Independent targeted literal and arithmetic audit; no all-target coverage claim',
        'exterior_profile_at_y':exterior_profile,
        'repaired_sources':bcontrols,'balanced_children':balanced,
        'base_m1':{'generic_232_H':normal1,'exception_220_H':exception},
        'arithmetic_controls_m':[2,97],'subfamily_49_residues':residues,
        'note_sha256':hashlib.sha256((HERE/'NOTE.md').read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
