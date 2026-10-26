#!/usr/bin/env python3
"""Bounded, independent tc scope checks. Standard library; no producer imports.
Graph automorphisms use generic false-twin compression and edge-colour backtracking.
Incoming n=10 source files are hashed, never imported or run.
"""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter
from functools import lru_cache
from math import factorial, prod
from fractions import Fraction
import hashlib, json

HERE=Path(__file__).resolve().parent

def need(test, why):
    if not test: raise RuntimeError(why)

def tournament(n,b):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(b>>k)&1; a[j][i]=1-a[i][j]
    return a

def words(a):
    return [p for p in permutations(range(len(a))) if all(a[x][y] for x,y in zip(p,p[1:]))]

def aut(a):
    n=len(a)
    return [p for p in permutations(range(n)) if all(a[i][j]==a[p[i]][p[j]] for i in range(n) for j in range(n))]

def key(a):
    return min(tuple(a[p[i]][p[j]] for i,j in combinations(range(len(a)),2)) for p in permutations(range(len(a))))

def codec(a,joint):
    n=len(a); w=7 if joint else 9; g=[[0]*(w*n) for _ in range(w*n)]
    def edge(x,y,c):g[x][y]=g[y][x]=c
    for v in range(n):
        if joint:
            for x,y in ((0,1),(0,2),(2,3),(0,4),(4,5),(5,6)):edge(w*v+x,w*v+y,1)
        else:
            for x,y in ((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(2,7),(2,8)):edge(w*v+x,w*v+y,1)
    for u in range(n):
        for v in range(n):
            if a[u][v]:edge(w*u+1,w*v+(3 if joint else 2),2 if joint else 1)
    return g

def decode(g,joint):
    N=len(g); tails=[];heads=[]
    if joint:
        e=[{j for j in range(N) if g[i][j]&1} for i in range(N)]
        # Only centers have frame degree 3; distinguish arms by their lengths.
        for c in range(N):
            if len(e[c])!=3:continue
            arms={}
            for v in e[c]:
                arm=[];prev=c
                while True:
                    arm.append(v);nxt=e[v]-{prev}
                    if not nxt:break
                    need(len(nxt)==1,'bad frame arm');prev,v=v,next(iter(nxt))
                arms[len(arm)]=arm
            need(set(arms)=={1,2,3},'asymmetric frame')
            tails.append(arms[1][-1]);heads.append(arms[2][-1])
    else:
        e=[{j for j in range(N) if g[i][j]} for i in range(N)]
        cores=[i for i in range(N) if len(e[i])>1]
        leafcount={i:sum(len(e[j])==1 for j in e[i]) for i in cores}
        for c in cores:
            if leafcount[c]!=1:continue
            ts=[i for i in e[c] if leafcount.get(i)==2];hs=[i for i in e[c] if leafcount.get(i)==3]
            need(len(ts)==len(hs)==1,'port types')
            tails+=ts;heads+=hs
    return [[int(bool(g[t][h]&(2 if joint else 1))) for h in heads] for t in tails]

def graph_aut(g):
    # Equal full rows are false twins. Every permutation in such a class is
    # an automorphism; all automorphisms induce a weight-preserving quotient map.
    classes={}
    for i,row in enumerate(g):classes.setdefault(tuple(row),[]).append(i)
    banks=list(classes.values());w=[len(c) for c in banks];r=len(banks)
    q=[[g[a[0]][b[0]] for b in banks] for a in banks]
    colors=[(w[i],tuple(sorted(Counter((q[i][j],w[j]) for j in range(r) if q[i][j]).items()))) for i in range(r)]
    def compress(xs):
        names={x:k for k,x in enumerate(sorted(set(xs)))};return [names[x] for x in xs]
    colors=compress(colors)
    while True:
        new=compress([(colors[i],tuple(sorted(Counter((q[i][j],colors[j]) for j in range(r) if q[i][j]).items()))) for i in range(r)])
        if len(set(new))==len(set(colors)):colors=new;break
        colors=new
    domains=[{j for j in range(r) if colors[j]==colors[i]} for i in range(r)]
    def search(mapping,unused):
        if not unused:return 1
        options={i:[j for j in domains[i] if j not in mapping.values() and all(q[i][k]==q[j][l] for k,l in mapping.items())] for i in unused}
        i=min(unused,key=lambda x:len(options[x]));ans=0
        for j in options[i]:ans+=search(mapping|{i:j},unused-{i})
        return ans
    return search({},set(range(r)))*prod(factorial(x) for x in w)

def omega(a):
    n=len(a);cs=[]
    for m in range(3,n+1,2):
        for support in combinations(range(n),m):
            for tail in permutations(support[1:]):
                p=(support[0],)+tail
                if all(a[p[i]][p[(i+1)%m]] for i in range(m)):cs.append(set(p))
    edges=[sum(1<<j for j,d in enumerate(cs) if i!=j and c&d) for i,c in enumerate(cs)]
    @lru_cache(None)
    def ind(s):
        if not s:return 1
        bit=s&-s;i=bit.bit_length()-1
        return ind(s^bit)+2*ind((s^bit)&~edges[i])
    return ind((1<<len(cs))-1)

def main():
    pins=json.loads((HERE/'SOURCES.json').read_text());verified=[]
    for f in pins['files']:
        data=(HERE/f['local_path']).read_bytes()
        need(hashlib.sha256(data).hexdigest()==f['sha256'],'source sha256')
        if 'git_blob' in f:need(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==f['git_blob'],'source Git blob')
        verified.append(f['local_path'])
    reps={};summary=[];orbitchecks=0
    for n in range(6):
        spectra=set();count=0;nontrivial=0;reps[n]={}
        for b in range(1<<(n*(n-1)//2)):
            a=tournament(n,b);hp=words(a);gs=aut(a);H=len(hp);A=len(gs)
            need(H%A==0,'path count divisibility');need(H>0,'positive path count')
            remaining=set(hp);orbits=0
            while remaining:
                p=min(remaining);orb={tuple(g[v] for v in p) for g in gs}
                need(len(orb)==A,'free path action');need(orb<=remaining,'orbits partition');remaining-=orb;orbits+=1;orbitchecks+=1
            need(orbits==H//A,'orbit quotient');need(omega(a)==H,'OCF small control')
            for joint in (False,True):need(decode(codec(a,joint),joint)==a,'exact unpermuted codec decode')
            need(H//A not in (7,21),'bounded tc gaps only')
            spectra.add(H//A);count+=1;nontrivial+=A>1;reps[n].setdefault(key(a),a)
        summary.append({'n':n,'labelled_tournaments':count,'isomorphism_classes':len(reps[n]),'nontrivial_automorphism_cases':nontrivial,'tc_spectrum':sorted(spectra)})
    need([len(reps[n]) for n in range(6)]==[1,1,1,2,4,12],'known small class counts')
    codec_tests=[]
    for n in range(6):
        for k,a in reps[n].items():
            A=len(aut(a));H=len(words(a));record={'n':n,'key':list(k),'H':H,'AutT':A,'tc':H//A}
            for joint in (True,False):
                g=codec(a,joint);observed=graph_aut(g);factor=1 if joint else 12**n
                need(observed==factor*A,'generic graph automorphism count')
                need(Fraction(H*factor,observed)==H//A,'correct graph quotient')
                # Reverse all graph labels as a decoder control, without role hints.
                scrambled=[list(reversed(row)) for row in reversed(g)]
                need(key(decode(scrambled,joint))==k,'unlabelled decoder class')
                record['AutJoint' if joint else 'AutGamma']=observed
            codec_tests.append(record)
    c3=[[0,1,0],[0,0,1],[1,0,0]];h=len(words(c3));a=len(aut(c3));ga=graph_aut(codec(c3,False))
    need((h,a,ga)==(3,3,5184),'C3 hostile')
    incoming=json.loads((HERE/'sources/tc_fibres_n10.json').read_text())
    # This checks only what the stored output says, NOT correctness/completeness of its census.
    incoming_digest={name:{'tournaments':row['tournaments'],'strong':row['strong'],'candidates_checked':row['candidates_checked'],'reported_hits':len(row['tc_7_or_21_found'])} for name,row in incoming.items()}
    out={'status':'PASS','independent_scope':'all 1100 labelled tournaments on orders 0..5, and generic graph automorphism counts on all 21 isomorphism classes at those orders; no n=10 execution','sources_verified':verified,'small_orders':summary,'path_orbits_explicitly_checked':orbitchecks,'graph_codec_checks':codec_tests,'wrong_denominator_hostile':{'tournament':'C3','H':h,'AutT':a,'AutGamma':ga,'tc':h//a,'naive_H_over_AutGamma':str(Fraction(h,ga))},'incoming_n10_status':'REPORTED ONLY; source bytes pinned; no replay or certification of reduction/coverage','incoming_n10_stored_output':incoming_digest}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':main()
