"""All-order capped closure and exact-safe rooted charge certificate."""
import importlib.util,json,hashlib,random
from collections import deque
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def load(name,relative,digest):
    p=ROOT.parent/relative
    if hashlib.sha256(p.read_bytes()).hexdigest()!=digest:raise RuntimeError(('dependency hash',str(p)))
    spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
f=load('frontier','ribbon/transfer.py','2ef89174390c503494d4e88e81f1b7ee149f985b6746f2912430ae88e46a289d')
t=load('graph_counter','ribbon-terminals/verify.py','9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e')
require=t.require
W=3;M=13;C=27
states,mat=f.automaton(W)
HH=[i for i,s in enumerate(states) if len(s)==1]
EE=[i for i,s in enumerate(states) if len(s)==1 and s[0][1]==2]
RR=[i for i,s in enumerate(states) if 1<=len(s)<=2 and any(a==2 for a,b in s)]
initial=(1,)+(0,)*(len(states)-1)

def read(v):return tuple(sum(v[i] for i in ids) for ids in (HH,EE,RR))

def build():
    seen={initial:''};q=deque([initial]);witness={}
    while q:
        v=q.popleft();word=seen[v];h,e,R=read(v)
        if h and e+R==2*h-2 and str(h) not in witness:witness[str(h)]=dict(word=word,profile=[h,e,R])
        for a in range(8):
            nv=f.advance(v,mat,a,C)
            if sum(nv[i] for i in HH)>M:continue
            if nv not in seen:seen[nv]=word+str(a);q.append(nv)
    return dict(width=W,bound=M,cap=C,states=[list(v) for v in seen],witnesses=witness)

def main():
    cert=json.loads((ROOT/'closure.json').read_text());fresh=build();require(cert==fresh,'closure differs from regeneration')
    bank={tuple(v) for v in cert['states']};require(initial in bank,'initial missing');edges=0;pruned=0
    for v in bank:
        require(len(v)==len(states) and all(type(x) is int and 0<=x<=C for x in v),'invalid vector')
        h,e,R=read(v);require(h<=M,'invalid retained count')
        require(not(h==13 and e+R==2*h-2),'forbidden count present')
        for a in range(8):
            nv=f.advance(v,mat,a,C)
            if sum(nv[i] for i in HH)>M:pruned+=1
            else:require(nv in bank,'closure failed');edges+=1
    physical=0
    def walk(word,v):
        nonlocal physical
        n=len(word)
        if n:
            rows=t.wordgraph(word,W);actual=t.profile(rows,n-1)
            require(actual==read(v),('root readout',word,actual,read(v)));physical+=1
        if n==6:return
        for x in range(1<<min(n,W)):
            a=x<<max(0,W-n);walk(word+str(a),f.advance(v,mat,a))
    walk('',initial)
    for h,z in cert['witnesses'].items():
        rows=t.wordgraph(z['word'],W);require(list(t.profile(rows,len(rows)-1))==z['profile'],('witness',h))
    require(sorted(map(int,cert['witnesses']))==[3,5,9,11],'unexpected spectrum')
    # Algebraic stress checks include very large coordinates and explicitly
    # verify that saturation cannot fabricate a successful charge identity.
    rng=random.Random(197);satchecks=0
    for _ in range(1000):
        v=tuple(rng.randrange(1000) if rng.randrange(4)==0 else 0 for _ in states);a=rng.randrange(8)
        exact=f.advance(v,mat,a);clipped=f.advance(tuple(min(C,x) for x in v),mat,a,C)
        require(tuple(min(C,x) for x in exact)==clipped,'saturation identity')
        rawR=sum(v[i] for i in RR);clipR=sum(min(C,v[i]) for i in RR)
        require(clipR>=min(rawR,C),'readout headroom')
        satchecks+=1
    alt=t.wordgraph('00113',2);altp=t.profile(alt,4);require(altp==(13,2,18),'ordinary count13 witness')
    small13=sorted({(e,R) for v in bank for h,e,R in [read(v)] if h==13 and R<C})
    print(json.dumps(dict(width=W,bound=M,cap=C,frontier_coordinates=len(states),closed_vectors=len(bank),retained_edges=edges,pruned_edges=pruned,charge_counts=sorted(map(int,cert['witnesses'])),physical_readout_checks=physical,saturation_stress_checks=satchecks,low_R_profiles_at_H13=small13,ordinary_H13_witness=dict(word='00113',width=2,profile=altp),certificate_sha256=hashlib.sha256((ROOT/'closure.json').read_bytes()).hexdigest()),sort_keys=True,indent=2))

if __name__=='__main__':main()
