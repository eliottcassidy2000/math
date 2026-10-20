"""Small complete blue-presentation probe for local strong-edge detours."""
from collections import Counter
from importlib.util import module_from_spec,spec_from_file_location
from itertools import permutations
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
spec=spec_from_file_location('outer_base',HERE.parent/'spectrum_basins_20261020e/outer_probe.py')
O=module_from_spec(spec);spec.loader.exec_module(O);G=O.G


def paths(a):
    full=(1<<len(a))-1
    def walk(p,unused):
        if not unused:yield tuple(p);return
        allowed=unused if not p else a[p[-1]]&unused
        while allowed:
            bit=allowed&-allowed;allowed-=bit
            yield from walk(p+[bit.bit_length()-1],unused^bit)
    yield from walk([],full)


def replay(a,words):
    b=a
    for p in words:
        if not G.is_path(b,p) or G.blue(b,p):return False
        b=G.complement(b,p)
    return b==G.complement(a,tuple(range(len(a))))


def isomorphic(a,b):
    if sorted(x.bit_count() for x in a)!=sorted(x.bit_count() for x in b):return False
    choices=[[v for v in range(len(b)) if b[v].bit_count()==row.bit_count()] for row in a]
    order=sorted(range(len(a)),key=lambda u:len(choices[u]));mapping={};used=set()
    def walk(i):
        if i==len(a):return True
        u=order[i]
        for v in choices[u]:
            if v in used or any(bool(a[u]>>x&1)!=bool(b[v]>>y&1) for x,y in mapping.items()):continue
            mapping[u]=v;used.add(v)
            if walk(i+1):return True
            used.remove(v);del mapping[u]
        return False
    return walk(0)


def hub(a,black):
    n=len(a);p=tuple(range(n));d=[x.bit_count() for x in a]
    for q in black:
        if (q[0],q[-1])!=(0,n-1):continue
        for j,(u,v) in enumerate(zip(q,q[1:])):
            if v!=u+1:continue
            rq=q[j+1:]+q[:j+1];rp=p[v:]+p[:v]
            words=[q,rq,rp]
            if replay(a,words):
                return dict(kind='hub',paths=words,
                            score_certified=d[u]+d[v]-(u==0)+(v==n-1)!=n-1)


def square(a,black):
    n=len(a);p=tuple(range(n));edges=G.edge_set(p);d=[x.bit_count() for x in a]
    for q in black:
        if G.edge_set(q)&edges:continue
        words=[q,p[::-1],q[::-1]]
        if replay(a,words):
            u,v=q[0],q[-1];ds=d[u]+d[v]
            shifted=ds-(u==0)-(v==0)+(u==n-1)+(v==n-1)
            return dict(kind='square',paths=words,unpaired_ends=u+v!=n-1,
                        score_certified=ds!=n-1 and shifted!=n-1 and u+v!=n-1)


def suffix_templates(a):
    n=len(a);z=n-1
    for s in range(1,n-3):
        F=tuple(range(s));b=s+1;C=tuple(range(s+2,z))
        ear=[F+(b,)+C+(z,s),C+F+(b,z,s),C+F+(s,b,z)]
        if replay(a,ear):return dict(kind='suffix_ear',start=s,paths=ear)
        for i in range(s+1,n-2):
            L=tuple(range(s+1,i));b,c,d=i,i+1,i+2;R=tuple(range(i+3,n))
            words=[F+L+(b,c,d)+R+(s,), (c,)+F+L+(b,d)+R+(s,),
                   (c,)+F+(s,)+L+(b,d)+R]
            if replay(a,words):return dict(kind='suffix_triangle',start=s,triangle=i,paths=words)


def extraction(a):
    n=len(a)
    for s in range(1,n-3):
        A=tuple(range(s))
        for u in range(s+2,n-1):
            C=tuple(range(s+1,u))
            for k in range(u,n-1):
                K=tuple(range(u,k+1));D=tuple(range(k+1,n))
                # R paths converted to K by reversing the second word.
                words=[A+K+(s,)+C+D,C+D+A+K+(s,),C+K+D+A+(s,)]
                if replay(a,words):return dict(kind='extraction',start=s,u=u,k=k,paths=words)


def main(max_n=7):
    records=[];hostiles=[];witnesses={}
    for n in range(4,max_n+1):
        stat=Counter()
        for a in G.blue_tilings(n):
            if not(a[0]>>(n-1)&1):continue
            p=tuple(range(n));b=G.complement(a,p)
            if not O.strong(a) or not O.strong(b):continue
            hp=list(paths(a));black=[q for q in hp if not G.blue(a,q)]
            if not black or not any(not G.blue(b,q) for q in paths(b)):continue
            stat['strong_mixed_open_blue']+=1
            if isomorphic(a,b):stat['blue_loops']+=1;continue
            stat['nonloop']+=1
            certs=[f(a,black) for f in (hub,square)]
            certs += [suffix_templates(a),extraction(a)]
            for c in certs:
                if c:
                    stat[c['kind']]+=1
                    witnesses.setdefault(c['kind'],dict(rows=a,certificate=c))
            if not any(certs):
                stat['uncovered']+=1
                hostiles.append(dict(n=n,rows=a,scores=[x.bit_count() for x in a],
                                     hp=len(hp),black_hp=len(black)))
        records.append(dict(n=n,counts=dict(stat)))
        print(records[-1],flush=True)
    out=dict(universe='all standard blue presentations, n4..7; open side, strong mixed endpoints; ordinary blue loops excluded',
             records=records,uncovered=hostiles,witnesses=witnesses)
    (HERE/'detours_probe.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':main()
