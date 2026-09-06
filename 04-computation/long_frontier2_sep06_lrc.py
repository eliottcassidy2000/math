#!/usr/bin/env python3
"""Joint interval positions close a declared 84-scale grid class and 81 actual entries."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import gcd,prod
from pathlib import Path
import hashlib
import json
import subprocess

Q=91**6
DENOMS=(128,101,113,131,149,167,227)
H=prod(DENOMS)
U=tuple(sorted(H//q for q in DENOMS))
LEAVES=(287,277,271,265,263)
P=prod(LEAVES)
V=tuple(sorted((P,)+tuple(69*(P//q) for q in LEAVES)))
PAIRS=((101,128),(113,128))
GATES=0


def check(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def atlas(a,b):
    c=(a+b)//gcd(a,b)
    if c>356:return False
    p=2
    while p*p<=c:
        e=0
        while c%p==0:c//=p;e+=1
        if e and (p%3!=2 or e>2):return False
        p+=1
    return c==1 or c%3==2


def graph(row):
    unseen=set(range(len(row)));groups=[]
    while unseen:
        first=min(unseen);unseen.remove(first);seen={first};pending=[first]
        while pending:
            i=pending.pop()
            for j in list(unseen):
                if atlas(row[i],row[j]):unseen.remove(j);seen.add(j);pending.append(j)
        groups.append(sorted(seen))
    return groups


def intervals(p,q):
    """Literal open arc intersections, coordinates in units1/(14pq)."""
    L=14*p*q
    aa=[(0,q)]+[(14*k*q-q,14*k*q+q)for k in range(1,p)]+[(L-q,L)]
    bb=[(0,p)]+[(14*k*p-p,14*k*p+p)for k in range(1,q)]+[(L-p,L)]
    i=j=0;out=[]
    while i<len(aa) and j<len(bb):
        lo,hi=max(aa[i][0],bb[j][0]),min(aa[i][1],bb[j][1])
        if lo<hi:out.append((lo,hi))
        if aa[i][1]<bb[j][1]:i+=1
        elif bb[j][1]<aa[i][1]:j+=1
        else:i+=1;j+=1
    check(out[0][0]==0 and out[-1][1]==L,'artificial origin split present')
    out=[(out[-1][0],out[0][1]+L)]+out[1:-1]
    cap=(p+q+13)//14
    expected=sorted([2*p]+[min(2*p,p+q-14*k)for k in range(1,cap)for _ in range(2)])
    check(sorted(b-a for a,b in out)==expected,'independent clipped-length formula')
    check(len(out)==2*cap-1,'exact circle component count')
    return L,out


def joint_minimum(p,q,t):
    L,arcs=intervals(p,q)
    events=defaultdict(lambda:[0,0]);before=0
    for lo,hi in arcs:
        check(0<t*(hi-lo)<L,'each projected interval shorter than one circle')
        a,b=(t*lo)%L,(t*hi)%L
        check(a!=b,'nonzero projected arc')
        events[a][0]+=1;events[b][1]+=1;before+=a>b
    value=before;best=len(arcs)+1;arg=None;trace=[]
    for a,(starts,ends) in sorted(events.items()):
        at=value-ends
        check(at>=0,'nonnegative open-endpoint coverage')
        trace.append((a,at,at+starts))
        if at<best:best,arg=at,a
        value=at+starts
    check(value==before,'cyclic coverage closes')
    return best,F(arg,L),L,trace


def literal_pair_count(p,q,t,alpha_num,alpha_den):
    den=alpha_den*t;count=0
    for j in range(t):
        x=alpha_num+alpha_den*j
        rp,rq=p*x%den,q*x%den
        if 14*min(rp,den-rp)<den and 14*min(rq,den-rq)<den:count+=1
    return count


def audit_all_cells(p,q,t):
    best,_,L,trace=joint_minimum(p,q,t)
    vals=[]
    for i,(a,at,after) in enumerate(trace):
        b=trace[(i+1)%len(trace)][0]+(L if i+1==len(trace) else 0)
        at_direct=literal_pair_count(p,q,t,a,L)
        middle_direct=literal_pair_count(p,q,t,a+b,2*L)
        check(at_direct==at,'literal grid at every critical phase')
        check(middle_direct==after,'literal grid in every open phase cell')
        vals.extend((at_direct,middle_direct))
    check(min(vals)==best,'independent full-cell minimum')


def ceildiv(a,b):return -((-a)//b)


def support(a,b,y):
    a,b=sorted((a,b));D=gcd(a,b);aa,bb=a//D,b//D
    delta=gcd(D,y);c,x=D//delta,y//delta
    r=pow(aa,-1,bb)*x%bb;s=(x-aa*r)//bb
    low=max(ceildiv(-Q-r,bb),ceildiv(s-Q,aa))
    high=min((Q-r)//bb,(s+Q)//aa)
    return c,x,aa+bb,c<=Q and low<=high


def profiles():
    root=Path(__file__).resolve().parent.parent
    rel='05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json'
    path=root/rel
    raw=path.read_bytes()if path.exists()else subprocess.check_output(['git','show','HEAD:'+rel],cwd=root)
    check(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','full profile bank pin')
    return {int(d):{(c,tuple(gs))for c,gs in row['profiles']}for d,row in json.loads(raw)['levels'].items()}


def main():
    generic=[t for t in range(620,791)if gcd(t,128*101*113)==1]
    check(len(generic)==84,'complete generic scale universe')
    counts=[];failures=[]
    for t in generic:
        excess=7*ceildiv(t,7)-t
        pair_data=[joint_minimum(p,q,t)for p,q in PAIRS]
        c1,c2=[r[0]for r in pair_data]
        check(max(c1,c2)>=excess+5,'uniform located overlap certificate')
        if c1<=excess:failures.append(t)
        counts.append([t,excess,c1,c2])
    check(failures==[633,667,687,741,761],'complete primary-pair failure set')
    check(min(max(c1,c2)-e for t,e,c1,c2 in counts)==5,'reported minimum surplus')
    for t in (621,633,645,687,761):
        for p,q in PAIRS:audit_all_cells(p,q,t)

    check(V==(360993824985,374026093035,382307113545,390963123663,393936227265,1501525040155),'actual six-shape')
    check(U==(4761938937472,6472815202432,7254766032256,8251604113024,8445001084423,9566018927488,10702575631744),'actual seven-shape')
    check(gcd(*V)==gcd(*U)==1 and all(v%2 for v in V),'primitive shapes and odd smaller shape')
    check(1 not in V+U and all(gcd(v,u)==1 for v in V for u in U),'unitless separated prime support')
    check(graph(V)==[list(range(6))] and graph(U)==[list(range(7))],'actual connected component graphs')
    dmin=min(gcd(v,w)for v,w in combinations(V,2))
    smax=max((u+w)//gcd(u,w)for u,w in combinations(U,2))
    check(dmin==1303226805 and smax==394,'complete entry extrema')
    check(max(Q//dmin+1,Q*smax//min(V)+1)==620,'exact sufficient entry threshold')
    check(790*sum(V)+sum(U)<Q*Q,'whole actual family physical box')
    check(620*min(V)+min(U)>356,'cross-coprime pairs cannot be atlas edges')
    for u,w in combinations(U,2):
        D=gcd(u,w);p,q=u//D,w//D
        check(790<7*q,'all individual interval credits vanish for whole band')
        R=Q*(p+q)-(p-1)*(q-1)
        for v in V:
            check(F(6*R,56*max(U)*v)<1,'all126 forward native widths fail')
            if w==max(U):check(56*D*v>6*Q,'all36 endpoint widths fail')
    for a,b in combinations(V,2):
        for u in U:check(F(6*Q,56*max(U)*gcd(a,b))<1,'all105 dual native widths fail')
    check(6*790<56*max(U)and 7<49*max(V),'both whole-arc sufficient gates fail')
    check(790<7*min(U),'all nested-origin credits vanish')

    bank=profiles()
    for d in range(1,7):check((1,(1,)*d)in bank[d],'complete inherited all-one complement word')
    actual=[t for t in range(620,791)if gcd(t,H)==1]
    check(len(actual)==81 and set(actual)<=set(generic),'complete actual scale universe')
    physical=[]
    for t in actual:
        row=tuple(t*v for v in V)+U
        check(gcd(*row)==1 and len(set(row))==13 and sum(row)<Q*Q,'actual primitive distinct boxed row')
        check(graph(row)==[list(range(6)),list(range(6,13))],'literal actual six-plus-seven graph')
        safe=[]
        for j in range(t):
            nums=[u*(2*j+1)%(2*t)for u in U]
            if all(7*min(r,2*t-r)>=t for r in nums):safe.append(j)
        check(len(safe)>=5,'literal actual half-grid safe count')
        j=safe[0];phase=F(2*j+1,2*t)
        clearance=min(min((v*phase)%1,1-(v*phase)%1)for v in row)
        check(clearance>F(1,14),'literal strict whole physical phase')
        check(sum((u%2)==1 for u in U)==1,'only one possible endpoint tail on odd half-grid')
        if t in (621,633,645,687,761):
            for same,other,kind in ((row[:6],U,'coefficient'),(U,row[:6],'amplitude')):
                for a,b in combinations(same,2):
                    for y in other:
                        c,x,width,positive=support(a,b,y)
                        check(not positive,'complete native mixed-support control')
                        check(c>Q if kind=='coefficient'else x>Q*width,'independent higher-multiple obstruction')
            for d in range(1,7):
                for indices in combinations(range(13),13-d):
                    chosen=set(indices);c=gcd(*(row[i]for i in chosen))
                    gs=tuple(sorted(gcd(c,row[i])for i in range(13)if i not in chosen))
                    check((c,gs)in bank[d]and c==1,'every full inherited body profile')
        physical.append([t,len(safe),str(phase),str(clearance)])
    hostile=[]
    for t in failures:
        c1=literal_pair_count(101,128,t,1,2);c2=literal_pair_count(113,128,t,1,2)
        check(c1==0 and c2>7*ceildiv(t,7)-t,'actual first-pair half-grid failure and second-pair repair')
        hostile.append([t,c1,c2])
    check([r[2]for r in hostile]==[18,16,14,16,24],'reported actual repair counts')

    manifest=dict(Q=Q,V=V,U=U,scale_counts=dict(generic=84,actual=81),joint_counts=counts,
                  primary_failures=hostile,physical_controls=physical,
                  scope='atleast5 weak-safe points on every translated grid in the stated ratio class; actual81-entry family strictly safe')
    encoded=json.dumps(manifest,sort_keys=True,separators=(',',':'))
    print('PROVED finite-scale joint-position certificate; actual entry and hereditary predicates retained')
    print(encoded)
    print('EXPLICIT_GATES',GATES)
    print('SEMANTIC_SHA256',hashlib.sha256(encoded.encode()).hexdigest())


if __name__=='__main__':main()
