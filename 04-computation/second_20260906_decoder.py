"""Exact actual balanced-entry controls beyond the current unit minimum gate."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
Q=91**6
QS=(179,181,183,185,187)
P=prod(QS)
V=tuple(sorted((P,)+tuple((356-q)*(P//q) for q in QS)))
U=(1,49,109,331,331**2,331**3,331**4)
TS=(36883259177,36883259192)
GATES=0


def need(ok,message):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(message)


def allowed_sum(n):
    if n>356:
        return False
    p=2
    while p*p<=n:
        e=0
        while n%p==0:
            n//=p
            e+=1
        if e and (p%3!=2 or e>2):
            return False
        p+=1
    return n==1 or n%3==2


def graph(row):
    return [{j for j,b in enumerate(row) if j!=i and allowed_sum((a+b)//gcd(a,b))} for i,a in enumerate(row)]


def components(adj):
    unseen=set(range(len(adj)))
    result=[]
    while unseen:
        todo=[min(unseen)]
        seen=set()
        while todo:
            u=todo.pop()
            if u in seen:
                continue
            seen.add(u)
            todo.extend(adj[u]-seen)
        unseen-=seen
        result.append(tuple(sorted(seen)))
    return result


def rank(rows):
    A=[list(map(F,row)) for row in rows]
    r=0
    for c in range(len(A[0])):
        pivot=next((k for k in range(r,len(A)) if A[k][c]),None)
        if pivot is None:
            continue
        A[r],A[pivot]=A[pivot],A[r]
        v=A[r][c]
        A[r]=[a/v for a in A[r]]
        for k in range(len(A)):
            if k!=r and A[k][c]:
                v=A[k][c]
                A[k]=[a-v*b for a,b in zip(A[k],A[r])]
        r+=1
        if r==len(A):
            break
    return r


def crossing(A,B,Y):
    """Independent exact signed-box support test; pair height is checked first."""
    A,B=sorted((A,B))
    D=gcd(A,B)
    a,b=A//D,B//D
    need(b<=Q,"inherited internal-pair height satisfied")
    delta=gcd(D,Y)
    c,x=D//delta,Y//delta
    if c>Q:
        return False
    r=pow(a,-1,b)*x%b
    s=(x-a*r)//b
    ceil=lambda n,d:-((-n)//d)
    lower=max(ceil(-Q-r,b),ceil(s-Q,a))
    upper=min((Q-r)//b,(s+Q)//a)
    return lower<=upper


def clearance_numerator(row,k,d):
    return min(min(v*k%d,d-v*k%d) for v in row)


def main():
    raw=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    inherited_hash=hashlib.sha256(raw).hexdigest()
    need(inherited_hash=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f',"pinned full joint-shadow profiles")
    data=json.loads(raw)
    profiles={int(k):{(c,tuple(word)) for c,word in level['profiles']} for k,level in data['levels'].items()}
    need(V==(185370716505,189592176609,193905909065,198314972625,202822562745,205114343115),"exact six-star cofactors")
    need(reduce(gcd,V)==reduce(gcd,U)==1 and 1 not in V and 1 in U,"unit placement and primitive shapes")
    need(all(gcd(a,b)==1 for a,b in combinations(QS,2)),"pairwise coprime omitted factors")
    need(all(gcd(q,356-q)==1 for q in QS),"each star edge is primitive")
    need(all(allowed_sum((P+v)//gcd(P,v)) for v in V if v!=P),"all five star edges are actual")
    need(all(allowed_sum((a+b)//gcd(a,b)) for a,b in ((1,49),(1,109),(1,331),(331,331**2),(331**2,331**3),(331**3,331**4))),"actual seven-shape spanning tree")
    need(min(V)>3*Q//28,"balanced automatic minimum cutoff fails")
    need(56*max(U)*min(V)>6*Q*(max(U)+1),"larger-unit native comparison also fails with g=1")
    # Positive-real cofactor bound newly proposed by the root, used only as
    # a cross-check on consistency, not as an input to any actual-entry proof.
    cofactor_bound=F(356**5*4**4,5**5)
    need(min(V)<cofactor_bound<355**5,"control is consistent with sharpened six-vertex tree bound")
    reports=[]
    for t in TS:
        core=tuple(t*v for v in V)
        row=core+U
        need(len(set(row))==13 and min(row)>0 and reduce(gcd,row)==1,"distinct primitive positive physical row")
        need(gcd(t,prod(U))==1,"physical scale is coprime to every unit-shape label")
        need(sum(row)<Q**2,"strict actual physical box")
        need(t*min(V)>Q*(U[-1]+U[-2]),"one-core/two-unit-shape dominance")
        min_c=min(t*gcd(a,b)//gcd(t*gcd(a,b),u) for a,b in combinations(V,2) for u in U)
        need(min_c>Q,"all two-core/one-unit-shape cleared coefficients exceed budget")
        adj=graph(row)
        comps=components(adj)
        need(set(comps)=={tuple(range(6)),tuple(range(6,13))},"actual graph has exactly claimed six-plus-seven components")
        edges=[(i,j) for i in range(13) for j in adj[i] if i<j]
        rows=[]
        for i,j in edges:
            d=gcd(row[i],row[j])
            r=[0]*13
            r[i],r[j]=row[j]//d,-row[i]//d
            need(max(map(abs,r))<=Q and sum(a*b for a,b in zip(r,row))==0,"literal bounded decoder edge row")
            rows.append(r)
        need(rank(rows)==11,"independent rational decoder rank")
        count=0
        for A,B in combinations(core,2):
            for Y in U:
                need(not crossing(A,B,Y),"first mixed orientation has no bounded relation")
                count+=1
        for A,B in combinations(U,2):
            for Y in core:
                need(not crossing(A,B,Y),"reverse mixed orientation has no bounded relation")
                count+=1
        need(count==231,"complete two mixed-support orientations")
        maxima={}
        profile_count=0
        for size in range(7,13):
            maximum=0
            for selected in combinations(range(13),size):
                c=reduce(gcd,(row[i] for i in selected))
                maximum=max(maximum,c)
                word=tuple(sorted(gcd(c,row[j]) for j in range(13) if j not in selected))
                need(c==1,"every large-subset gcd equals one")
                need((c,word) in profiles[13-size],"full inherited profile, not just scalar cap")
                profile_count+=1
            maxima[size]=maximum
        need(profile_count==4095,"complete seven-through-twelve subset universe")
        need(6*t<56*max(U) and 7<49*max(V),"both uniform Lipschitz-grid sufficient comparisons fail")
        endpoint_failures=[]
        for u in U[:-1]:
            D=gcd(u,max(U))
            delta=gcd(D,t*min(V))
            A,B=u//D,max(U)//D
            R=Q*(A+B)-(A-1)*(B-1)
            need(delta==1 and D//delta<=Q,"all maximum-endpoint coefficient gates pass")
            need(6*delta*R<56*max(U)*min(V),"all maximum-endpoint native phase comparisons fail")
            endpoint_failures.append({'u':u,'D':D,'delta':delta,'c':D//delta,'R':R})
        need(6*Q*(max(U)+109)<56*max(U)*min(V),"uniform arithmetic upper bound excludes every endpoint partner")
        phase=(1,2) if t==TS[0] else (11,23)
        numerator=clearance_numerator(row,*phase)
        expected=F(1,2) if t==TS[0] else F(3,23)
        need(F(numerator,phase[1])==expected>F(1,14),"literal full-row strict safe phase")
        if t==TS[1]:
            need(t%16==8,"mixed-parity control has dyadic valuation three")
            need(all(clearance_numerator(row,k,16)==1 for k in range(1,16,2)),"all odd sixteenth phases miss target")
            for d in range(2,23):
                for k in range(1,d):
                    need(14*clearance_numerator(row,k,d)<d,"no denominator below twenty-three supplies weak safety")
        reports.append({'t':t,'g':1,'physical_row':row,'sum':sum(row),'component_sizes':list(map(len,comps)),'decoder_rank':11,'decoder_edges':len(edges),'mixed_supports':count,'minimum_cleared_coefficient':min_c,'native_endpoint_failures':endpoint_failures,'subset_profile_count':profile_count,'subset_maxima':maxima,'phase':str(F(*phase)),'clearance':str(expected)})
    print('STATUS: PASS; actual balanced equality entries violate automatic minimum/native unit comparisons')
    print('Q:',Q,'Q_SQUARED:',Q**2)
    print('PRIMITIVE_SIX_SHAPE:',V,'MINIMUM:',min(V),'BALANCED_CUTOFF:',3*Q//28)
    print('PRIMITIVE_SEVEN_SHAPE:',U)
    print('INHERITED_PROFILE_SHA256:',inherited_hash)
    for report in reports:
        print(json.dumps(report,sort_keys=True))
    print('LOGICAL SCOPE: refutes forcing the balanced minimum/native inequalities from actual entry plus all retained gcd profiles; both displayed rows are safe')
    print('OPEN: uniform balanced larger-unit closure and general unit-free actual entry')
    print('ACTIVE GATES:',GATES)


if __name__=='__main__':
    main()
