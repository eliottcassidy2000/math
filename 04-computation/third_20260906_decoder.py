"""Rooted gcd products + inherited leaf bounds force actual decoder closure."""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, prod
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
Q=91**6
MINIMUM_BOUNDS={2:177,3:31684,4:6684150,5:1694040507,6:468424857663}
PAIR_BOUNDS={2:1,3:177,4:35483,5:8350122,6:2170260903}
HEIGHT_BOUNDS={2:13520696477,3:124998734,4:914513,5:5397,6:28}
SCALE_CAPS={11:2,10:4,9:9,8:30,7:90}
GATES=0


def need(ok,message):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(message)


def gcd_product_controls():
    cases=0
    rooted=0
    for n in range(2,7):
        for w in product(range(1,7),repeat=n):
            h=reduce(gcd,w)
            for i,x in enumerate(w):
                pairgcds=[gcd(x,y) for j,y in enumerate(w) if j!=i]
                P=prod(pairgcds)
                need((h*x**(n-2))%P==0,"exact rooted gcd-product divisibility at arbitrary common scale")
                if h==1:
                    need(min(pairgcds)**(n-1)<=x**(n-2),"minimum partner from product bound")
                    if n>=3 and x>1:
                        need(min(pairgcds)**(n-1)<x**(n-2),"strict primitive boundary away from a unit")
                rooted+=1
            cases+=1
    primes=(2,3,5,7,11,13,17)
    for n in range(3,8):
        P=prod(primes[:n])
        w=tuple(P//p for p in primes[:n])
        need(reduce(gcd,w)==1,"complementary-prime equality control is primitive")
        for i,x in enumerate(w):
            need(prod(gcd(x,y) for j,y in enumerate(w) if j!=i)==x**(n-2),"exact rooted product equality")
    need(gcd(2,4)*gcd(2,6)>2,"erasing common scale fails at the first three-label hostile")
    return cases,rooted


def cutoffs():
    rows=[]
    for a,M in MINIMUM_BOUNDS.items():
        b=13-a
        D=PAIR_BOUNDS[a]
        if a==2:
            need(D==1,"primitive two-label gcd")
        else:
            need(D**(a-1)<M**(a-2)<=(D+1)**(a-1),"exact strict-power integer pair bound")
        H=HEIGHT_BOUNDS[a]
        need(H==min(Q//SCALE_CAPS[b],a*Q//(7*(b+1)*D)),"exact automatic larger-height ceiling")
        need(SCALE_CAPS[b]*H<=Q and 7*(b+1)*D*H<=a*Q,"both uniform native gates are paid")
        old=a*Q//(7*(b+1)*355**(a-2))
        need(H>=old,"comparison with inherited outgoing-edge envelope")
        rows.append((a,b,M,D,H,old))
    return rows


def inert_sums():
    prime=[True]*357
    prime[0]=prime[1]=False
    for p in range(2,19):
        if prime[p]:
            for k in range(p*p,357,p):
                prime[k]=False
    result={1}
    for p in range(2,357):
        if prime[p] and p%3==2:
            result={s*p**e for s in result for e in range(3) if s*p**e<=356}
    return result


def graph(row,allowed):
    adj=[set() for _ in row]
    for i,j in combinations(range(len(row)),2):
        if (row[i]+row[j])//gcd(row[i],row[j]) in allowed:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def components(adj):
    unseen=set(range(len(adj)))
    out=[]
    while unseen:
        seen={min(unseen)}
        queue=list(seen)
        for i in queue:
            for j in sorted(adj[i]-seen):
                seen.add(j)
                queue.append(j)
        unseen-=seen
        out.append(tuple(sorted(seen)))
    return out


def matrix_rank(rows):
    A=[list(map(F,row)) for row in rows]
    r=0
    for col in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][col]),None)
        if pivot is None:
            continue
        A[r],A[pivot]=A[pivot],A[r]
        value=A[r][col]
        A[r]=[v/value for v in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][col]:
                value=A[i][col]
                A[i]=[v-value*w for v,w in zip(A[i],A[r])]
        r+=1
        if r==len(A):
            break
    return r


def crossing(A,B,Y):
    D=gcd(A,B)
    a,b=sorted((A//D,B//D))
    need(b<=Q,"every internal primitive pair fits the inherited budget")
    delta=gcd(D,Y)
    c,x=D//delta,Y//delta
    if c>Q:
        return False
    r=pow(a,-1,b)*x%b
    s=(x-a*r)//b
    ceil=lambda n,d:-((-n)//d)
    return max(ceil(-Q-r,b),ceil(s-Q,a))<=min((Q-r)//b,(s+Q)//a)


def actual_control():
    qs=(215,251,257,263,273)
    P=prod(qs)
    V=tuple(sorted((P,)+tuple(83*P//q for q in qs)))
    U=(2,3,4,8,16,19,28)
    t,g=3251,1
    row=tuple(t*v for v in V)+U
    need(V==(302746510145,314257784295,321594541905,329282060835,384417661719,995780689995),"exact high-minimum cofactor shape")
    need(len(set(row))==13 and reduce(gcd,row)==reduce(gcd,V)==reduce(gcd,U)==1,"primitive distinct physical row and shapes")
    need(1 not in V and 1 not in U,"both shapes are unitless")
    need(all(gcd(a,b)==1 for a,b in combinations(qs,2)) and all(gcd(83,q)==1 for q in qs),"primitive common-numerator star")
    need(3*Q//28<min(V)<=MINIMUM_BOUNDS[6] and max(U)==HEIGHT_BOUNDS[6],"control reaches new height ceiling while old minimum gate fails")
    need(sum(row)==8608905638154474<Q**2,"physical box retained")
    allowed=inert_sums()
    atlas=[(a,b) for a in range(1,178) for b in range(a+1,357-a) if gcd(a,b)==1 and a+b in allowed]
    need(len(atlas)==5855,"strict actual atlas rebuilt")
    adj=graph(row,allowed)
    need(components(adj)==[tuple(range(6)),tuple(range(6,13))],"literal actual balanced graph")
    center=V.index(P)
    need(len(adj[center])==5 and all(len(adj[i])==1 for i in range(6) if i!=center),"actual six-graph is exactly the common-numerator star")
    numerators=tuple(min(P,w)//gcd(P,w) for w in V if w!=P)
    need(numerators==(83,)*5 and all(gcd(a,b)>1 for a,b in combinations(numerators,2)),"incoming coprime-numerator star subclass does not apply")
    edges=[(i,j) for i in range(13) for j in adj[i] if i<j]
    relations=[]
    for i,j in edges:
        d=gcd(row[i],row[j]);r=[0]*13
        r[i],r[j]=row[j]//d,-row[i]//d
        need(max(map(abs,r))<=355 and sum(a*b for a,b in zip(r,row))==0,"literal bounded decoder edge")
        relations.append(r)
    need(matrix_rank(relations)==11,"independent weighted decoder rank")
    min_factor=min(gcd(a,b)//gcd(gcd(a,b),u) for a,b in combinations(V,2) for u in U)
    need(gcd(t,prod(U))==1 and min_factor==174684705,"retained distinguished divisor")
    need(t*min_factor>Q,"all two-V one-U distinguished coefficients exceed Q")
    need(t*min(V)>Q*(U[-1]+U[-2]),"all one-V two-U supports fail by amplitude")
    mixed=0
    for inside,outside in ((row[:6],row[6:]),(row[6:],row[:6])):
        for A,B in combinations(inside,2):
            for Y in outside:
                need(not crossing(A,B,Y),"independent complete mixed signed-box test")
                mixed+=1
    need(mixed==231,"both complete mixed orientations")
    raw=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    profile_hash=hashlib.sha256(raw).hexdigest()
    need(profile_hash=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f',"pinned inherited full profiles")
    data=json.loads(raw)
    profiles={int(k):{(c,tuple(word)) for c,word in level['profiles']} for k,level in data['levels'].items()}
    profile_count=0
    for size in range(7,13):
        for selected in combinations(range(13),size):
            d=reduce(gcd,(row[i] for i in selected))
            word=tuple(sorted(gcd(d,row[j]) for j in range(13) if j not in selected))
            need(d==1 and (d,word) in profiles[13-size],"complete seven-through-twelve profile survives")
            profile_count+=1
    need(profile_count==4095,"full retained subset universe")
    v=min(V)
    pair=(v,min(V[1:],key=lambda w:gcd(v,w)))
    D=gcd(*pair)
    need(D==1151127415 and D<=PAIR_BOUNDS[6],"explicit successful partner of the smallest label")
    need(prod(gcd(v,w) for w in V if w!=v)==v**4,"primitive product lemma has exact equality on the control")
    need(56*D*max(U)<=6*Q and 56*max(pair)>=6*90,"inherited dual selected-pair theorem applies")
    need(6*Q*(U[-1]+U[-2])<56*max(U)*min(V),"all opposite-orientation native phase inequalities fail")
    for a,b in combinations(U,2):
        d=gcd(a,b)
        for w in V:
            need(d*w//gcd(d,w)>3*Q//28,"all cross-divisor scores fail before any local-to-global loss")
    x=F(1301,6502)
    clearance=min(min((n*x)%1,(-n*x)%1) for n in row)
    need(clearance==F(1289,6502)>F(1,14),"literal full-row strict safe time")
    return {'V':V,'U':U,'t':t,'g':g,'sum':sum(row),'decoder_edges':len(edges),'mixed_supports':mixed,'profile_count':profile_count,'all_large_subset_gcds':1,'selected_pair':pair,'selected_pair_gcd':D,'minimum_cleared_factor':min_factor,'phase':str(x),'clearance':str(clearance),'inherited_profile_sha256':profile_hash}


def main():
    cases,rooted=gcd_product_controls()
    table=cutoffs()
    sharp=(177*178,177*179,178*179)
    need(reduce(gcd,sharp)==1 and min(gcd(a,b) for a,b in combinations(sharp,2))==177,"sharp actual three-component pair-gcd control")
    need(len(components(graph(sharp,inert_sums())))==1,"sharp three-component control has actual connected atlas")
    control=actual_control()
    print('STATUS: PASS; rooted gcd divisibility forces uniform automatic decoder classes')
    print('FINITE_ROOTED_PRODUCT_UNIVERSE:',cases,'positive tuples of lengths2..6 over1..6;',rooted,'distinguished-root cases')
    for a,b,M,D,H,old in table:
        print('AUTOMATIC_CLASS:',a,b,'minimum_bound',M,'forced_pair_gcd',D,'maxU',H,'old_maxU',old)
    print('SHARP_THREE_COMPONENT:',sharp,'pair_gcds',(177,178,179))
    print('ACTUAL_BOUNDARY_CONTROL:',json.dumps(control,sort_keys=True))
    print('MECHANISM: primitive rooted divisor-product bound plus inherited leaf cofactors and dual coefficient escape')
    print('SCOPE: automatic balanced maxU<=28; universal decoder entry and LRC14 remain OPEN')
    print('ACTIVE GATES:',GATES)


if __name__=='__main__':
    main()
