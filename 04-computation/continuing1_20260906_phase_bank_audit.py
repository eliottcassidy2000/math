"""Independent exact phase-bank referee: literal residues and rational arcs."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, ceil
import sys

sys.stdout.reconfigure(newline="\n")
gates = 0

def need(ok, label):
    global gates
    gates += 1
    if not ok:
        raise ArithmeticError(label)

def distance(x):
    a = x % 1
    return min(a,1-a)

def clearance(shape, phase):
    return min(distance(v*phase) for v in shape)

def gap(points):
    pts = sorted(set(x % 1 for x in points))
    return max([b-a for a,b in zip(pts,pts[1:])]+[pts[0]+1-pts[-1]])

A = {a%23 for v in (1,3,5,6,10,11) for a in (v,-v)}
B = [b for b in range(23) if clearance(A,F(b,23)) >= F(2,23)]
need(B == [3,5,6,10,11,12,13,17,18,20],"complete residue bank")
need(min(clearance(A,F(b,23)) for b in B)==F(2,23),"uniform bank minimum clearance sharp")
need(set(B) == set(range(1,23))-{pow(a,-1,23) for a in A},"independent inverse-residue complement")
profile = [23*gap([F(g*b,23) for b in B]) for g in range(23)]
need(profile == [23,6,4,10,5,5,6,4,5,7,8,6,6,8,7,5,4,6,5,5,10,4,6],"all residue gap profiles")
need(all(profile[g]==profile[-g%23] for g in range(23)),"sign symmetry")
need(gap([F(b,23) for b in B])==F(6,23),"native g1 gap")
print("BANK",B,"clearance",F(2,23),"margin",F(2,23)-F(1,14))
print("GAP_NUMERATORS g0..22",','.join(map(str,profile)))

lift_controls=0
for t in range(1,15):
    for g in range(1,31):
        common=gcd(g,t)
        effective_t,effective_g=t//common,g//common
        images = [(g*(F(b,23)+k)/t)%1 for b in B for k in range(t)]
        m=profile[effective_g%23]
        need(gap(images)==F(m,23*effective_t),"literal full lifted union effective-order gap")
        need(len(set(images))==(effective_t if effective_g%23==0 else len(B)*effective_t),"collapsed residue count")
        lift_controls+=1
need(gap([F(2*(b+23*k),46) for b in B for k in range(2)]) != F(profile[2],46),"noncoprime wrong-gap hostile")

def proper_phase(U):
    # Every weak-safe component has an endpoint from these literal rational walls.
    walls={F(0),F(1)}
    for u in U:
        for k in range(u+1):
            for sign in (-1,1):
                x=F(8*k+sign,8*u)
                if 0<=x<=1:
                    walls.add(x)
    safe=[x for x in sorted(walls) if clearance(U,x)>=F(1,8)]
    need(bool(safe),"literal proper seven-body phase exists")
    return safe[0]

def bank_hit(beta,t,g):
    inv=pow(g,-1,t) if t>1 else 0
    candidates=[]
    for b in B:
        raw=F(g*b,23)
        integral=raw.numerator//raw.denominator
        phase=raw-integral
        z=t*beta-phase
        floor=z.numerator//z.denominator
        for k in (floor,floor+1):
            y=(k+phase)/t
            j=(inv*(k-integral))%t
            x=(j+F(b,23))/t
            need((g*x-y).denominator==1,"literal lift address")
            candidates.append((abs(y-beta),x,y,b,j))
    return min(candidates)

def strictify(V,U,t,g,beta,x,y):
    need(clearance(V,t*x)>=F(2,23),"bank strictly protects V")
    radius=F(3,56*max(U))
    need(abs(y-beta)<=radius,"bank hits closed proper-U arc")
    need(clearance(U,g*x)>=F(1,14),"closed arc gives weak U safety")
    row=[t*v for v in V]+[g*u for u in U]
    if clearance(row,x)>F(1,14):
        return x
    need(beta!=y,"weak endpoint differs from safe center")
    step=min(F(5,322)/(2*t*max(V)),abs(beta-y)/(2*g))
    need(step>0,"positive strictification displacement")
    xp=x+(step if beta>y else -step)
    need(clearance(row,xp)>F(1,14),"literal strict phase after endpoint perturbation")
    return xp

V=tuple(a+23*(100+i) for i,a in enumerate((1,3,5,6,10,11)))
need(gcd(*V)==1 and {v%2 for v in V}=={0,1},"primitive genuinely mixed-parity V control")
constructive_controls=0
for U in combinations(range(1,11),7):
    beta=proper_phase(U)
    for g in (1,2,3,7,10,23):
        m=profile[g%23]
        t=ceil(F(28*m*max(U),69))
        while gcd(t,g)!=1:
            t+=1
        need(28*m*max(U)<=69*t,"exact bank gate")
        dist,x,y,b,j=bank_hit(beta,t,g)
        need(dist<=F(m,46*t),"nearest literal bank hit within half gap")
        phase=strictify(V,U,t,g,beta,x,y)
        row=[t*v for v in V]+[g*u for u in U]
        need(len(set(row))==13 and clearance(row,phase)>F(1,14),"actual thirteen distinct speeds strict-safe")
        constructive_controls+=1

# Equality and a real weak-safe bank endpoint, with a mixed-parity primitive V.
Ve=(1,3,5,6,10,11)
Ue=(1,2,3,4,5,6,23)
te,ge,beta=56,1,F(1,8)
need(56*max(Ue)==23*te,"native bank gate equality")
need(clearance(Ue,beta)==F(1,8),"proper supplier center attained")
xe=beta+F(3,56*23)
need((te*xe)%1==F(3,23),"weak endpoint is an actual bank lift")
need(clearance(Ue,xe)==F(1,14) and clearance(Ve,te*xe)>=F(2,23),"actual endpoint has exactly weak full clearance")
xep=strictify(Ve,Ue,te,ge,beta,xe,xe)
print("EQUALITY V",Ve,"U",Ue,"t",te,"g",ge,"endpoint",xe,"strict phase",xep)

# Incoming actual-entry controls are inherited here only as physical safety controls.
Vi=(185370716505,189592176609,193905909065,198314972625,202822562745,205114343115)
Ui=(1,49,109,331,109561,36264691,12003612721)
need(all(v%23 in A for v in Vi),"incoming V residue hypothesis")
need(all(v%2 for v in Vi+Ui),"both incoming primitive shapes are all odd")
for t in (36883259177,36883259192):
    need(56*max(Ui)<=23*t,"incoming bank gate")
    dist,x,y,b,j=bank_hit(F(1,2),t,1)
    strict=strictify(Vi,Ui,t,1,F(1,2),x,y)
    need(clearance([t*v for v in Vi]+list(Ui),strict)>F(1,14),"incoming bank strict witness")
    # Actual-clearance single-phase control, not only the universal proper-LRC radius.
    k=(t-1)//2
    xhalf=(F(1,2)+k)/t
    need(clearance(Vi,t*xhalf)==F(1,2),"single-phase V actual clearance")
    need(7*max(Ui)<=6*t,"actual U half-phase wider-arc condition")
    need(clearance([t*v for v in Vi]+list(Ui),xhalf)>F(1,14),"incoming row already has actual-clearance single-phase proof")
    print("INHERITED t",t,"bank",strict,"single-phase",xhalf,"clearance",clearance([t*v for v in Vi]+list(Ui),xhalf))

# Separate producer physical control; decoder equality is intentionally not inferred.
t=36883259192
Unew=(2,3,4,5,6,7,max(Ui)-2)
xnew=F(1,8)+F(3,23*t)
need(gcd(*Unew)==1 and min(Unew)>1 and {u%2 for u in Unew}=={0,1},"mixed-parity nonunit primitive U")
need(clearance(Unew,F(1,8))==F(1,8),"new control proper phase")
claimed=F(1,8)-F(3*max(Unew),23*t)
need(clearance([t*v for v in Vi]+list(Unew),xnew)==claimed>F(1,14),"literal new physical-control clearance")
print("MIXED_NONUNIT_CONTROL",Unew,"phase",xnew,"clearance",claimed)

# Final actual 6+7 control, with both primitive components of mixed parity.
Q=91**6
Va=(185370716505,189592176609,193905909065,198314972625,201676672560,205114343115)
Ua=(2,3,6,331,109561,36264691,12003612721)
actual=[t*v for v in Va]+list(Ua)
need(gcd(*Va)==gcd(*Ua)==1 and min(Ua)>1,"new entry primitive unitless shapes")
need({v%2 for v in Va}=={0,1} and {u%2 for u in Ua}=={0,1},"both primitive shapes genuinely mixed")
need(all(v%23 in A for v in Va),"actual entry residue class")
need(gcd(*actual)==1 and len(set(actual))==13 and sum(actual)<Q*Q,"actual physical box")

def atlas(a,b):
    d=gcd(a,b)
    z=a//d+b//d
    if z>356:
        return False
    p=2
    while p*p<=z:
        e=0
        while z%p==0:
            e+=1
            z//=p
        if e and (p%3!=2 or e>2):
            return False
        p+=1
    return z==1 or z%3==2

edges=[]
matrix=[]
for i,j in combinations(range(13),2):
    if atlas(actual[i],actual[j]):
        d=gcd(actual[i],actual[j])
        row=[0]*13
        row[i],row[j]=actual[j]//d,-actual[i]//d
        need(sum(x*n for x,n in zip(row,actual))==0 and max(map(abs,row))<=Q,"literal bounded atlas relation")
        edges.append((i,j))
        matrix.append(row)
seen=set()
components=[]
for seed in range(13):
    if seed in seen:
        continue
    todo=[seed]
    comp=set()
    while todo:
        i=todo.pop()
        if i in comp:
            continue
        comp.add(i)
        todo.extend(j if i==a else a for a,j in edges if i in (a,j))
    seen|=comp
    components.append(comp)
need(components==[set(range(6)),set(range(6,13))],"complete actual atlas components")
rank=0
for col in range(13):
    pivot=next((r for r in range(rank,len(matrix)) if matrix[r][col]),None)
    if pivot is None:
        continue
    matrix[rank],matrix[pivot]=matrix[pivot],matrix[rank]
    for r in range(rank+1,len(matrix)):
        if matrix[r][col]:
            a,b=matrix[rank][col],matrix[r][col]
            matrix[r]=[a*x-b*y for x,y in zip(matrix[r],matrix[rank])]
            content=gcd(*matrix[r])
            if content:
                matrix[r]=[x//content for x in matrix[r]]
    rank+=1
need(rank==11,"literal atlas relation rank")
cleared=[]
for a,b in combinations(Va,2):
    d=t*gcd(a,b)
    for u0 in Ua:
        c=d//gcd(d,u0)
        need(c>Q,"all105 two-core one-outside coefficient obstructions")
        cleared.append(c)
ratios=[]
cleared_ratios=[]
for v0 in Va:
    for a,b in combinations(Ua,2):
        ratio=F(t*v0,Q*(a+b))
        need(ratio>1,"all126 reverse-orientation amplitude obstructions")
        ratios.append(ratio)
        d=gcd(a,b)
        c=d//gcd(d,t*v0)
        cleared_ratios.append(c*ratio)
need(min(cleared)==36446911794129044,"exact minimum cleared coefficient")
need(min(ratios)==F(778001386374677670,778001386051180109),"independent raw amplitude minimum")
need(min(cleared_ratios)==F(3112005545498710680,3102632035080254131),"exact minimum cleared amplitude ratio")
need(sum(actual)==43300016482122890820283,"exact physical sum")
need(clearance(actual,xnew)==F(26509924715,212078740354)>F(1,14),"new actual entry literal strict phase")
need(28*max(Ua)>3*t and 56*max(Ua)<23*t,"generic single supplier misses while bank passes")
maxima=tuple(max(gcd(*(actual[i] for i in subset)) for subset in combinations(range(13),size)) for size in range(12,6,-1))
need(maxima==(1,1,1,1,2,3),"all scalar large-subset gcd maxima; no profile claim")
print("ACTUAL_MIXED_6+7 rank",rank,"edges",len(edges),"sum",sum(actual),"phase",xnew,"clearance",clearance(actual,xnew))

odd12=[b for b in range(12) if all(distance(F(a*b,12))>F(1,14) for a in range(1,12,2))]
need(odd12==[1,2,3,5,6,7,9,10,11] and gap(F(b,12) for b in odd12)==F(1,6),"odd-only12 bank beats23 on odd shapes")
producer_boundary=F(41,322)
producer_epsilon=F(1,2958897468039744)
producer_endpoint_row=[56*v for v in Vi]+[2,3,4,5,6,7,23]
need(clearance(producer_endpoint_row,producer_boundary)==F(1,14),"producer large equality endpoint")
need(clearance(producer_endpoint_row,producer_boundary-producer_epsilon)==F(9189122571553,128647716001728)>F(1,14),"producer exact endpoint perturbation")
for ts in range(1,100):
    if ts%7==0:
        continue
    k=(3*pow(ts,-1,7))%7
    z=F(ts*k,7)
    j=z.numerator//z.denominator
    x=(j+F(1,2))/ts
    need(abs(x-F(k,7))==F(1,14*ts),"inherited seven-clock exact address")
    need(all(distance(u*x)>F(1,14) for u in range(1,ts) if u%7),"seven-clock full bounded unit U control")

print("UNIVERSE",lift_controls,"literal lifted gap controls;",constructive_controls,"new mixed-parity constructive controls")
print("PASS",gates,"always-active exact gates; no producer imports or floating point")
