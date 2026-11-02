#!/usr/bin/env python3
"""Finite exact controls, not a finite-computation proof of infinite results."""
import hashlib, itertools, json, pathlib

CHECKS=0
def check(p,msg):
    global CHECKS
    CHECKS+=1
    if not p: raise RuntimeError(msg)

def C(n): return n//2 if n%2==0 else 3*n+1

def U(n):
    v=3*n+1;a=0
    while v%2==0: v//=2;a+=1
    return v,a

def predecessors(m):
    ans=[]
    for a in range(1,7):
        if (pow(2,a)*m)%9 in (4,7):
            n=(pow(2,a)*m-1)//3
            check(n>0 and n%2==1 and n%3!=0,'odd inverse domain')
            check(U(n)==(m,a),'exact odd inverse')
            v=n
            for _ in range(a+1): v=C(v)
            check(v==m,'standard clock')
            ans.append(n)
    check(len(set(ans))==2,'two branches')
    return ans

# Six complete residue classes certify the finite residue statement.
residue_table={str(r):[a for a in range(1,7) if pow(2,a,9)*r%9 in (4,7)]
               for r in (1,2,4,5,7,8)}
for x in residue_table.values(): check(len(x)==2,'complete mod9 classes')
check({pow(2,a,9) for a in range(1,7)}=={1,2,4,5,7,8},'order six')
level={1}; depth_receipt=[]
for k in range(1,13):
    nxt=set()
    for m in level:
        kids=predecessors(m)
        check(not(nxt&set(kids)),'disjoint fibers')
        nxt.update(kids)
    level=nxt
    check(len(level)==2**k,'level cardinality')
    depth_receipt.append({'depth':k,'vertices':len(level),'standard_radius_bound':7*k})

# Every map [n]->[n], n<=5: underlying simple graph is 2-degenerate.
functional_maps=0
for n in range(1,6):
    for f in itertools.product(range(n),repeat=n):
        edges={tuple(sorted((i,j))) for i,j in enumerate(f) if i!=j}
        left=set(range(n))
        while left:
            degrees={v:sum(v in e and set(e)<=left for e in edges) for v in left}
            v=min(degrees,key=degrees.get)
            check(degrees[v]<=2,'functional degeneracy')
            left.remove(v)
        functional_maps+=1
check({tuple(sorted((n,C(n)))) for n in (1,2,4)}=={(1,2),(1,4),(2,4)},'Collatz triangle')

# Eisenstein convention omega^2+omega+1=0, Q=a^2-ab+b^2.
units={(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)}
def norm(v): a,b=v;return a*a-a*b+b*b
def c7(v): a,b=v;return (a+2*b)%7
def c3(v): a,b=v;return (a+b)%3
check(all(norm(u)==1 for u in units),'six exact units')
check({c7(u) for u in units}==set(range(1,7)),'quotient K7')
check(all(c3(u)!=0 for u in units),'source proper 3 colors')
cycle_steps=((1,0),(1,0),(0,-1))
cycle_sum=tuple(map(sum,zip(*cycle_steps)))
check(cycle_sum==(2,-1) and norm(cycle_sum)==7,'nonzero kernel displacement')
check(c7(cycle_sum)==0,'modular cycle closure')
check([c7(u) for u in cycle_steps]==[1,1,5],'bad lifted triangle edge labels')
# Exhaustively distinguish all rooted K7 triangles and their uniquely lifted steps.
byres={c7(u):u for u in units}; closed=0; nonclosed=0
for a,b in itertools.permutations(range(1,7),2):
    ds=[byres[a],byres[(b-a)%7],byres[-b%7]]
    sm=tuple(map(sum,zip(*ds)))
    check(c7(sm)==0,'all quotient triangles close modulo ideal')
    if sm==(0,0): closed+=1
    else: nonclosed+=1
check(closed==12 and nonclosed==18,'complete triangle lift census')
# Exact triangular ball check in coordinates ±(1,0),±(0,1),±(1,1).
ball={(0,0)}
for r in range(1,13):
    ball|={(a+x,b+y) for a,b in ball for x,y in units}
    check(len(ball)==1+3*r*(r+1),'triangular ball formula')
check(2**16>1+3*112*113,'unit routing contradiction threshold')

# Independent small-field exact replay of the spectral count and polygon base.
finite_field_receipts=[]
for p in (5,7,11,13):
    def eta(x):
        x%=p
        if x==0:return 0
        return 1 if pow(x,(p-1)//2,p)==1 else -1
    d=next(x for x in range(2,p) if eta(x)==-1)
    field=list(itertools.product(range(p),repeat=2))
    def N(z):a,b=z;return (a*a-d*b*b)%p
    torus=[z for z in field if N(z)==1]
    check(len(torus)==p+1,'norm torus degree')
    for c in range(1,p):
        norm_shell=[z for z in field if N(z)==c]
        for t in range(p):
            tc=sum((2*z[0])%p==t for z in norm_shell)
            kc=sum((x+c*pow(x,-1,p))%p==t for x in range(1,p))
            check(tc==1-eta(t*t-4*c),'torus trace exact distribution')
            check(kc==1+eta(t*t-4*c),'Kloosterman trace exact distribution')
            check(tc+kc==2,'coefficient proof of minus Kloosterman')
    counts={(0,0):1};Ws=[];R3=None
    for k in range(1,5):
        nxt={}
        for (a,b),v in counts.items():
            for x,y in torus:
                z=((a+x)%p,(b+y)%p);nxt[z]=nxt.get(z,0)+v
        counts=nxt;Ws.append(counts.get((0,0),0))
        if k==3:R3=counts
    check(Ws[1]==p+1,'W2 formula')
    check(Ws[2]==(p+1)*(1-eta(-3)),'W3 formula')
    check(Ws[3]==3*p*(p+1),'W4 formula')
    smooth=0
    for c in range(1,p):
        if c in (1,9%p):continue
        vals={R3.get(z,0) for z in field if N(z)==c}
        quartic=p+2+sum(eta((t*t-4*c)*(t-c-1)*(t-c+3)) for t in range(p))
        check(vals=={quartic},'R3 is smooth quartic point count')
        check(quartic%6==0,'six divides smooth R3')
        smooth+=1
    finite_field_receipts.append({'p':p,'nonsquare':d,'closed_walks_1_to_4':Ws,'smooth_quartics':smooth})

# Corrected conjugation-fixed unramified residue classes mod55.
fixed=[];even_not_fixed=[]
for r in range(1,55):
    if __import__('math').gcd(r,55)!=1: continue
    subgroup=[];x=1
    while True:
        x=x*r%55;subgroup.append(x)
        if x==1:break
    if 54 in subgroup:fixed.append([r,len(subgroup)])
    elif len(subgroup)%2==0:even_not_fixed.append(r)
check(fixed==[[19,10],[24,10],[29,10],[39,10],[54,2]],'correct fixed-prime classes')
check(pow(2,10,55)==34 and 2 in even_not_fixed,'false degree20 hostile')

# Exact sign/nullness demonstration for Q(sqrt5) square oracle repair.
# (0+1*sqrt5)^2=5+0*sqrt5; rational-only roots would miss it.
check((0*0+5*1*1,2*0*1)==(5,0),'irrational square root hostile')

root=pathlib.Path(__file__).resolve().parent
pins=json.loads((root/'SOURCE_PINS.json').read_text());source_total=0
for pin in pins:
    loc=pin.get('local_source','sources/'+pathlib.Path(pin['path']).name)
    data=(root/loc).read_bytes()
    check(hashlib.sha256(data).hexdigest()==pin['sha256'],'source SHA256 '+pin['path'])
    check(len(data)==pin['size'],'source size '+pin['path'])
    source_total+=len(data)
# Recover the actual arithmetic-obstructed Gamma9 as a low-chromatic graph.
from fractions import Fraction
jpath=root/'sources'/'05-knowledge__results__cubic_descent_audit_20261016e.json'
hub=json.loads(jpath.read_text())['even_extension_hub']
coords=[tuple(map(Fraction,row)) for row in hub['coordinates']]
actual=[]
for i,j in itertools.combinations(range(9),2):
    a,b,c,d=(x-y for x,y in zip(coords[i],coords[j]))
    nd=(a*a+7*b*b+c*c+7*d*d,2*a*b+2*c*d)
    if nd==(1,0):actual.append([i,j])
    check(any(v!=0 for v in (a,b,c,d)),'Gamma9 injectivity')
check(actual==hub['edges'],'Gamma9 exact unit and nonunit distances')
check(all(hub['three_colouring'][i]!=hub['three_colouring'][j] for i,j in actual),'Gamma9 proper3')
# h-a-w_ab-w_cd-d-h is a 5-cycle.
oddcycle=[0,1,5,8,4]
check(all(sorted((oddcycle[t],oddcycle[(t+1)%5])) in actual for t in range(5)),'Gamma9 odd cycle')

print(json.dumps({'status':'PASS' ,'scope':'Finite exact controls; infinite claims have written proofs; inherited Selmer/modularity proofs not replayed',
 'checks':CHECKS,'functional_maps_exhausted_n_le_5':functional_maps,
 'odd_inverse_residue_table':residue_table,'inverse_levels':depth_receipt,
 'rooted_quotient_triangles':{'closed_exactly':closed,'closed_only_modulo_ideal':nonclosed},
 'Gamma9':{'vertices':9,'edges':len(actual),'chromatic_number':3,'coordinate_field':'Q(sqrt7,i)'},'finite_field_replay':finite_field_receipts,'conjugation_fixed_classes':fixed,'source_files':len(pins),'source_bytes':source_total},indent=2,sort_keys=True))
