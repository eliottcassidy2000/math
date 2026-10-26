"""Recompute 81/27 unnamed response laws and the cyclotomic chirality controls.
Only canonical labelling is imported from the separately audited pinned tcanon.
"""
from pathlib import Path
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import combinations,permutations
import sys,json
from audit import need,paths,partner,arc,relabel
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'04-computation/experiments/joint_port_bisimulation_20261026/audit_hyp9416'))
from tcanon import canon as ircanon
def build(word):
    r=[0]*12
    for b,minus in enumerate(word):
        for i in range(3):r[4*b+i]|=1<<(4*b+(i+1)%3)
        for i in range(3):
            if minus:r[4*b+3]|=1<<(4*b+i)
            else:r[4*b+i]|=1<<(4*b+3)
        for i in range(4):r[4*b+i]|=15<<(4*((b+1)%3))
    return tuple(r)
def poly_mul(a,b):
    v=0
    while b:
        if b&1:v^=a
        a<<=1;b>>=1
    return v
def gf_mul(a,b,mod,deg):
    v=poly_mul(a,b)
    for k in range(v.bit_length()-1,deg-1,-1):
        if v>>k&1:v^=mod<<(k-deg)
    return v
def main():
    controls=[];supports=[];families=[];trip=[]
    for w in [(0,0,0),(0,0,1)]:
        t=build(w);aut=ircanon(t)[1];pres={};multip=Counter();phase=Counter();joint=Counter();H=0
        for p in paths(t):
            H+=1;i=(p[0]//4-2)%3;j=(p[-1]//4-2)%3;phase[i]+=1;joint[i,j]+=1
            code=sum(arc(t,p[a],p[b])<<k for k,(a,b) in enumerate(( (a,b) for a in range(12) for b in range(a+2,12) )))
            multip[code]+=1;pres.setdefault(code,p)
        need(set(multip.values())=={aut},'free path action')
        support=set();fam=defaultdict(set)
        for p in pres.values():
            c=ircanon(partner(t,p))[0];support.add(c);fam[(p[0]//4-2)%3].add(c)
        need(len(support)==len(pres),'81/27 partner injectivity')
        tri={v for v in combinations(range(12),3) if all(sum(arc(t,a,b) for b in v)==1 for a in v)};trip.append(tri)
        supports.append(support);families.append(fam)
        controls.append({'word':w,'H':H,'Aut':aut,'presentations':len(pres),'target_classes':len(support),'entrance_phase_mass':[phase[i] for i in range(3)],'joint_phase_mass':[[joint[i,j] for j in range(3)] for i in range(3)],'phase_target_sizes':[len(fam[i]) for i in range(3)] if any(w) else None})
    need(trip[0]==trip[1] and len(trip[0])==67,'labelled C3 equality')
    need([r['Aut'] for r in controls]==[81,27] and [r['H'] for r in controls]==[462915]*2,'hostile counts')
    shared=supports[0]&supports[1];need(len(shared)==74,'shared supports')
    need(all(not families[1][i]&families[1][j] for i in range(3) for j in range(i+1,3)),'distinct period-three phase supports')
    sh=[len(shared&families[1][i]) for i in range(3)];need(sh==[74,0,0],'shared support phase')
    tv=1-Fraction(len(shared),max(map(len,supports)))
    # Phi_7 mod2 splits into reciprocal cubics; Phi_5 mod2 is irreducible.
    need(poly_mul(0b1011,0b1101)==0b1111111,'Phi7 factorization')
    need(all(gf_mul(a,b,0b11111,4)!=0 for a in range(1,16) for b in range(1,16)),'F16 no zero divisors')
    phases=[]
    for p in [3,5,7,11,13,17,19,23,31]:
        orb=[];v=1
        while v not in orb:orb.append(v);v=2*v%p
        d=len(orb);possible=(p-1 not in orb);need(possible==(d%2==1),'odd-order chirality criterion')
        count=0
        for bits in range(1<<((p-1)//2)):
            D={i if bits>>(i-1)&1 else p-i for i in range(1,(p+1)//2)}
            count+=({2*x%p for x in D}==D)
        expected=2**((p-1)//(2*d)) if possible else 0;need(count==expected,'all inverse selectors checked')
        phases.append({'p':p,'doubling_orbit':orb,'degree':d,'inverse_selector_possible':possible,'selectors':count})
    # The F8 orbit {1,2,4} is the exponent form of the Paley orientation.
    t=tuple(sum(1<<j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7))
    need(sum(1 for _ in paths(t))==189,'Paley H')
    aut=sum(relabel(t,p)==t for p in permutations(range(7)));need(aut==21,'Paley automorphisms')
    # Frobenius normal basis, then the phase cocycle for all words and shifts.
    F=lambda x:gf_mul(x,x,0b1011,3)
    beta=[3,F(3),F(F(3))];need(beta==[3,5,7],'normal basis')
    def ff(x,k):
        for _ in range(k%3):x=F(x)
        return x
    for word in range(8):
        x=0
        for i in range(3):
            if word>>i&1:x^=beta[i]
        for i in range(3):
            for delta in range(3):need(ff(ff(x,-i),-delta)==ff(x,-i-delta),'Frobenius phase transport')
    print(json.dumps({'81_27':controls,'cyclic_triples':67,'shared_targets':74,'shared_by_minus_phase':sh,'unnamed_target_law_total_variation':str(tv),'cyclotomic_chirality':phases,'F8':{'normal_basis':beta,'Paley_H':189,'Paley_Aut':aut,'phase_cocycle_checks':72}},sort_keys=True,indent=2))
if __name__=='__main__':main()
