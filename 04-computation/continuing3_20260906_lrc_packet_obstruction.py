"""Actual decoder entry with a normalized six-component first packet at29.
No external modules or producer imports. Full inherited profiles read as data.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd,prod,lcm
from pathlib import Path
import argparse,hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
N=0

def check(ok,label):
    global N
    N+=1
    if not ok:raise ArithmeticError(label)

def inert_sum(s):
    if s>356:return False
    p=2
    while p*p<=s:
        e=0
        while s%p==0:s//=p;e+=1
        if e and (p%3!=2 or e>2):return False
        p+=1
    return s==1 or s%3==2

def edges(row):return tuple((i,j) for i,j in combinations(range(len(row)),2) if inert_sum((row[i]+row[j])//gcd(row[i],row[j])))

def components(n,es):
    todo=set(range(n));out=[]
    while todo:
        c={min(todo)}
        while True:
            old=len(c)
            for i,j in es:
                if i in c or j in c:c.update((i,j))
            if len(c)==old:break
        todo-=c;out.append(tuple(sorted(c)))
    return tuple(out)

def packet(row,d):return tuple(k for k in range(d) if all(14*min(v*k%d,(-v*k)%d)>=d for v in row))
def clearance(row,x):return min(min((v*x)%1,(-v*x)%1) for v in row)

def main():
    parser=argparse.ArgumentParser()
    here=Path(__file__).resolve()
    default=here.parents[1] if here.parent.name=='04-computation' else Path('C:/w/s0905')
    parser.add_argument('--repo',type=Path,default=default)
    args=parser.parse_args()
    Q=91**6;q=(16,25,27,161,187,247);P=prod(q)
    V=tuple(sorted(P//r for r in q));U=(1,43,67,79,97,109,127);t=412164
    row=tuple(t*v for v in V)+U
    check(P==lcm(*range(1,29))==80313433200,'exact distributed covering product')
    check(all(gcd(a,b)==1 for a,b in combinations(q,2)),'pairwise coprime cofactor groups')
    check(V==(325155600,429483600,498841200,2974571600,3212537328,5019589575),'literal normalized shape')
    check(gcd(*V)==gcd(*U)==gcd(*row)==1 and len(set(row))==13,'primitive distinct physical row')
    check(min(V)>1 and set(v%2 for v in V)=={0,1},'small component is unitless mixed parity')
    check(sum(row)==5135637177376615<Q*Q,'actual bounded physical box')
    ve,ue,re=edges(V),edges(U),edges(row)
    check(components(6,ve)==(tuple(range(6)),),'full normalized six-component atlas connected')
    check(components(7,ue)==(tuple(range(7)),),'full normalized seven-component atlas connected')
    check(components(13,re)==(tuple(range(6)),tuple(range(6,13))),'actual atlas exactly6+7')
    check(len(re)==len(ve)+len(ue)==14,'every actual graph edge retained')
    check(13-len(components(13,re))==11,'actual weighted edge-span rank')
    for v in V:
        for u in U:check(gcd(t*v,u)==1,'literal cross-coprimality')
    dmin=min(gcd(v,w) for v,w in combinations(V,2))
    usum=max((u+w)//gcd(u,w) for u,w in combinations(U,2))
    check(dmin==1738800 and usum==236,'uniform orientation extrema')
    tmin=max(Q//dmin+1,Q*usum//min(V)+1)
    tmax=(Q*Q-sum(U))//sum(V)
    check(tmin==412164 and tmax==25880486141010,'exact bounded actual scale class endpoints')
    check(tmin*dmin>Q and tmin*min(V)>Q*usum,'uniform entire-class relation margins')
    check((tmin-1)*min(V)<=Q*usum,'lower endpoint is exact for this uniform amplitude gate')
    check(tmax*sum(V)+sum(U)<=Q*Q<(tmax+1)*sum(V)+sum(U),'upper endpoint exact physical box')
    check(prod(U)==305613336829,'finite excluded prime modulus for family')
    print('ACTUAL_SCALE_CLASS',tmin,'<=t<=',tmax,'gcd(t,305613336829)=1; primitive shapes fixed')

    coeff=[];ratios=[]
    for v,w in combinations(V,2):
        D=t*gcd(v,w)
        for u in U:
            c=D//gcd(D,u);coeff.append(c)
            check(c>Q,'two V labels: outside coefficient above Q')
    for u,w in combinations(U,2):
        D=gcd(u,w)
        for v in V:
            delta=gcd(D,t*v);r=F(t*v//delta,Q*((u+w)//D));ratios.append(r)
            check(r>1,'two U labels: cleared target above complete signed amplitude')
    check(len(coeff)==105 and len(ratios)==126,'all231 mixed support orientations')
    check(min(coeff)==716670763200 and min(ratios)==F(4786336882800,4786326552917),'exact entry margins')
    # No bounded crossing: internal actual edges span both weighted kernels;
    # the complete decoder criterion gives W=V_dec, not rank alone.
    rel=Path('05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json')
    # The inherited pin is the canonical LF blob; Windows checkout may use CRLF.
    raw=(args.repo/rel).read_bytes().replace(b'\r\n',b'\n')
    pin=hashlib.sha256(raw).hexdigest()
    check(pin=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','inherited full-profile data pin')
    bank={int(d):{(c,tuple(gs)) for c,gs in level['profiles']} for d,level in json.loads(raw)['levels'].items()}
    for d in range(1,7):check((1,(1,)*d) in bank[d],'all-one signature explicitly retained')
    bodies=0
    for tails in range(1,7):
        for ids in combinations(range(13),13-tails):
            c=gcd(*(row[i] for i in ids));complement=set(range(13))-set(ids)
            gs=tuple(sorted(gcd(c,row[i]) for i in complement))
            check(c==1,'every large actual body gcd one')
            check((c,gs) in bank[tails],'every inherited joint-shadow profile retained')
            bodies+=1
    print('ACTUAL_ROW','t',t,'V',V,'U',U,'sum',sum(row))
    print('ATLAS_EDGES_V',ve,'ATLAS_EDGES_U',ue,'rank',11)
    print('ENTRY_MARGINS','min_distinguished',min(coeff),'min_cleared_amplitude',min(ratios),'mixed_supports',231)
    print('ALL_PROFILE_BODIES',bodies,'gcd one; all six all-one signatures retained','JSON_SHA256',pin)
    owners=[]
    for d in range(2,29):
        ds=tuple(i for i,v in enumerate(V) if v%d==0)
        check(bool(ds),'divisibility zero owner for every denominator2..28')
        check(packet(V,d)==() and packet(row,d)==(),'complete primitive and physical small packets empty')
        owners.append((d,ds[0]))
    pv=packet(V,29);pu=packet(U,29);pr=packet(row,29)
    check(pv==(5,6,7,9,10,11,12,17,18,19,20,22,23,24),'complete first primitive packet29')
    check(pr==(5,12,17,24),'complete first physical packet29')
    # Same-denominator unit multiplier join is the simplest inherited CRT
    # compatibility check. Retain the exact residue t, not just packet size.
    check(gcd(t,29)==1 and t%29==16,'effective mod29 multiplier')
    joined=tuple(k for k in range(29) if (t*k)%29 in pv and k in pu)
    check(joined==pr,'literal component-packet compatibility joins physical packet')
    check(clearance(row,F(5,29))==F(3,29)>F(1,14),'small exact physical certificate')
    check(clearance(row,F(166010,332021))==F(7,29),'separate lifted certificate')
    # The same shapes and inherited profiles do not force positive compatibility
    # at their first useful V denominator. A19-unit scale perturbation kills it.
    counts=tuple(sum((mu*k)%29 in pv for k in pu) for mu in range(29))
    check(counts==(0,8,8,2,8,2,0,8,6,4,2,6,2,4,10,10,4,2,6,2,4,6,8,0,2,8,2,8,8),'complete multiplier join profile')
    t2=t+19;row2=tuple(t2*v for v in V)+U
    check(t2==412183 and t2%29==6,'located scale perturbation')
    check(gcd(t2,prod(U))==1 and all(gcd(v,u)==1 for v in V for u in U),'second row cross-coprimality')
    check(sum(row2)<Q*Q and len(set(row2))==13 and gcd(*row2)==1,'second row physical box')
    check(t2*dmin>Q and t2*min(V)>Q*usum,'second row all231 supports excluded by uniform extrema')
    check(components(13,edges(row2))==(tuple(range(6)),tuple(range(6,13))),'second actual atlas exactly6+7')
    for tails in range(1,7):
        for ids in combinations(range(13),13-tails):
            c=gcd(*(row2[i] for i in ids));gs=(1,)*tails
            check(c==1 and (c,gs) in bank[tails],'second row every inherited profile all-one')
    for den in range(2,31):check(packet(row2,den)==(),'second complete physical denominator exclusion')
    check(counts[t2%29]==0 and packet(row2,29)==(),'nonempty primitive packets have zero compatibility')
    check(packet(row2,31)==(3,28),'second full physical first packet31')
    check(clearance(row2,F(3,31))==F(3,31)>F(1,14),'second literal physical phase')
    print('LOCATED_MULTIPLIER_JOIN_PROFILE_29',counts)
    print('SECOND_ACTUAL_ROW','t',t2,'same_V_U_and_all_one_profiles','packet29_empty','first_packet31',(3,28),'clearance_at3/31',F(3,31))
    check(packet(U,2)==(1,),'opposite-component small packet survives: adaptive selector not refuted')
    check(min(V)<=3*Q//28,'inherited larger-unit closure already applies')
    print('ZERO_OWNER_BY_DENOMINATOR',owners)
    print('FIRST_V_PACKET',29,pv,'FIRST_PHYSICAL_PACKET',29,pr)
    print('INHERITED_JOIN','t mod29',16,'U_PACKET29',pu,'COMPATIBLE',joined)
    print('PHYSICAL_PHASE',F(5,29),'clearance',F(3,29),'lifted_phase',F(166010,332021),'clearance',F(7,29))
    print('SCOPE: actual each-component/full-row denominator<=28 forcing refuted; adaptive choice of U denominator2 NOT refuted.')
    print('SCOPE: larger-unit theorem already proves safety; no new general CRT theorem or LRC closure.')
    print('PASS',N,'always-active exact gates')

if __name__=='__main__':main()
