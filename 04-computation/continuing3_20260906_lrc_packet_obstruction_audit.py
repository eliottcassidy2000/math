"""Independent standard-library referee for both packet rows and the bounded scale class.

No producer imports. Atlas sums are generated multiplicatively; weighted ranks
use rational elimination; safe packets are complements of unioned unsafe masks.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd,prod,lcm
from pathlib import Path
import argparse,hashlib,json,sys
sys.stdout.reconfigure(newline='\n')
N=0


def need(ok,reason):
    global N
    N+=1
    if not ok:raise ArithmeticError(reason)


def atlas_sums():
    primes=[]
    for n in range(2,357):
        if all(n%p for p in primes if p*p<=n):primes.append(n)
    sums={1}
    for p in primes:
        if p%3!=2:continue
        old=set(sums)
        sums.update(a*p for a in old if a*p<=356)
        sums.update(a*p*p for a in old if a*p*p<=356)
    return sums


def graph(row,sums):
    es=[];parents=list(range(len(row)))
    def root(i):
        while parents[i]!=i:i=parents[i]
        return i
    for i,j in combinations(range(len(row)),2):
        d=gcd(row[i],row[j]);height=(row[i]+row[j])//d
        if height in sums:
            es.append((i,j));parents[root(i)]=root(j)
    groups={}
    for i in range(len(row)):groups.setdefault(root(i),[]).append(i)
    return tuple(es),tuple(sorted(tuple(g) for g in groups.values()))


def weighted_rank(row,es):
    matrix=[]
    for i,j in es:
        d=gcd(row[i],row[j]);v=[F(0)]*len(row);v[i]=F(row[j]//d);v[j]=-F(row[i]//d)
        need(sum(v[k]*row[k] for k in range(len(row)))==0,'literal weighted edge relation')
        matrix.append(v)
    rank=0
    for col in range(len(row)):
        pivot=next((j for j in range(rank,len(matrix)) if matrix[j][col]),None)
        if pivot is None:continue
        matrix[pivot],matrix[rank]=matrix[rank],matrix[pivot]
        a=matrix[rank][col];matrix[rank]=[v/a for v in matrix[rank]]
        for j in range(rank+1,len(matrix)):
            a=matrix[j][col]
            if a:matrix[j]=[v-a*w for v,w in zip(matrix[j],matrix[rank])]
        rank+=1
    return rank


def packet(row,d):
    unsafe=0
    for v in row:
        mask=0
        for k in range(d):
            r=(v*k)%d
            if 14*min(r,d-r)<d:mask|=1<<k
        unsafe|=mask
    return tuple(k for k in range(d) if not (unsafe>>k)&1)


def distance(row,x):
    q=x.denominator;p=x.numerator
    return F(min(min((p*v)%q,q-(p*v)%q) for v in row),q)


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path,default=here.parent if here.name=='04-computation' else Path('C:/w/s0905'));args=parser.parse_args()
    name='continuing3_20260906_lrc_packet_obstruction'
    pp=here/(name+'.py');po=here/(name+'.out')
    if not po.exists():po=here.parent/'05-knowledge'/'results'/(name+'.out')
    need(hashlib.sha256(pp.read_bytes()).hexdigest()=='66eacc9647b765c392fb080ed41a031fd50a307b6d58b49629d2cdf342832080','frozen primary source pin')
    need(hashlib.sha256(po.read_bytes()).hexdigest()=='b3e340d2c5731d9e9b0bf48b78ff8f8b84db887cf1b5c375e129a364a63ee0f4','frozen primary output pin')
    Q=91**6;qs=(16,25,27,161,187,247);P=lcm(*range(1,29));V=tuple(sorted(P//q for q in qs));U=(1,43,67,79,97,109,127)
    need(P==prod(qs)==80313433200,'independent lcm/cofactor normalization')
    need(all(gcd(a,b)==1 for a,b in combinations(qs,2)) and gcd(*V)==1,'disjoint prime groups and primitive cofactor shape')
    need(all(gcd(v,u)==1 for v in V for u in U),'fixed shape cross-coprimality')
    sums=atlas_sums();ve,vc=graph(V,sums);ue,uc=graph(U,sums)
    need(ve==((0,3),(0,5),(1,3),(1,4),(2,3),(4,5)),'complete independent V edge set')
    need(ue==((0,1),(0,2),(0,5),(1,2),(1,6),(2,4),(3,5),(5,6)),'complete independent U edge set')
    need(vc==(tuple(range(6)),) and uc==(tuple(range(7)),),'both primitive components connected')
    need(tuple((V[i]+V[j])//gcd(V[i],V[j]) for i,j in ve)==(274,263,214,212,188,41),'published V reduced sums')
    need(tuple((U[i]+U[j])//gcd(U[i],U[j]) for i,j in ue)==(44,68,110,110,170,164,188,236),'published U reduced sums')
    dmin=min(gcd(v,w) for v,w in combinations(V,2));maxsum=max((u+w)//gcd(u,w) for u,w in combinations(U,2))
    lower1=Q//dmin+1;lower2=Q*maxsum//min(V)+1;tmin=max(lower1,lower2);tmax=(Q*Q-sum(U))//sum(V)
    need((dmin,maxsum,tmin,tmax)==(1738800,236,412164,25880486141010),'entire class exact thresholds')
    need(tmin*dmin>Q and tmin*min(V)>Q*maxsum,'both inequalities hold uniformly from class lower endpoint')
    need((tmin-1)*min(V)<=Q*maxsum,'lower endpoint sharp for stated amplitude sufficient gate')
    need(tmax*sum(V)+sum(U)<=Q*Q<(tmax+1)*sum(V)+sum(U),'integer upper endpoint of physical box')
    need(tmin*min(V)>355*max(U),'all cross ratios uniformly outside atlas')
    need(prod(U)==305613336829 and gcd(tmin,prod(U))==1,'explicit class coprimality modulus and nonvacuity')
    raw=(args.repo/'05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json').read_bytes().replace(b'\r\n',b'\n')
    need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','canonical LF inherited profile pin')
    data=json.loads(raw)['levels']
    profiles={int(d):{(c,tuple(gs)) for c,gs in rec['profiles']} for d,rec in data.items()}
    for d in range(1,7):need((1,(1,)*d) in profiles[d],'each exact inherited all-one profile present')
    pv=packet(V,29);pu=packet(U,29)
    need(pv==(5,6,7,9,10,11,12,17,18,19,20,22,23,24),'complete six-shape packet29')
    need(pu==(5,9,10,12,14,15,17,19,20,24),'complete seven-shape packet29')
    need(packet(U,2)==(1,),'adaptive opposite-component packet survives')
    for d in range(2,29):
        need(P%d==0 and any(gcd(q,d)==1 for q in qs),'distributed owner proof at every claimed denominator')
        need(any(v%d==0 for v in V) and packet(V,d)==(),'all-numerator six-shape exclusions')
    join=tuple(len(set(pu)&{k for k in range(29) if (m*k)%29 in pv}) for m in range(29))
    need(join==(0,8,8,2,8,2,0,8,6,4,2,6,2,4,10,10,4,2,6,2,4,6,8,0,2,8,2,8,8),'entire located multiplier join profile')
    summary=[]
    for t,den,expected,clear in [(412164,29,(5,12,17,24),F(3,29)),(412183,31,(3,28),F(3,31))]:
        row=tuple(t*v for v in V)+U
        need(tmin<=t<=tmax and gcd(t,prod(U))==1,'literal row belongs to proved bounded class')
        need(len(set(row))==13 and gcd(*row)==1 and sum(row)<Q*Q,'primitive distinct physical row within box')
        es,groups=graph(row,sums)
        need(groups==(tuple(range(6)),tuple(range(6,13))) and len(es)==14,'literal entire physical graph')
        need(weighted_rank(row,es)==11,'direct rational weighted relation rank')
        need(all(gcd(t*v,u)==1 for v in V for u in U),'literal scaled cross gcds')
        count1=count2=0
        for v,w in combinations(V,2):
            D=t*gcd(v,w)
            for u in U:
                need(D//gcd(D,u)>Q,'first mixed orientation cannot have a nonzero outside coefficient')
                count1+=1
        for u,w in combinations(U,2):
            D=gcd(u,w)
            for v in V:
                need(F(t*v,gcd(D,t*v))>Q*((u+w)//D),'second mixed orientation exceeds every signed box sum')
                count2+=1
        need((count1,count2)==(105,126),'both complete crossing orientations retained')
        bodycount=0
        for k in range(7,13):
            for body in combinations(range(13),k):
                c=gcd(*(row[i] for i in body));outside=[i for i in range(13) if i not in body]
                actual=(c,tuple(sorted(gcd(c,row[i]) for i in outside)))
                need(c==1 and actual in profiles[13-k],'actual full body profile reconstructed')
                bodycount+=1
        need(bodycount==4095,'entire inherited large-body universe for this row')
        for d in range(2,den):need(packet(row,d)==(),'every numerator excluded below claimed first physical denominator')
        need(packet(row,den)==expected,'complete first physical packet')
        phase=F(expected[0],den)
        need(distance(row,phase)==clear>F(1,14),'literal physical safe phase and exact clearance')
        literal29=packet(row,29)
        need(literal29==tuple(k for k in pu if (t*k)%29 in pv),'direct mask packet equals component join')
        if t==412164:
            need(sum(row)==5135637177376615,'published first physical sum')
            need(t*dmin==716670763200 and F(t*min(V),Q*maxsum)==F(4786336882800,4786326552917),'published strict crossing margins')
            need(distance(row,F(166010,332021))==F(7,29),'separate lifted phase control')
        else:need(literal29==() and join[t%29]==0 and any(v%30==0 for v in row),'zero29 alignment and explicit30 owner')
        summary.append((t,den,expected,str(clear)))
    need(min(V)<=3*Q//28,'inherited larger-unit closure scope')
    print('BOUNDED_CLASS',tmin,tmax,'gcd(t,305613336829)=1; analytic monotonic margins and full physical box')
    print('LITERAL_ROWS',summary)
    print('EXACT_ENTRY both weighted ranks11, actual6+7 graph, all231 crossing supports excluded in each row')
    print('PROFILES all4095 body/complement profiles reconstructed per row; canonical inherited JSON checked')
    print('PACKETS six-shape first29; full rows first29 and31; all29 scale multipliers checked; U denominator2 survives')
    print('PASS',N,'always-active gates; no primary imports, multiplicative atlas and rational-rank/unsafe-mask referee')


if __name__=='__main__':main()
