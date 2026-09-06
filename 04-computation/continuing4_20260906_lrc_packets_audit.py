"""Independent exact audit: no producer imports or source-derived constants.

Atlas product sieve; fraction-free weighted rank; complement-mask profiles;
unsafe-residue bitsets; literal center-pair overlap enumeration. All gates
remain active under -O. Only the pinned inherited profile JSON is input data.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod, isqrt
from functools import reduce
from pathlib import Path
import argparse
import hashlib
import json
import sys
sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    if not ok:raise RuntimeError(label)
    GATES+=1


def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def atlas_sums():
    out={1}
    for p in range(2,357):
        if prime(p) and p%3==2:
            old=tuple(out)
            for v in old:
                for power in (p,p*p):
                    if v*power<=356:out.add(v*power)
    return out-{1}


def graph(v,allowed):
    return [(i,j) for i,j in combinations(range(len(v)),2)
            if (v[i]+v[j])//gcd(v[i],v[j]) in allowed]


def components(n,edges):
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:i=parent[i]
        return i
    for i,j in edges:parent[find(j)]=find(i)
    return sorted(sorted(i for i in range(n) if find(i)==r)
                  for r in {find(j) for j in range(n)})


def integer_rank(a):
    a=[row[:] for row in a];r=0
    for col in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][col]),None)
        if k is None:continue
        a[r],a[k]=a[k],a[r]
        for k in range(r+1,len(a)):
            if a[k][col]:
                pivot,factor=a[r][col],a[k][col]
                a[k]=[pivot*x-factor*y for x,y in zip(a[k],a[r])]
                d=reduce(gcd,map(abs,a[k]))
                if d:a[k]=[x//d for x in a[k]]
        r+=1
        if r==len(a):break
    return r


def packet(row,d):
    allowed=(1<<d)-1
    for v in row:
        bad=0;r=0
        for j in range(d):
            if 14*min(r,d-r)<d:bad|=1<<j
            r=(r+v)%d
        allowed &= ~bad
    return tuple(j for j in range(d) if (allowed>>j)&1)


def clearance(row,a,b):
    return F(min(min(v*a%b,(-v*a)%b) for v in row),b)


def center_pair_lengths(p,q):
    p,q=sorted((p,q));den=p*q;out=[]
    for i in range(p):
        for j in range(q):
            k=(q*i-p*j)%den
            distance=min(k,den-k)
            if 14*distance<p+q:
                out.append(F(min(2*p,p+q-14*distance),14*den))
    return sorted(out)


def ceil(q):return -(-q.numerator//q.denominator)


def main():
    ap=argparse.ArgumentParser()
    default_profiles=Path(__file__).resolve().parent/'overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    if not default_profiles.exists():
        default_profiles=Path('C:/w/s0905/04-computation')/default_profiles.name
    ap.add_argument('--profiles',type=Path,default=default_profiles)
    args=ap.parse_args()
    here=Path(__file__).resolve().parent
    source=here/'continuing4_20260906_lrc_packets.py'
    transcript=here/'continuing4_20260906_lrc_packets.out'
    if not transcript.exists():transcript=here.parent/'05-knowledge'/'results'/transcript.name
    need(hashlib.sha256(source.read_bytes()).hexdigest()=='22db0faeddf433aed8d1daeb664976813e4ed8a0b5a0fec57fdeaac027a23b08','frozen producer source pin')
    need(hashlib.sha256(transcript.read_bytes()).hexdigest()=='46c71d72aefab9d3e120c58317cc616c8ad37c1db3b7de9e2335c887a8fae151','frozen producer output pin')
    raw=args.profiles.read_bytes()
    need(hashlib.sha256(raw).hexdigest()=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f','inherited profile bytes')
    profiles={int(k):{(int(d),tuple(w)) for d,w in v['profiles']} for k,v in json.loads(raw)['levels'].items()}
    Q=91**6
    leaves=(179,191,251,257,263);groups=(81,101,125,128,161,209,221)
    P=prod(leaves);H=prod(groups)
    V=tuple(sorted([P]+[83*(P//a) for a in leaves]))
    U=tuple(sorted(H//q for q in groups))
    need(Q==567869252041 and P==580028043449 and H==973398810384000,'literal product reconstruction')
    for xs in (leaves,groups):
        for a,b in combinations(xs,2):need(gcd(a,b)==1,'cofactor coprimality')
    need(reduce(gcd,V)==reduce(gcd,U)==1 and gcd(prod(V),H)==1,'primitive shapes and disjoint supports')
    allowed=atlas_sums()
    ratio_count=sum(gcd(a,n-a)==1 for n in allowed for a in range(1,(n+1)//2))
    need(ratio_count==5855,'independent complete inert quotient atlas size')
    ve,ue=graph(V,allowed),graph(U,allowed)
    need(components(6,ve)==[list(range(6))] and len(ve)==5,'literal six-star topology')
    need(components(7,ue)==[list(range(7))] and len(ue)==8,'literal seven-body topology')
    dmin=min(gcd(a,b) for a,b in combinations(V,2))
    threshold=max(Q//dmin+1,430*Q//min(V)+1)
    need(dmin==712259437 and threshold==1334,'all-parameter crossing dominance threshold')
    need(1334*dmin>Q and 1334*min(V)>430*Q,'strict lower endpoint')
    need(1333*min(V)<=430*Q,'preceding integer fails this selected amplitude gate')
    need(97096*sum(V)+sum(U)<Q*Q,'entire bounded class in physical box')
    need(all(v%2 for v in V) and {u%2 for u in U}=={0,1},'half-phase and mixed parity')
    need(min(V)>F(3*Q,28) and min(U)>1 and min(V)>1,'no unit or small-minimum shortcut asserted')
    for r in range(1,7):need((1,(1,)*r) in profiles[r],'full all-one word supplier membership')
    for d in range(2,29):
        need(any(u%d==0 for u in U),'literal divisor obstruction')
        need(packet(U,d)==(),'complete U numerator exclusion')
    a29,b29=packet(V,29),packet(U,29)
    need(a29==(3,7,8,9,11,13,14,15,16,18,20,21,22,26),'V29 exact packet')
    need(b29==(1,2,4,5,24,25,27,28),'U29 exact packet')
    need(not packet(U,30) and packet(U,31),'first two useful U clocks')
    joins=tuple(len(set(b29)&{j for j in range(29) if t*j%29 in a29}) for t in range(29))
    need(joins==(0,0,2,4,6,2,0,4,8,8,6,8,0,6,2,2,6,0,8,6,8,8,4,0,2,6,4,2,0),'full scale-residue join profile')
    lengths={pair:center_pair_lengths(*pair) for pair in combinations(groups,2)}
    selected=lengths[(81,101)]
    need(len(selected)==25 and sum(selected)==F(389,19089),'literal selected center intersections')
    need(selected==sorted([F(1,707)]+[F(min(162,182-14*k),14*81*101) for k in range(1,13) for _ in range(2)]),'full individual length multiset')
    # Correct central width is 1/(7*101)=1/707.
    cmin=sum(ceil(1334*l)-1 for l in selected)
    need(cmin==13 and cmin>6,'uniform unpooled class closure')
    need(center_pair_lengths(1,2)==[F(1,14)],'contained-interval hostile control')
    need(sum(14*min(7*j%7,(-7*j)%7)<7 for j in range(7))>1,'noncoprime sheet defeats the coprime marginal count')
    results=[]
    for t,first,aa,cl,nsafe in [(1369,31,4,F(4,31),462),(1373,29,2,F(3,29),458),(1583,37,1,F(5,37),530)]:
        row=tuple(t*v for v in V)+U
        need(1334<=t<=97096 and gcd(t,H)==1,'scale class membership')
        need(reduce(gcd,row)==1 and len(set(row))==13 and sum(row)<Q*Q,'actual primitive distinct boxed row')
        ee=graph(row,allowed)
        need(components(13,ee)==[list(range(6)),list(range(6,13))],'complete actual decoder graph')
        mat=[]
        for i,j in ee:
            d=gcd(row[i],row[j]);r=[0]*13;r[i]=row[j]//d;r[j]=-row[i]//d
            need(sum(r[k]*row[k] for k in range(13))==0 and max(map(abs,r))<=355,'literal bounded weighted row')
            mat.append(r)
        need(integer_rank(mat)==11,'independent integer weighted rank')
        crossing=[0,0]
        for ids in combinations(range(13),3):
            small=[i for i in ids if i<6];large=[i for i in ids if i>=6]
            if len(small)==2 and len(large)==1:
                i,j=small;k=large[0];d=gcd(row[i],row[j])
                need(d//gcd(d,row[k])>Q,'full first crossing orientation')
                crossing[0]+=1
            if len(small)==1 and len(large)==2:
                k=small[0];i,j=large;d=gcd(row[i],row[j]);delta=gcd(d,row[k])
                need(row[k]//delta>Q*((row[i]+row[j])//d),'full reverse orientation, every outside multiple')
                crossing[1]+=1
        need(crossing==[105,126],'all 231 mixed supports')
        total=0
        for r in range(1,7):
            for omitted in combinations(range(13),r):
                mask=sum(1<<i for i in omitted)
                d=reduce(gcd,(row[i] for i in range(13) if not (mask>>i)&1))
                word=tuple(sorted(gcd(d,row[i]) for i in omitted))
                need(d==1 and (d,word) in profiles[r],'complete complement profile')
                total+=1
        need(total==4095,'proper large-body profile universe')
        for d in range(2,first):need(not packet(row,d),'every preceding physical denominator')
        bank=packet(row,first)
        need(bank and aa==bank[0] and clearance(row,aa,first)==cl,'first physical packet and exact phase')
        for d in (29,31):
            pv,pu=packet(V,d),packet(U,d)
            need(pv and pu and gcd(t,d)==1,'both component banks and invertible scale')
            need(packet(row,d)==tuple(j for j in pu if t*j%d in pv),'located physical CRT identity')
        E=7*((t+6)//7)-t
        unpooled=sum(ceil(t*l)-1 for l in selected)
        need(unpooled>E and 0<=E<=6,'interval count versus full marginal excess')
        need(gcd(t,7)==1 and all(gcd(t,u)==1 for u in U),'individual sheet and strict-endpoint coprimality')
        pooled={p:t*sum(ls)-len(ls) for p,ls in lengths.items()}
        positive={p:v for p,v in pooled.items() if v>0}
        rounded={p:ceil(t*sum(ls))-len(ls) for p,ls in lengths.items()}
        if t in (1369,1373):
            need(set(positive)=={(81,101)} and sum(positive.values())<E,'all 21 pooled pair credits fail even summed')
            need(sum(max(0,v) for v in rounded.values())<=E,'pooled ceiling still fails')
            need(unpooled==15,'unpooled control credit')
        safe=[]
        den=2*t
        for j in range(t):
            num=2*j+1
            margin=min(min(u*num%den,(-u*num)%den) for u in U)
            if 14*margin>=den:safe.append((j,F(margin,den)))
        need(len(safe)==nsafe and len(safe)>=unpooled-E,'independent integer half-phase grid count')
        for j,c in safe:
            need(c>F(1,14) and clearance(row,2*j+1,den)==c,'literal full-row strict grid phase')
        results.append({'t':t,'sum':sum(row),'edges':ee,'crossings':crossing,'profiles':total,'first':first,
                        'packet':bank,'phase':str(F(aa,first)),'clearance':str(cl),'E':E,'unpooled':unpooled,
                        'positive_pooled':{str(p):str(v) for p,v in positive.items()},'safe_lifts':len(safe)})
    print('PASS proof audit: actual entries, full profiles, packet failure, and inherited whole-class grid closure.')
    print('Atlas ratios',ratio_count,'entry threshold',threshold,'uniform separate interval credit',cmin)
    print('Primitive V',V,'U',U)
    print('Compatibility profile modulo29',joins)
    for r in results:print(json.dumps(r,sort_keys=True))
    print('SCOPE: fixed-clock packet obstruction; no new unresolved row or LRC closure claim.')
    print('PASS',GATES,'always-active exact gates; no producer imports')


if __name__=='__main__':main()
