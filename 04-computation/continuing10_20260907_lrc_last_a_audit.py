#!/usr/bin/env python3
"""Independent complete-atlas and native endpoint referee for the last A word.

Original-boundary partitions and floor-quotient order statistics; no producer
imports. All finite tests raise explicitly and remain active under -O.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
from bisect import bisect_left,bisect_right
from hashlib import sha256
from math import gcd
import argparse
import json
import sys
sys.stdout.reconfigure(encoding='utf-8',newline='\n')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent if HERE.name=='04-computation' else Path('C:/w/s0905')
ap=argparse.ArgumentParser()
ap.add_argument('--primary-dir',type=Path,default=HERE if HERE.name=='04-computation' else Path('C:/w/continuing10_20260907_root_lrc'))
ap.add_argument('--certificate',type=Path)
args=ap.parse_args()
G=0;T=7200;WORD=(1,9,16,18,24,32,60)
def need(v,s):
    global G
    G+=1
    if not v:raise ArithmeticError(s)
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()

# Complete atlas by a sieve over the allowed prime-power factorizations.
primes=[p for p in range(2,357) if all(p%d for d in range(2,int(p**.5)+1))]
allowed={1}
for p in primes:
    if p%3==2:allowed|={n*p**k for n in tuple(allowed) for k in [1,2] if n*p**k<=356}
ATLAS=sorted((p,s-p) for s in allowed if s>=3 for p in range(1,(s+1)//2) if gcd(p,s)==1)
need(len(ATLAS)==5855,'all strict-atlas coprime pairs from permitted prime products')

def components(p,q):
    L=14*p*q;cuts={0,L}
    for v,scale in [(p,q),(q,p)]:
        for j in range(v+1):
            for a in [(14*j-1)*scale,(14*j+1)*scale]:
                if 0<a<L:cuts.add(a)
    cuts=sorted(cuts);pieces=[]
    for a,b in zip(cuts,cuts[1:]):
        r=p*(a+b)%(2*L);s=q*(a+b)%(2*L)
        if 14*max(min(r,2*L-r),min(s,2*L-s))<2*L:pieces.append((a,b))
    need(pieces[0][0]==0 and pieces[-1][1]==L,'origin circular danger component')
    return L,[(pieces[-1][0]-L,pieces[0][1])]+pieces[1:-1]

CACHE={}
def profile(n,p,q):
    p,q=sorted((p,q));key=(n,p,q)
    if key in CACHE:return CACHE[key]
    L,I=components(p,q)
    left=sorted(n*a%L for a,b in I);right=sorted(n*b%L for a,b in I)
    M=sum(n*b//L-n*a//L for a,b in I);walls=sorted(set(left+right))
    lo=hi=M-bisect_right(right,0);arglo=arghi=F(0)
    lc=Counter(left);rc=Counter(right);rows=[]
    for k,w in enumerate(walls):
        nxt=walls[k+1] if k+1<len(walls) else walls[0]+L
        at=M-bisect_right(right,w)+bisect_left(left,w)
        after=M-bisect_right(right,w)+bisect_right(left,w)
        rows.append([w,lc[w],rc[w],at,after])
        for value,phase in [(at,F(w,L)),(after,F(w+nxt,2*L))]:
            if value<lo:lo=value;arglo=phase
            if value>hi:hi=value;arghi=phase
    rec=dict(n=n,p=p,q=q,L=L,components=len(I),walls=len(walls),minimum=lo,maximum=hi,
        minimizer=[arglo.numerator,arglo.denominator],maximizer=[arghi.numerator,arghi.denominator],
        intervals_sha256=sha256(canon(I)).hexdigest(),events_sha256=sha256(canon(rows)).hexdigest())
    CACHE[key]=rec;return rec

def native(n,p,q,phase,closed=False):
    A,B=phase.numerator,phase.denominator;den=B*n;count=0
    for j in range(n):
        r=p*(B*j+A)%den;s=q*(B*j+A)%den
        cost=14*max(min(r,den-r),min(s,den-s))
        count+=cost<=den if closed else cost<den
    return count

raw=(ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
need(sha256(raw).hexdigest()=='580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d','inherited full minimum-tree certificate')
OLD=next(r for r in json.loads(raw)['clocks'] if r['t']==T)
oldweights={tuple(d):v[0] for d,v in OLD['weights']}
TABLE={};CANDIDATES={}
for a,b in combinations(WORD,2):
    e=gcd(a,b);n=T//e;found=[]
    for p,q in ATLAS:
        margins=sorted((e*gcd(n,p),e*gcd(n,q)))
        if margins==[a,b]:found.append((e*profile(n,p,q)['minimum'],p,q))
    need(bool(found),'locally compatible atlas edge exists')
    cost,p,q=min(found);TABLE[a,b]=cost;CANDIDATES[a,b]=found
    need(cost==oldweights[a,b],'independent entire-atlas minimum for each of21 pairs')
    rec=profile(n,p,q);need(e*native(n,p,q,F(*rec['minimizer']))==cost,'native first attaining minimum of each margin pair')
need(TABLE=={(1,9):114,(1,16):135,(1,18):133,(1,24):126,(1,32):128,(1,60):136,
 (9,16):142,(9,18):0,(9,24):126,(9,32):142,(9,60):135,
 (16,18):142,(16,24):0,(16,32):0,(16,60):136,
 (18,24):0,(18,32):142,(18,60):102,(24,32):0,(24,60):0,(32,60):144},'complete literal21-edge table')
need(all(TABLE[1,b]>116 for b in WORD[2:]),'margin1 forces only9-neighbor')
need(min(c for e,c in TABLE.items() if c and e!=(1,9))==102,'any second positive minimum with114 closes')
zero={e for e,c in TABLE.items() if c==0}
need(zero=={(9,18),(16,24),(16,32),(18,24),(24,32),(24,60)},'complete remaining zero graph')
remaining=zero|{(1,9)}
# Enumerate every subgraph of this seven-edge generous graph, and impose
# connectedness plus the already-proved forbidden16--24--18 wedge.
survivors=[]
def connected(edges):
    seen={1}
    for _ in range(6):seen|={b for a,b in edges if a in seen}|{a for a,b in edges if b in seen}
    return len(seen)==7
for bits in product([0,1],repeat=len(remaining)):
    edges={e for bit,e in zip(bits,sorted(remaining)) if bit}
    if connected(edges) and not ({(16,24),(18,24)}<=edges):survivors.append(edges)
tree={(1,9),(9,18),(18,24),(24,32),(16,32),(24,60)}
need(survivors==[tree],'unique entire remaining actual connected graph')
need(gcd(32,24)==8 and gcd(18,24)==6 and 114+6>116 and 114+8>116,'actual nonzero arm multiplicity exceeds remaining budget2')
print('TABLE all21 weights independently minimized over complete5855 atlas; unique residual connected six-edge tree PASS',flush=True)

ARMS={};ALL={};NEW=set()
for a,expected_size,expected_low in [(32,317,14),(18,231,5)]:
    e=gcd(a,24);n=T//e;full=[]
    for p,q in ATLAS:
        for u,v in [(p,q),(q,p)]:
            if [e*gcd(n,u),e*gcd(n,v)]==[a,24]:
                rec=profile(n,p,q);full.append([u,v,e*rec['minimum']]);NEW.add((n,p,q))
    full.sort();low=[r for r in full if r[2]<=2]
    need(len(full)==expected_size and len(low)==expected_low and all(r[2]==0 for r in low),'complete directed arm alphabet and low-credit set')
    ALL[a]=full;ARMS[a]=low
found=[]
for (u,v,c1),(w,z,c2) in product(ARMS[32],ARMS[18]):
    R=[u*z,v*z,w*v];g=gcd(*R);R=[r//g for r in R]
    need(len(set(R))==3 and gcd(*R)==1,'actual distinct primitive triple')
    need([gcd(T,2*r) for r in R]==[32,24,18],'actual scale-gcd2 margins')
    d=gcd(R[0],R[2]);p,q=sorted((R[0]//d,R[2]//d));e=gcd(T,2*d);n=T//e
    rec=profile(n,p,q);NEW.add((n,p,q));credit=e*rec['minimum']
    need(e==2 and credit==142 and credit>116,'every actual composite uniform minimum142')
    found.append(dict(arms=[[u,v,c1],[w,z,c2]],primitive=R,valid=True,e=e,p=p,q=q,credit=credit))
need(len(found)==70,'complete14 times5 product universe')
for key in sorted(NEW):
    rec=CACHE[key]
    for label,phase in [('minimum','minimizer'),('maximum','maximizer')]:
        need(native(*key,F(*rec[phase]))==rec[label],'new profile literal attaining grid '+label)
need(native(7,1,1,F(1,2))==0 and native(7,1,1,F(1,2),True)==2,'strict versus closed endpoint hostile')
print('ARMS full317/231; low14/5 allzero; every70 actual primitive endpoint minimum142>116 PASS',flush=True)
print('PROFILE_UNIVERSE',len(CACHE),'complete-atlas pair profiles;',len(NEW),'native arm/endpoint profiles with literal extrema')
source=args.primary_dir/'continuing10_20260907_lrc_last_a.py'
need(sha256(source.read_bytes()).hexdigest()=='cc2668c8249bd1d1c804c65af93a3e55a2663a6c4891bdf797983771c98eb325','frozen primary source; never imported')
cp=args.certificate or ((ROOT/'05-knowledge/results') if HERE.name=='04-computation' else args.primary_dir)/'continuing10_20260907_lrc_last_a_certificate.json'
raw=cp.read_bytes();need(sha256(raw).hexdigest()=='bce97b849036ba89ca8e8d1593b2ab03ff7aa9ecc647bdc72d55bd4962383274','frozen primary certificate')
J=json.loads(raw)
need(tuple(J['word'])==WORD and J['E']==sum(d*((T+7*d-1)//(7*d)) for d in WORD)-T==116,'exact physical marginal excess')
need(J['pair_table']==[[i,j,TABLE[WORD[i],WORD[j]]] for i,j in combinations(range(7),2)],'entire primary21-edge table')
index={a:i for i,a in enumerate(WORD)}
need(J['cheap_edges']==sorted([index[a],index[b]] for a,b in remaining),'complete primary affordable graph')
need(J['connected_graphs']==[sorted([index[a],index[b]] for a,b in tree)],'unique primary entire connected graph')
for block in J['arms']:
    a,b=block['margins'];need(b==24 and sorted(block['rows'])==ALL[a] and sorted(block['low'])==ARMS[a],'all primary directed arms')
key=lambda r:tuple(r['arms'][0][:2]+r['arms'][1][:2])
ind={key(r):dict(arms=r['arms'],primitive=r['primitive'],scale_gcd=2,e=r['e'],p=r['p'],q=r['q'],credit=r['credit']) for r in found}
need(ind=={key(r):r for r in J['products']},'every primary primitive path product')
need(len(J['profiles'])==618 and {(r['n'],r['p'],r['q']) for r in J['profiles']}==NEW,'complete primary native profile universe')
for rec in J['profiles']:need(CACHE[rec['n'],rec['p'],rec['q']]==rec,'all primary moment-free native interval/event hashes and extrema')
control=J['full13_control'];A=control['body'];B=control['tails'];U=A+B
need(A==[T*i for i in range(1,7)] and len(U)==len(set(U))==13 and min(U)>1 and gcd(*U)==1,'actual distinct primitive nonunit13 control')
need(tuple(gcd(T,b) for b in B)==WORD,'actual control ordered clock margins')
edges=[]
for i,j in combinations(range(7),2):
    D=gcd(B[i],B[j]);pq=tuple(sorted((B[i]//D,B[j]//D)))
    if pq in ATLAS:edges.append([i,j])
seen={0}
for _ in range(6):seen|={j for i,j in edges if i in seen}|{i for i,j in edges if j in seen}
need(edges==control['edges'] and len(seen)==7,'entire actual strict graph connected')
den=7*T;safe=[j for j in range(T) if 14*min(min(u*(7*j+1)%den,den-u*(7*j+1)%den) for u in U)>=den]
need(len(safe)==control['safe_count']==2579 and sha256(canon(safe)).hexdigest()==control['safe_hash'],'all actual thirteen-row safe lifts')
print('PRIMARY all618 native records,548 directed arms,70 products and21 entries match; literal full13 safe2579 PASS')
print('PASS',G,'always-active independent exact gates; raw LF')
