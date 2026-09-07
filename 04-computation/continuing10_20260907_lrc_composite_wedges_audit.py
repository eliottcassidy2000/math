#!/usr/bin/env python3
"""Independent referee: boundary partitions, floor quotients, native residues.

No producer module is imported or executed. All checks remain active under -O.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, lcm, comb
from itertools import combinations, combinations_with_replacement, product
from bisect import bisect_left, bisect_right
from collections import Counter
from hashlib import sha256
import argparse
import json
import sys
sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent if HERE.name=='04-computation' else Path('C:/w/s0905')
parser=argparse.ArgumentParser()
parser.add_argument('--profiles', type=Path, default=ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json')
parser.add_argument('--inherited', type=Path, default=ROOT/'05-knowledge/results/continuing9_20260907_lrc_ratio_tree_certificate.json')
parser.add_argument('--certificate', type=Path, default=(ROOT/'05-knowledge/results' if HERE.name=='04-computation' else HERE)/'continuing10_20260907_lrc_composite_wedges_certificate.json')
ARGS=parser.parse_args()
GATES=0
T=7200
def need(ok, label):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(label)
def canonical(x): return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def pinned(path, pin):
    raw=path.read_bytes();need(sha256(raw).hexdigest()==pin, str(path)+' hash')
    return json.loads(raw)
J=pinned(ARGS.certificate,'890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9')
SOURCE=HERE/'continuing10_20260907_lrc_composite_wedges.py'
need(sha256(SOURCE.read_bytes()).hexdigest()=='1d2f52c100802ab1d3b330c920a3f494eca662b7a51a7c5fedc6f1fe9a287cb5','frozen producer source pin; never imported')
OLD=pinned(ARGS.inherited,'f93d2c8b4f112027bbba23fa529742cd12ac8107e7c9bfdf57f15d6666799533')
SUP=pinned(ARGS.profiles,'935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f')

def native(n,p,q,beta,closed=False):
    den=n*beta.denominator; A=beta.numerator; B=beta.denominator
    total=0
    for j in range(n):
        r=p*(B*j+A)%den;s=q*(B*j+A)%den
        lhs=max(14*min(r,den-r),14*min(s,den-s))
        total+=lhs<=den if closed else lhs<den
    return total

def boundary_partition(p,q):
    """Partition at EVERY original endpoint, test exact cell midpoint.

    This uses neither pairwise interval intersection nor two pointers.
    """
    L=14*p*q
    cuts={0,L}
    for speed,scale in [(p,q),(q,p)]:
        for j in range(speed+1):
            for a in [14*j-1,14*j+1]:
                if 0<a*scale<L: cuts.add(a*scale)
    cuts=sorted(cuts); pieces=[]
    for a,b in zip(cuts,cuts[1:]):
        twice=a+b;den=2*L
        r=p*twice%den;s=q*twice%den
        if 14*max(min(r,den-r),min(s,den-s))<den: pieces.append((a,b))
    need(pieces[0][0]==0 and pieces[-1][1]==L,'origin pieces')
    return L,[(pieces[-1][0]-L,pieces[0][1])]+pieces[1:-1]

PROF={}
def audit_profile(rec):
    n,p,q=rec['n'],rec['p'],rec['q'];L,I=boundary_partition(p,q)
    need(L==rec['L'] and len(I)==rec['components'],'native partition size')
    need(sha256(canonical(I)).hexdigest()==rec['intervals_sha256'],'all partition intervals')
    # Independent closed formula: count(beta) = sum(floor(nb/L)-floor(na/L))
    #                  - #{nb mod L <= L beta} + #{na mod L < L beta}.
    # Thus strict wall and chamber values follow from order statistics,
    # without updating a running event count.
    left=sorted(n*a%L for a,b in I);right=sorted(n*b%L for a,b in I)
    base=sum((n*b)//L-(n*a)//L for a,b in I)
    lc=Counter(left);rc=Counter(right);walls=sorted(lc.keys()|rc.keys())
    rows=[];values=[base-bisect_right(right,0)]
    for w in walls:
        at=base-bisect_right(right,w)+bisect_left(left,w)
        after=base-bisect_right(right,w)+bisect_right(left,w)
        rows.append([w,lc[w],rc[w],at,after]);values.extend((at,after))
    need(len(walls)==rec['walls'],'complete wall universe')
    need(sha256(canonical(rows)).hexdigest()==rec['events_sha256'],'every strict wall and chamber value')
    need((min(values),max(values))==(rec['minimum'],rec['maximum']),'whole-circle extrema')
    for label in ['minimum','maximum']:
        beta=F(*rec['minimizer' if label=='minimum' else 'maximizer'])
        need(native(n,p,q,beta)==rec[label],'literal attaining grid '+label)
    PROF[n,p,q]=(rec['minimum'],rec['maximum'])

for rec in J['profiles']: audit_profile(rec)
need(len(PROF)==1154,'all exact profile records')
need(sum(r['walls'] for r in J['profiles'])==2149609,'total complete wall universe')
print('PROFILES 1154; independent original-boundary partition and quotient-order statistics; all 2149609 walls plus chambers and literal extrema PASS',flush=True)

# Reconstruct the strict sum atlas by a smallest-prime-factor sieve, iterating
# sums first instead of iterating the producer's p,q order.
spf=list(range(357))
for p in range(2,19):
    if spf[p]==p:
        for j in range(p*p,357,p):
            if spf[j]==j:spf[j]=p
def good_sum(s):
    while s>1:
        p=spf[s];k=0
        while s%p==0:s//=p;k+=1
        if p%3!=2 or k>2:return False
    return True
ATLAS={(p,s-p) for s in range(3,357) if good_sum(s) for p in range(1,(s+1)//2) if gcd(p,s)==1}
need(ATLAS=={tuple(r) for r in J['atlas']} and len(ATLAS)==5855,'complete sieve atlas')

P={int(k):{(r[0],tuple(r[1])) for r in v['profiles']} for k,v in SUP['levels'].items()}
STATES=sorted(set(SUP['levels']['6']['gcds']) & {d for d in range(1,T+1) if T%d==0})
need(len(STATES)==26 and comb(len(STATES)+3,4)==23751,'four-free multiset universe')
MASKS=range(1,127)
def valid(d):
    # Dynamic subset gcds retain all positional masks, including equal labels.
    if gcd(*d)!=1:return False
    gs=[0]*128
    for mask in MASKS:
        bit=mask&-mask;gs[mask]=gcd(gs[mask-bit],d[bit.bit_length()-1])
        c=gs[mask];comp=tuple(sorted(gcd(c,d[i]) for i in range(7) if not mask>>i&1))
        if (c,comp) not in P[len(comp)]:return False
    return True
def excess(d):return sum(r*(((T-1)//(7*r))+1) for r in d)-T
for block in J['conditional_words']:
    found=[];attempted=0
    for free in combinations_with_replacement(STATES,4):
        attempted+=1;word=tuple(sorted(list(free)+block['forced']))
        if valid(word):found.append([list(word),excess(word)])
    need(attempted==block['attempted']==23751,'full conditional attempts')
    need(sorted(found)==block['words'],'every conditional word and excess')
    need(max(e for d,e in found)==block['maximum'],'attained conditional ceiling')
need(valid(J['dropped_middle_profile_hostile']) and excess(J['dropped_middle_profile_hostile'])==164,'middle-margin hostile')
print('PROPER_PROFILES all126 positional conditions; three complete23751 universes recover1033/108,949/116,162/107; dropped-middle hostile164 PASS',flush=True)

USED=set()
def credits(D,p,q):
    e=gcd(T,D);p,q=sorted((p,q));key=(T//e,p,q);USED.add(key)
    need(key in PROF,'required pair profile present')
    return tuple(e*a for a in PROF[key])
ARMS={}
for block in J['arm_domains']:
    a,b=block['margins'];e=gcd(a,b);n=T//e;found={}
    for p,q in ATLAS:
        for u,v in [(p,q),(q,p)]:
            if [e*gcd(n,u),e*gcd(n,v)]==[a,b]:found[u,v]=credits(e,p,q)[0]
    need(found=={(u,v):c for u,v,c in block['arms']},'complete directed arm set and credits')
    ARMS[a,b]=found
for block in J['wedges']:
    a,mid,b=block['margins'];E=block['ceiling'];first=0;rest={};zero=0
    for (u,v),cu in ARMS[a,mid].items():
        for (w,z),cw in ARMS[b,mid].items():
            if cu+cw>E:first+=1;continue
            # Integer 2x2 path products, no rational propagation/LCM.
            R=[u*z,v*z,w*v];g=gcd(*R);R=[r//g for r in R]
            need(len(set(R))==3 and gcd(*R)==1,'primitive actual triple')
            need([gcd(T,2*r) for r in R]==[a,mid,b],'complete scale-gcd2 witness')
            D=gcd(R[0],R[2]);p,q=sorted((R[0]//D,R[2]//D));e=gcd(T,2*D)
            cr=credits(2*D,p,q)[0];forest=cu+cw+cr-min(cu,cw,cr)
            need(forest>E,'residual actual forest closes')
            rest[u,v,w,z]=(R,[e,p,q,cr],forest);zero+=cu==cw==0
    frozen={tuple(r['arms'][0][:2]+r['arms'][1][:2]):(r['primitive'],r['pair'],r['forest']) for r in block['residual']}
    need(rest==frozen,'every unpruned residual certificate')
    need(first==block['arm_sum_closed'] and first+len(rest)==block['universe']==72996,'complete wedge partition')
    need(min(r[2] for r in rest.values())==block['minimum_residual_credit'] and zero==block['zero_wedges'],'residual sharp and zero-arm controls')
third=J['third_pattern'];R=third['primitive'];found=[]
for i,j in combinations(range(3),2):
    D=gcd(R[i],R[j]);p,q=sorted((R[i]//D,R[j]//D));e=gcd(T,6*D)
    found.append([i,j,e,p,q,credits(6*D,p,q)[0]])
need(found==third['costs'] and [r[-1] for r in found]==[0,96,66],'named actual triple all pair counts')
need(96<=107<162==sum(sorted(r[-1] for r in found)[-2:]),'same-triple forest repairs pair near miss')
print('NATIVE_WEDGES 2*72996; residual57/214 with minima114/142; named third forest162>107 PASS',flush=True)

# Reconstruct each inherited tree by its integer kernel; propagate ratios by
# eliminating leaves and then normalize. The displayed equations are checked
# directly on every recovered primitive row, including rejected rows.
def kernel_tree(edges):
    vals={0:F(1)}
    while len(vals)<7:
        before=len(vals)
        for i,j,num,den in edges:
            if i in vals and j not in vals:vals[j]=vals[i]*num/den
            if j in vals and i not in vals:vals[i]=vals[j]*den/num
        need(len(vals)>before,'connected exact equations')
    L=lcm(*(r.denominator for r in vals.values()));u=[int(vals[i]*L) for i in range(7)];g=gcd(*u);u=[v//g for v in u]
    need(all(u[j]*den==u[i]*num for i,j,num,den in edges),'literal homogeneous equations')
    return u
selected={(r['word_index'],r['orientation_index']):r for r in J['selected_bank']}
total=good=0
for wi,old in enumerate(OLD['fixed_tree_bank']):
    choices=[];d=old['word']
    for i,j,e,p,q in old['tree']:
        local=[];n=T//e
        need((p,q) in ATLAS,'every selected tree edge belongs to strict atlas')
        if [e*gcd(n,p),e*gcd(n,q)]==[d[i],d[j]]:local.append((i,j,q,p))
        if [e*gcd(n,q),e*gcd(n,p)]==[d[i],d[j]]:local.append((i,j,p,q))
        choices.append(local)
    orientations=list(product(*choices));need(len(orientations)==len(old['orientations']),'all inherited orientations')
    for oi,edges in enumerate(orientations):
        total+=1;u=kernel_tree(edges)
        need(u==old['orientations'][oi]['row'],'inherited primitive row')
        ok=len(set(u))==7 and [gcd(T,a) for a in u]==d
        need(ok==old['orientations'][oi]['valid'],'joint clipped valuations and distinctness')
        if not ok:need((wi,oi) not in selected,'rejected realization not promoted');continue
        good+=1;rec=selected[wi,oi];need(u==rec['row'] and valid(d) and excess(d)==rec['E'],'selected actual row/full profiles')
        costs={}
        for r in rec['pairs']:
            i,j=r['i'],r['j'];D=gcd(u[i],u[j]);p,q=sorted((u[i]//D,u[j]//D));e=gcd(T,D)
            need((e,p,q)==(r['e'],r['p'],r['q']) and r['evaluated']==(p+q<=10000),'all21 actual ratios and cutoff')
            if r['evaluated']:
                lo,hi=credits(D,p,q);need((lo,hi)==(r['minimum'],r['maximum']),'selected complete pair profile')
                costs[i,j]=lo
        need(max(costs.values())>rec['E'],'single actual composite pair closes selected row')
        # Verify the recorded forest directly; no claim of complete graph optimum.
        components=[{i} for i in range(7)];weight=0
        for cost,i,j in rec['forest']:
            need(costs[min(i,j),max(i,j)]==cost,'forest actual certified edge')
            A=next(s for s in components if i in s);B=next(s for s in components if j in s)
            need(A is not B,'forest acyclic');components.remove(A);components.remove(B);components.append(A|B);weight+=cost
        need(weight==rec['forest_credit'] and weight>rec['E'],'selected forest credit')
        # Independent Prim maximum-tree calculation, rather than Kruskal.
        seen={0};prim=0
        while len(seen)<7:
            crossing=[(c,j if i in seen else i) for (i,j),c in costs.items() if (i in seen)!=(j in seen)]
            if not crossing:seen.add(next(i for i in range(7) if i not in seen));continue
            c,j=max(crossing);seen.add(j);prim+=c
        need(prim==weight,'maximum forest on evaluated actual edges only')
        graph={frozenset((i,j)) for i,j,*_ in old['tree']}
        for tag,i,m,j in rec['families']:
            if tag=='uniform_wedge':
                need(frozenset((i,m)) in graph and frozenset((m,j)) in graph and d[m]==24 and sorted((d[i],d[j])) in [[4,18],[16,18]],'actual uniform theorem tag')
            else:
                h=gcd(u[i],u[m],u[j]);need([u[i]//h,u[m]//h,u[j]//h]==R and gcd(T,h)==6,'actual named theorem tag')
        need(bool(rec['families']),'selected realization has theorem tag')
need(total==29 and good==20 and len({r['word_index'] for r in selected.values()})==10,'exact selected29/20/10 scope')

need(native(7,1,1,F(1,2))==0 and native(7,1,1,F(1,2),True)==2,'strict versus closed wall hostile')
USED.add((7,1,1))
u=OLD['actual_row'];den=2*T
bad=[[14*min(v*(2*j+1)%den,den-v*(2*j+1)%den)<den for v in u] for j in range(T)]
need(sum(not any(r) for r in bad)==2020,'native zero-tree control safe count')
need(all(sum(row[i] and row[j] for row in bad)==0 for i,j,*_ in OLD['actual_graph']),'all six minima coexist')
for rec in J['full13_controls']:
    A=rec['body'];B=rec['tails'];U=A+B;d=[gcd(T,v) for v in B]
    need(len(U)==len(set(U))==13 and min(U)>1 and gcd(*U)==1 and gcd(*A)==T,'actual13 primitive distinct nonunit and selected six gcd')
    need(d==rec['margins'] and valid(d) and excess(d)==rec['E'],'actual13 every proper profile')
    edges=[]
    for i,j in combinations(range(7),2):
        D=gcd(B[i],B[j]);pq=tuple(sorted((B[i]//D,B[j]//D)))
        if pq in ATLAS:edges.append([i,j])
    need(edges==rec['edges']==[[0,1],[1,2]],'entire strict graph, four isolates')
    den=7*T;safe=[]
    for j in range(T):
        if min(min(v*(7*j+1)%den,den-v*(7*j+1)%den) for v in U)*14>=den:safe.append(j)
    need(len(safe)==rec['safe_count'] and sha256(canonical(safe)).hexdigest()==rec['safe_indices_sha256'],'all actual13 safe lifts')
need(USED==set(PROF),'all1154 profiles used in declared universes, none unaccounted')
print('SELECTED29 orientations,20 valid on10 words; all actual ratios/forests/family tags; omitted large ratios stay unevaluated PASS')
print('LITERAL_CONTROLS full13 safe2541/2504/2476, six simultaneous zero edges with2020 safe, strict0 versus closed2 PASS')
print('PASS',GATES,'always-active independent exact gates; raw LF')
