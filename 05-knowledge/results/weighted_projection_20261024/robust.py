"""Exact replay of incoming fractional loads and their weighted envelope.

No optimization solver. Incoming native menus and charges are hash-pinned.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import Counter
import json,subprocess,hashlib

ROOT=Path(__file__).resolve().parents[3]
REV='8ca0f5fe460490874c15fbd8758ec3cce1593b97'
FILES={};checks=0
def read(path):
    b=subprocess.check_output(['git','show',REV+':'+path],cwd=ROOT)
    FILES[path]=hashlib.sha256(b).hexdigest();return json.loads(b)
def need(v,msg):
    global checks
    checks+=1
    if not v:raise RuntimeError(msg)
def add(a,b):return (((a>>4)+(b>>4))%128)*16|((a^b)&15)
def incoming_to_here(g):
    x,y=g>>4,g&15
    return (x<<4)|((y<<1)&15)|((x&1)^(y>>3))
def mul(values):
    z=1
    for a in values:z*=a
    return z
def main():
    cert=read('05-knowledge/results/session47_resilient_native_wiring_20261024/helper-resilience-v3/fractional_weights.json')
    menus=read('05-knowledge/results/three_target570_20261024/menus.json')
    E=read('05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json')['order8']
    X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4']
    need(len(X)==3528,'original full carrier')
    for s in range(len(X)):
        old=(X[s]%128)*16+sum((Y[j][s]%2)<<j for j in range(4))
        new=(X[s]%128)*16+sum((Y[j][s]%2)<<(j-1) for j in range(1,5))
        need(Y[0][s]%2==(X[s]%2)^(Y[4][s]%2),'Y0 = X parity xor Y4 at every original node')
        need(incoming_to_here(new)==old,'native parity-basis conversion')
    need(len({incoming_to_here(g) for g in range(2048)})==2048,'charge conversion is bijective')
    for g in range(2048):
        for h in (1,2,4,8,16):need(incoming_to_here(add(g,h))==add(incoming_to_here(g),incoming_to_here(h)),'basis conversion is additive on generators')
    rows={r['source']:r for r in menus['rows']};ids=cert['source_ids'];need(len(ids)==53 and len(set(ids))==53,'53 actual slots')
    labels=[rows[s]['label'][1]*16+rows[s]['label'][0] for s in ids]
    need(all(rows[s]['role']=='helper' for s in ids),'only original helper supply')
    odds=[1+s%3 for s in ids];counts=Counter();weighted=Counter();edgecharge={}
    for I in combinations(range(53),4):
        m=sum(1<<i for i in I);g=0
        for i in I:g=add(g,labels[i])
        counts[g]+=1;weighted[g]+=mul(odds[i] for i in I);edgecharge[m]=g
    lower=[];ratios=[];totaledges=0
    need(len(cert['charges'])==2048 and len(counts)==2048,'all exact charge families')
    for incoming_g,item in enumerate(cert['charges']):
        g=incoming_to_here(incoming_g)
        den=item['denominator'];E=item['edges'];loads=[0]*53;S=0;M=Q(0);vmin=None;seen=set()
        for e,n in E:
            need(e in edgecharge and edgecharge[e]==g and e not in seen,'each certificate edge is a distinct actual charge witness')
            need(isinstance(n,int) and n>0,'positive integer certificate coefficient')
            seen.add(e);I=[i for i in range(53) if e>>i&1];v=mul(odds[i] for i in I)
            for i in I:loads[i]+=n
            S+=n;M=max(M,Q(n,v));vmin=v if vmin is None else min(vmin,v)
        need(max(loads)<=den and S>13*den,'exact fractional transversal certificate')
        ell=max(Q(vmin),Q(S-13*den)/M);lower.append(ell);ratios.append(Q(S,den));totaledges+=len(E)
        need(ell<=weighted[g],'envelope is consistent with unrestricted weighted coefficient')
    L=min(lower);U=max(weighted.values());UU=max(counts.values())
    return {'status':'PASS frozen fractional certificate and source-weighted envelope',
            'checks':checks,'frozen_revision':REV,'inputs':FILES,
            'quadruples':sum(counts.values()),'certificate_edges':totaledges,'minimum_fractional_mass':str(min(ratios)),
            'basis_conversion':'(x,b0,b1,b2,b3) -> (x, x mod2 xor b3, b0,b1,b2); checked on all3528 nodes and G generators',
            'unit_helper_odds':{'multiplicity_min':min(counts.values()),'multiplicity_max':UU,
                'all_13_restrictions_L':1,'all_13_restrictions_U':UU,'core28_gap_bound':str(Q(1,2*UU*378))},
            'explicit_weight_fixture':{'helper_odds':dict(zip(map(str,ids),odds)),'L':str(L),'U':U,
                'core28_gap_bound':str(L/(2*U*378))},
            'scope':'uniform bound for every fixed <=13 forced-U set; no adaptive-stationarity or simultaneous moving-pin claim'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
