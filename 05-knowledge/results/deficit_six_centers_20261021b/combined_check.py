"""Lightweight exact controls for the combined degree-35 reduction.

This pins and consumes the existing complete order-ten receipt; it does not
repeat its census or claim to classify any of the ten remaining orders.
"""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction
from collections import Counter
from itertools import combinations
import argparse
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHECKS = 0
PINS = {
 '05-knowledge/results/odd_shore_deficit_centers_20261021/session22-centers/PROOF.md': '6f2108a23a1268afa3fa55dfbcacd8b1d406cd1fa8156255779fee0f39729d13',
 '05-knowledge/results/odd_shore_deficit_centers_20261021/session22-centers-review/THIRD_REVIEW.md': '9b4f0004fab410315f57264a3133e160292772e9caefa498ffaf9a4beed01df1',
 '05-knowledge/results/paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md': 'f628ee105b839088c40c7144477ee667f8b333cabb9fc3e0a5d11227c53928a7',
 '05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md': '787ab6282c22e9f60ac734d2c3e96c0836d34b8e1d2b0c6fb11c05e5eefd4388',
 '05-knowledge/results/deficiency_three_20261020/centers.md': 'be5717f730c9325b6d68113eb6826b0b614c4e98af0e285eddd56fdaabd8ef3e',
 '05-knowledge/results/metagraph_classification_20261021a/README.md': 'a8fd5875952710217d4398eabcead0868899a335da9d7ff1feabb989487beb62',
 '05-knowledge/results/deficit_six_centers_20261021b/check.py': '52e9b7806c0d2502d81d2c96e02b761e1da44e894a1093d02a89c79a26af1a5e',
 '05-knowledge/results/deficit_six_centers_20261021b/receipt.json': '8ea09c7747176477426a89184df29e5d385fe75e9fdfe11d8ccc70f6586d39af',
}


def ck(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def source(t):
    ck(type(t) in (list,tuple) and all(type(x) is int for x in t), 'exact source rows')
    n = len(t)
    for a,r in enumerate(t):
        ck(0 <= r < 1 << n and not (r >> a & 1), 'source range/diagonal')
    for a,b in combinations(range(n),2):
        ck(((t[a] >> b) & 1)+((t[b] >> a) & 1) == 1, 'tournament pair')


def arc(t,a,b):
    return (t[a] >> b) & 1


def permutation(p,n):
    ck(type(p) in (list,tuple) and all(type(x) is int for x in p)
       and sorted(p) == list(range(n)), 'literal permutation')


def is_map(t,p,opposite=False):
    permutation(p,len(t))
    return all(arc(t,p[a],p[b]) == (arc(t,a,b)^opposite)
               for a,b in combinations(range(len(t)),2))


def path(t,p):
    permutation(p,len(t))
    ck(all(arc(t,a,b) for a,b in zip(p,p[1:])), 'literal Hamilton word')


def target(t,p):
    path(t,p); kept = {frozenset((a,b)) for a,b in zip(p,p[1:])}
    out = [0]*len(t)
    for a,b in combinations(range(len(t)),2):
        i,j = (a,b) if arc(t,a,b) ^ (frozenset((a,b)) not in kept) else (b,a)
        out[i] |= 1 << j
    return out


def load_pins():
    for name,h in PINS.items():
        raw = (ROOT/name).read_bytes().replace(b'\r\n',b'\n')
        ck(sha256(raw).hexdigest() == h, 'inherited proof/receipt pin: '+name)
    return json.loads((HERE/'receipt.json').read_text(encoding='utf-8'))


def consume_ten(old):
    c = old['order_ten_complete_normalized_census']; rows = c['frames']
    ck(old['checks'] == 7041 and c['frame_count'] == len(rows) == 128, 'frozen complete order-ten receipt')
    ck(sorted(x['mask'] for x in rows) == list(range(128)), 'complete exact mask domain')
    hist = Counter()
    for x in rows:
        mask,t = x['mask'],x['source']; source(t)
        ck(type(mask) is int and len(t) == 10, 'exact frame/source size')
        # Reconstruct all 45 arcs independently from the seven frame bits.
        connection = []
        for shore in (0,1):
            connection.append({j if mask >> (2*shore+j-1) & 1 else 5-j for j in (1,2)})
        f = {0} | {j for j in (1,2,3) if mask >> (j+3) & 1}
        for a,b in combinations(range(10),2):
            if a % 2 == b % 2:
                yes = (b//2-a//2) % 5 in connection[a%2]
            elif a % 2 == 0:
                yes = (b//2-a//2) % 5 in f
            else:
                yes = (a//2-b//2) % 5 not in f
            ck(arc(t,a,b) == yes, 'independent normalized frame incidence')
        p = x['hamilton_cycle']; path(t,p)
        ck(arc(t,p[-1],p[0]), 'normalizing cycle closes')
        ck(len(x['automorphisms']) == x['Aut'], 'full inherited automorphism count')
        for a in x['automorphisms']:
            ck(is_map(t,a), 'actual automorphism sidecar')
        for a in x['anti_involutions']:
            ck(is_map(t,a,True) and all(a[a[i]] == i for i in range(10)), 'actual anti-involution sidecar')
        for key in ('H','Aut','h','blue_H','b'):
            ck(type(x[key]) is int and x[key] >= 0, 'integer census count')
        ck(x['Aut'] > 0 and x['H'] == x['h']*x['Aut'] and x['blue_H'] == x['b']*x['Aut'], 'exact free divisions')
        ck(type(x['self_converse']) is bool and x['self_converse'] == bool(x['anti_involutions']), 'inherited SC type')
        ck(target(t,x['black_word']) == x['black_target'], 'source-bound complement target')
        if x['blue_word'] is not None:
            bp = x['blue_word']; path(t,bp)
            r = dict(zip(bp,bp[::-1]))
            ck(is_map(t,[r[i] for i in range(10)],True), 'actual blue word')
        if x['self_converse']:
            hist[x['h']-x['b']] += 1
    expected = {1528:8,1628:8,2908:12,2926:12,3064:12,3114:12}
    ck(hist == expected and sum(hist.values()) == 64 and 6 not in hist, 'complete imported order10 exclusion')
    initial = {x['n']:len(x['strong_delta6']) for x in old['inherited_census']}
    ck(initial == {3:0,4:0,5:1,6:1,7:1,8:0,9:0}, 'pinned inherited small center coefficients')
    return {'frames':128, 'self_converse_frames':64, 'histogram':dict(sorted(hist.items())),
            'scope':'consumption of pinned exact census; map/word sidecars replayed, not a new census'}


def cycle_choice_hostile():
    t = [sum(1 << ((i+j)%7) for j in (1,2,4)) for i in range(7)]
    source(t); first = list(range(7)); closing = [0,1,2,4,5,6,3]
    indices = []
    for p in (first,closing):
        path(t,p); ck(arc(t,p[-1],p[0]), 'actual Paley closing word')
        succ = dict(zip(p,p[1:]+p[:1])); power = list(range(7)); auto = 0
        for k in range(7):
            auto += is_map(t,power)
            power = [succ[x] for x in power]
        indices.append(7//auto)
    ck(indices == [1,7], 'chosen cycle, not every cycle, has q>1')
    r = dict(zip(closing,closing[::-1]))
    ck(not is_map(t,[r[i] for i in range(7)],True), 'nonblue closing sidecar')
    return {'source':t,'cycles':[first,closing],'exact_q':indices}


def support():
    pairs = []
    for q in range(2,9):
        for m in (1,3,5,7):
            n = q*m
            if n < 10 or (q in (6,8) and m >= 5) or (q == 7 and m == 7):
                continue
            ck(m % 2 == 1 and q <= 8, 'literal necessary cut index and odd shore')
            pairs.append([n,q,m])
    pairs.sort()
    before = sorted({n for n,q,m in pairs})
    ck(before == [10,12,14,15,18,20,21,24,25,28,35], 'odd-shore/score-run support')
    after = [n for n in before if n != 10]
    ck(after == [12,14,15,18,20,21,24,25,28,35], 'exact ten remaining orders')
    ck(max(after) == 35 and len(after) == 10, 'combined bound')
    return {'before_order10_exclusion':before, 'remaining_orders':after,
            'remaining_normalized_parameter_pairs':[p for p in pairs if p[0] != 10],
            'coefficient_status':'OPEN at every remaining order; no frame census here'}


def add(a,b):
    c = [0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i] += x
    for i,x in enumerate(b):c[i] += x
    while len(c)>1 and not c[-1]:c.pop()
    return c


def mul(a,b):
    c = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j] += x*y
    return c


def gcd_poly(a,b):
    a,b = list(map(Fraction,a)),list(map(Fraction,b))
    def trim(x):
        while x and not x[-1]:x.pop()
    trim(a);trim(b)
    while b:
        while len(a)>=len(b):
            k=len(a)-len(b);v=a[-1]/b[-1]
            for j,x in enumerate(b):a[k+j]-=v*x
            trim(a)
        a,b=b,a
    return [x/a[-1] for x in a] if a else []


def recurrence():
    p=[1,-1,0,-1];u=[1,0,-1,0,0,0,-1];pp=mul(p,p);uu=mul(u,u)
    den=mul(pp,uu);delta=add(uu,[-x for x in mul(pp,u)])
    ck(all(x%2==0 for x in delta), 'integral half-difference')
    r=[0]*36
    for k in (5,6,7,35):r[k]=1
    num=add([0]*5+[x//2 for x in delta],mul(mul(pp,u),r))
    num=add(num,[0]*10+mul(pp,[1,1,0,1]))
    ck(len(den)-1==18 and len(num)-1==47 and num[-1]==-1, 'new sharp formal degree bill')
    ck(gcd_poly(p,u)==[1], 'coprime denominator factors')
    ck(gcd_poly(p,[i*p[i] for i in range(1,len(p))])==[1], 'squarefree p')
    ck(gcd_poly(u,[i*u[i] for i in range(1,len(u))])==[1], 'squarefree u')
    ck(gcd_poly(u,[1,1,0,1])==[1], 'no center-independent cancellation')
    ck(gcd_poly(num,den)==[1], 'synthetic reduced denominator')
    coefficients=[]
    for n in range(80):
        coefficients.append((num[n] if n<len(num) else 0)-sum(den[j]*coefficients[n-j] for j in range(1,min(n,18)+1)))
    ck(sum(den[j]*coefficients[47-j] for j in range(19))==-1, 'fails at47 for synthetic top atom')
    for n in range(48,80):
        ck(sum(den[j]*coefficients[n-j] for j in range(19))==0, 'homogeneous from48')
    return {'denominator':den,'minimal_eventual_order':18,'safe_onset':48,
            'synthetic_R6_degrees':[5,6,7,35],'synthetic_numerator_degree':47,
            'synthetic_numerator_top':-1,'scope':'formal test, not existence of an order35 center'}


def verify():
    global CHECKS
    CHECKS=0
    old=load_pins();ten=consume_ten(old);choice=cycle_choice_hostile();remaining=support();rec=recurrence()
    return {'status':'PROVED combined bound35; ten residual orders OPEN',
            'incoming_commit':'6b12a62f6b','proof_and_receipt_lf_pins':PINS,
            'order_ten_consumption':ten,'cycle_choice_hostile':choice,
            'remaining_support':remaining,'trivalent_recurrence':rec,'checks':CHECKS}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=verify();data=json.dumps(result,sort_keys=True,indent=2)+'\n';p=HERE/'combined_receipt.json'
    if args.write:p.write_text(data,encoding='utf-8',newline='\n')
    if args.verify:ck(p.read_text(encoding='utf-8')==data,'immutable combined receipt equality')
    print(json.dumps({'checks':result['checks'],'sha256_lf':sha256(data.encode()).hexdigest(),
                      'order_bound':35,'remaining_orders':10,'safe_onset':48},sort_keys=True))


if __name__=='__main__':main()
