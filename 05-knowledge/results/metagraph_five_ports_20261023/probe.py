"""Rooted full-five census over all survivors of the complete v4 menu."""
from pathlib import Path
from itertools import combinations, permutations
from collections import Counter, defaultdict
import json, hashlib, argparse

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=0
def need(x,why):
    global checks
    checks+=1
    if not x: raise RuntimeError(why)

def matrix(n,code):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(code>>(n*(n-1)//2-1-k))&1;a[j][i]=1-a[i][j]
    return a

def word(a,order):
    z=0
    for i,j in combinations(order,2):z=2*z+a[i][j]
    return z

def root_lut():
    can=[]
    for m in range(1024):
        a=matrix(5,m)
        can.append([min(word(a,(v,)+p) for p in permutations(w for w in range(5) if w!=v)) for v in range(5)])
    types=sorted({x for row in can for x in row});index={c:i for i,c in enumerate(types)}
    return [[index[c] for c in row] for row in can],types

def features(a,lut,types):
    n=len(a);rooted=[[0]*len(types) for _ in range(n)];v4=[0]*n;roles=[[0]*5 for _ in range(n)]
    for S in combinations(range(n),4):
        if sorted(sum(a[u][v] for v in S) for u in S)==[1,1,2,2]:
            for v in S:v4[v]+=1
    for S in combinations(range(n),5):
        m=word(a,S);ds=[sum(a[u][v] for v in S) for u in S]
        for i,v in enumerate(S):rooted[v][lut[m][i]]+=1
        if sorted(ds)==list(range(5)):
            for i,v in enumerate(S):roles[v][4-ds[i]]+=1
    return rooted,v4,roles

def records(a,initial,r):
    rec=[tuple(x) for x in initial]
    for _ in range(r):rec=[(rec[v],tuple(sorted(rec[w] for w in range(len(a)) if a[v][w]))) for v in range(len(a))]
    return tuple(sorted(rec))

def histogram(xs):return [[list(k) if isinstance(k,tuple) else k,v] for k,v in sorted(Counter(xs).items())]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    lut,types=root_lut();need(len(types)==48,'rooted five orbit count');projection=[]
    dual=[]
    for c in types:
        a=matrix(5,c);b=[[a[j][i] for j in range(5)] for i in range(5)]
        dc=min(word(b,(0,)+p) for p in permutations(range(1,5)));dual.append(types.index(dc))
    for m in range(1024):
        for root in range(5):need(lut[m^1023][root]==dual[lut[m][root]],'converse rooted orbit transport')
    for c in types:
        a=matrix(5,c);projection.append(sum(sorted(sum(a[u][v] for v in S) for u in S)==[1,1,2,2] for S in combinations(range(5),4) if 0 in S))
    out={'rooted_five_types':len(types),'rooted_type_representatives':types,'orders':{},'drt23':[]}
    for n in (8,9):
        menu=json.loads((HERE/f'menu_n{n}.json').read_text());base=next(x for x in menu['menus'] if x['name']=='v4' and x['round']==1);need(base['all_collisions_recorded'],'collision list truncated')
        allcodes=sorted({c for g in base['collisions'] for c in g});features_by={};rows=[]
        for code in allcodes:
            a=matrix(n,code);r,v,t=features(a,lut,types);features_by[code]=(a,r,v,t)
            need(all(sum(x*y for x,y in zip(rr,projection))==(n-4)*vv for rr,vv in zip(r,v)),'rooted projection')
        counts=[]
        for rounds in (0,1,2):
            total=base['signatures']-base['collision_groups'];bad=[]
            for group in base['collisions']:
                bank=defaultdict(list)
                for code in group:
                    a,r,_,_=features_by[code];bank[records(a,r,rounds)].append(code)
                total+=len(bank);bad.extend(sorted(g) for g in bank.values() if len(g)>1)
            counts.append({'round':rounds,'augmented_signatures':total,'collisions':sorted(bad)})
        for code,(a,r,v,t) in features_by.items():
            rows.append({'code':code,'v4':v,'tt5_roles':t,'rooted5':r})
        out['orders'][str(n)]={'tested_classes':len(allcodes),'base_one_round_signatures':base['signatures'],'repairs':counts,'records':rows}
        print(n,'full rooted5 repairs',counts,flush=True)
    # Exact cardinality minimum within the declared 48-coordinate menu.
    pairs=[];record_bank={}
    for n in (8,9):
        records_n=out['orders'][str(n)]['records'];record_bank.update({(n,r['code']):r['rooted5'] for r in records_n})
        menu=json.loads((HERE/f'menu_n{n}.json').read_text());base=next(x for x in menu['menus'] if x['name']=='v4' and x['round']==1)
        pairs.extend(((n,a),(n,b)) for group in base['collisions'] for a,b in combinations(group,2))
    def separates(cols,a,b):
        return sorted(tuple(row[k] for k in cols) for row in record_bank[a])!=sorted(tuple(row[k] for k in cols) for row in record_bank[b])
    single_failures=[sum(not separates((k,),a,b) for a,b in pairs) for k in range(len(types))]
    winning_pairs=[list(cols) for cols in combinations(range(len(types)),2) if all(separates(cols,a,b) for a,b in pairs)]
    need(min(single_failures)>0 and len(winning_pairs)==14,'coordinate minimum changed')
    need(winning_pairs[0]==[6,24] and [types[i] for i in winning_pairs[0]]==[66,204],'motif addresses')
    out['coordinate_minimum']={'scope':'joint vertex histograms appended to one-round v4; complete n8/n9 inherited classes',
        'collision_pairs':len(pairs),'single_coordinate_missed_pairs':single_failures,'successful_coordinate_pairs':winning_pairs,
        'selected_indices':winning_pairs[0],'selected_rooted_words':[types[i] for i in winning_pairs[0]],'minimum_cardinality':2}
    need(all(set(cols)!={dual[k] for k in cols} for cols in winning_pairs),'converse-closed pair unexpectedly succeeds')
    need([types[dual[k]] for k in (6,24,35)]==[457,204,66],'selected converse types')
    for records_n in record_bank.values():
        for row in records_n:
            even=row[6]+row[35];odd=row[6]-row[35]
            need(even>=abs(odd) and (even-odd)%2==0 and ((even+odd)//2,(even-odd)//2)==(row[6],row[35]),'integer parity recovery')
    out['converse_closure']={'rooted_type_dual_indices':dual,'minimum_separating_converse_closed_cardinality':3,
        'selected_indices':[6,24,35],'selected_rooted_words':[66,204,457],
        'coordinate_rule':'E=c66+c457, O=c66-c457, S=c204; converse(E,O,S)=(E,-O,S); E>=abs(O), E=O modulo2',
        'scope':'all-size coordinate identities; minimum only in the declared 48 raw-coordinate finite n8/n9 menu'}
    src=ROOT/'04-computation/experiments/sierpinski_ocf_20261022/drt23_classes.json';data=json.loads(src.read_text())
    out['drt23_source_sha256']=hashlib.sha256(src.read_bytes()).hexdigest()
    for d in data[1:3]:
        a=matrix(23,int(d['code'],2));r,v,t=features(a,lut,types);need(sum(x[0] for x in t)==d['TT5'],'DRT TT5 agreement')
        out['drt23'].append({'label':d['label'],'inherited_H':d['H'],'inherited_census5':d['census5'],'rooted_tt5_role_histogram':histogram(map(tuple,t)),'source_histogram':histogram(x[0] for x in t),'rooted5_histogram':histogram(map(tuple,r))})
        print(d['label'],'source5',histogram(x[0] for x in t),flush=True)
    need(out['drt23'][0]['inherited_census5']==out['drt23'][1]['inherited_census5'],'wrong DRT hostile')
    out['checks']=checks
    if args.write:(HERE/'probe.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    else:need(json.loads((HERE/'probe.json').read_text())==out,'frozen probe receipt')

if __name__=='__main__':main()
