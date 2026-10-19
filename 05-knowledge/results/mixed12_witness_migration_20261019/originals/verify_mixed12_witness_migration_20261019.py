#!/usr/bin/env python3
"""Exact small-common certificates for mixed one-lower/two-upper roots.
Standard library only. No floating point, optimization, or external imports.
"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from collections import Counter
import json
import hashlib
import argparse

def close(pred):
    pred=list(pred)
    for k in range(len(pred)):
        for j in range(len(pred)):
            if pred[j]>>k&1:pred[j]|=pred[k]
    assert all(not(pred[j]>>j&1) for j in range(len(pred)))
    return tuple(pred)

def extensions(pred):
    n=len(pred);out=[]
    def visit(word,mask):
        if len(word)==n:out.append(tuple(word));return
        for v in range(n):
            if not(mask>>v&1) and not(pred[v]&~mask):visit(word+[v],mask|1<<v)
    visit([],0)
    return out

def cores(c):
    # x=0,y=1,b=2,Z=3..c+2, upper slaves c+3,c+4.
    n=c+5;Z=list(range(3,c+3));A=[c+3,c+4]
    edges=[(2,z) for z in Z]+[(2,a) for a in A]+[(z,a) for z in Z for a in A]+[tuple(A)]+list(combinations(Z,2))
    seen=set()
    for mask in range(1<<len(edges)):
        pred=[0]*n;pred[1]=4
        for a in A:pred[a]=1
        for i,(a,z) in enumerate(edges):
            if mask>>i&1:pred[z]|=1<<a
        pred=close(pred)
        if pred not in seen:seen.add(pred);yield pred

def generators(pred,words,inc):
    gens={};rankwords=[]
    for o in words:
        ranks=[o.index(v) for v in range(len(pred))];rankwords.append((o,ranks))
        h=min(ranks[:2]);k=max(ranks[:2])+1
        counts=[1]+[int(ranks[a]<ranks[b]) for a,b in inc]
        for j in range(h+1):
            for t in range(k,len(pred)+1):
                key=(o[:j],o[t:]);column=gens.setdefault(key,[0]*len(counts))
                for i,v in enumerate(counts):column[i]+=v
    keys=sorted(gens,key=lambda x:(len(x[0]),x[0],len(x[1]),x[1]))
    return keys,[gens[key] for key in keys],rankwords

def witness(rows):
    # A chosen row is 3N_forward-2E. A nonpositive combination is a certificate.
    for size in range(1,4):
        for ids in combinations_with_replacement(range(len(rows)),size):
            if all(sum(rows[i][j] for i in ids)<=0 for j in range(len(rows[0]))):
                values=[sum(rows[i][j] for i in ids) for j in range(len(rows[0]))]
                return Counter(ids),values
    return None,None

def exceptional_witness(pred,M,inc,rows):
    if pred==(0,4,0,4,13,5) and M==(2,0,3,1,5,4):
        pairs={(0,3):3,(1,3):4,(1,5):1,(4,5):1}
    elif pred==(0,4,0,4,5,13) and M==(2,0,3,1,4,5):
        pairs={(0,3):3,(1,3):4,(1,4):1,(4,5):1}
    else:return None,None
    weights={inc.index(pair):w for pair,w in pairs.items()}
    values=[sum(w*rows[i][j] for i,w in weights.items()) for j in range(len(rows[0]))]
    assert max(values)<=0
    return weights,values

def family(c):
    records=[];kinds=Counter();norders=0;total_words=0
    for pred in cores(c):
        words=extensions(pred);n=len(pred)
        inc=[(a,b) for a,b in combinations(range(n),2) if not(pred[b]>>a&1 or pred[a]>>b&1)]
        keys,columns,rankwords=generators(pred,words,inc)
        certificates=[];certificate_hash=hashlib.sha256();core_kinds=Counter();core_orders=0
        total_words+=len(words)
        row_templates=[([3*col[i+1]-2*col[0] for col in columns],[col[0]-3*col[i+1] for col in columns]) for i in range(len(inc))]
        for M,ranks in rankwords:
            span=ranks[1]-ranks[0]
            if span<1 or any(abs(ranks[a]-ranks[b])>span for a,b in inc):continue
            rows=[row_templates[i][int(ranks[a]>ranks[b])] for i,(a,b) in enumerate(inc)]
            weights,values=witness(rows)
            if weights is None:weights,values=exceptional_witness(pred,M,inc,rows)
            assert weights is not None,(c,pred,M)
            oriented=[]
            for i,w in sorted(weights.items()):
                a,b=inc[i]
                if ranks[a]>ranks[b]:a,b=b,a
                oriented.append([a,b,w])
            total=sum(weights.values());kinds[str(total)]+=1;norders+=1
            cert={'M':list(M),'weighted_forward_pairs':oriented,'weight_sum':total,'maximum_generator_margin':max(values),'minimum_generator_margin':min(values)}
            certificate_hash.update((json.dumps(cert,sort_keys=True,separators=(',',':'))+'\n').encode())
            core_kinds[str(total)]+=1;core_orders+=1
            if c<=2:certificates.append(cert)
        record={'predecessor_masks':list(pred),'extension_count':len(words),'generator_count':len(keys),'admissible_orders':core_orders,'weight_sum_counts':dict(sorted(core_kinds.items())),'certificate_stream_sha256':certificate_hash.hexdigest()}
        if c<=2:record['certificates']=certificates
        records.append(record)
    return {'common_count':c,'core_count':len(records),'total_bare_extensions':total_words,'admissible_orders':norders,'weight_sum_counts':dict(sorted(kinds.items())),'cores':records}

def exceptional_table():
    pred=(0,4,0,4,13,5);words=extensions(pred)
    pre=[(),(2,),(2,3)]
    suf=[(),(4,),(5,),(3,4),(4,5),(5,4),(3,4,5),(3,5,4),(5,3,4)]
    table=[]
    for u in pre:
        row=[]
        for v in suf:
            os=[o for o in words if (not u or o[:len(u)]==u) and(not v or o[-len(v):]==v)]
            row.append(sum(6-3*int(o.index(0)<o.index(3))-4*int(o.index(3)<o.index(1))-int(o.index(1)<o.index(5))-int(o.index(5)<o.index(4)) for o in os))
        table.append(row)
    assert table==[[0,1,6,7,9,1,6,3,3],[10,4,8,4,9,3,4,2,2],[7,1,3,0,2,0,0,0,0]]
    assert len(words)==37
    return {'predecessor_masks':list(pred),'extension_count':37,'prefixes':[list(x) for x in pre],'suffixes':[list(x) for x in suf],'deficit_coefficients':table,'bare_comparison_numerators':[sum(o.index(a)<o.index(b) for o in words) for a,b in [(0,3),(3,1),(1,5),(5,4)]]}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--max-common',type=int,choices=(2,3),default=3);args=parser.parse_args()
    families=[family(c) for c in range(args.max_common+1)]
    assert [f['core_count'] for f in families]==[7,40,357,4824][:len(families)]
    assert [f['admissible_orders'] for f in families]==[60,722,15899,609424][:len(families)]
    assert families[0]['weight_sum_counts']=={'1':60}
    assert families[1]['weight_sum_counts']=={'1':710,'2':3,'3':7,'9':2}
    assert families[2]['weight_sum_counts']=={'1':15672,'2':209,'3':18}
    if len(families)==4:assert families[3]['weight_sum_counts']=={'1':601284,'2':8008,'3':132}
    result={'status':'PASS','claim':'No maximum-span majority root with one lower slave, two upper slaves, and at most three common neighbors in a no-balanced poset.','families':families,'exceptional_sharp_inequality':exceptional_table()}
    out=Path(__file__).with_name('mixed12_witness_migration_audit_20261019.json' if args.max_common==3 else 'mixed12_witness_migration_c2_audit_20261019.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','families':[{k:v for k,v in f.items() if k!='cores'} for f in families],'exceptional_inequality':result['exceptional_sharp_inequality']},indent=2))
if __name__=='__main__':main()
