#!/usr/bin/env python3
"""Complete strong gap1-ready order11 fibre over H1429, and sharp order12 control."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import surgery_orders as grammar
import surgery_ready as ready

GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    need(grammar.U[10]==987,'independent one-start class-series count')
    hist=Counter();values=[];rows=[]
    for index in range(987):
        T=grammar.unrank_one(10,index)
        observed=ready.counts(T,(0,))
        dual=ready.counts(ready.converse(T))
        need(ready.strong(T) and dual['ends'][0]==1,'actual strong unique-start parent')
        H=observed['H']+observed['F'][0]-1
        hist[H]+=1;values.append(H);rows.append(T)
    need(sum(hist.values())==987 and hist[1429]==0,'complete newborn order11 scalar fibre empty')
    controls=[]
    chosen={0,493,986,values.index(min(values)),values.index(max(values))}
    chosen.add(min(range(987),key=lambda i:abs(values[i]-1429)))
    for index in sorted(chosen):
        target=ready.converse(ready.birth(rows[index],1))
        H=ready.counts(target)['H'];S=ready.counts(ready.converse(target),(10,))['F'][10]
        need(ready.strong(target) and H==values[index] and S-H==1,'independent literal newborn control')
        controls.append({'rank':index,'rows':target,'H':H,'aux':10})
    word=(0,0,1,1,1,1,1,0,1)
    positive,r,a=ready.ready_target((0,),(0,),word)
    obs=ready.counts(positive,(r,));dual=ready.counts(ready.converse(positive),(a,))
    need(len(positive)==12 and ready.strong(positive) and obs['H']==1429
         and dual['F'][a]==1430,'sharp actual order12 positive control')
    total=GATES+grammar.GATES+ready.GATES
    rawhist=(json.dumps(sorted(hist.items()),separators=(',',':'))+'\n').encode()
    result={'scope':'all987 strong gap1-ready rooted order11 classes, no H1429; actual order12 positive',
            'classes':987,'target':1429,'hits':0,'minimum':min(values),'maximum':max(values),
            'distinct_counts':len(hist),'histogram_sha256':sha256(rawhist).hexdigest(),
            'literal_controls':controls,'positive_order12':{'rows':positive,'root':r,'aux':a,
                'H':obs['H'],'F':obs['F'][r],'S':dual['F'][a]},'gates':total}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();args.output.write_bytes(raw)
    print(json.dumps({'gates':total,'classes':987,'hits':0,'receipt_bytes':len(raw),
                      'receipt_sha256':sha256(raw).hexdigest()},sort_keys=True))


if __name__=='__main__':main()
