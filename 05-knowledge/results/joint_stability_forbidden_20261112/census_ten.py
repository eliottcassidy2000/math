"""Complete order-ten marked-unit census from McKay's ordinary catalogue.

The large public input is cached outside the repository. Its URL, checksum,
record count, and triangle encoding are retained in the receipt.
"""
from pathlib import Path
import gzip,json,hashlib,time,sys,os
from collections import Counter
import numpy as np
from census_nine import tables,strong,automorphisms,move,need,save
HERE=Path(__file__).resolve().parent
URL='https://users.cecs.anu.edu.au/~bdm/data/tourn10.txt.gz'

def run():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(os.environ['TEMP'])/'math-tourn10-mckay.txt.gz'
    source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
    records=[];hist=Counter();strong_count=unit_sources=unit_masks=unit_orbits=0;start=time.time()
    with gzip.open(path,'rb') as f:
        for idx,line in enumerate(f,1):
            word=line.strip();need(len(word)==45 and set(word)<={48,49},'catalogue upper triangle encoding')
            rows=[0]*10;b=0
            for i in range(10):
                for j in range(i+1,10):
                    if word[b]==49:rows[i]|=1<<j
                    else:rows[j]|=1<<i
                    b+=1
            H,U,W,masks=tables(np.array(rows,dtype=np.int64));H=int(H);hist[H]+=1
            need(H%2==1,'Hamilton parity')
            isstrong=strong(rows);strong_count+=isstrong
            if masks:
                aut=automorphisms(rows);a=len(aut);need(H%a==0 and len(masks)%a==0,'free marked action')
                unit_sources+=1;unit_masks+=len(masks);unit_orbits+=len(masks)//a
                if isstrong:
                    reps=sorted({min(move(int(m),p) for p in aut) for m in masks});need(len(reps)==len(masks)//a,'marked quotient')
                    records.append({'code':int(word,2),'rows':rows,'H':H,'Aut':a,'marks':[{'selected':m,'U':sum(int(U[v]) for v in range(10) if m>>v&1)} for m in reps]})
            if idx%100000==0:print(json.dumps({'classes':idx,'strong_unit_classes':len(records),'seconds':round(time.time()-start,1)}),flush=True)
    need(idx==9733056,'complete ordinary order-ten catalogue')
    result={'scope':'All 9,733,056 ordinary order-ten tournament classes in McKay catalogue, every nonempty selected set, actual selection automorphisms, all strong unit seeds retained. Converse types distinct.','source_url':URL,'source_sha256':source_hash,'ordinary_classes':idx,'strong_classes':strong_count,'all_unit_sources':unit_sources,'all_unit_masks':unit_masks,'all_unit_orbits':unit_orbits,'strong_unit_sources':len(records),'strong_unit_orbits':sum(len(r['marks']) for r in records),'H_histogram':dict(sorted(hist.items())),'records':records}
    save('unit_ten.json.gz',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['records','H_histogram']}),flush=True)
if __name__=='__main__':run()
