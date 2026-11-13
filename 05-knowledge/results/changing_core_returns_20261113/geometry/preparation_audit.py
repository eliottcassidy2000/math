"""Independent recount of all supplied small-patch endpoints; no scan rerun."""
from pathlib import Path
from math import gcd
import hashlib,json
from neutral_scc import marked,flip,sccs,relabel


def need(test,message):
    if not test:raise RuntimeError(message)


def main():
    here=Path(__file__).resolve().parent;prep=here.parent/'preparation'
    sources=json.loads((here.parent/'selector/receipt.json').read_text())['sources']
    patch=json.loads((prep/'neutral_881_patch.json').read_text())
    three=json.loads((prep/'neutral_881_three.json').read_text())
    four=json.loads((prep/'neutral_881_four.json').read_text())
    checked=[]
    for source in sources:
        for rowkey in ('rows','prepared_rows'):
            z=marked(source[rowkey],source['root'])
            need(z[0]==881,'selector cube endpoint H')
            checked.append(z)
    for i,source in enumerate(patch['sources']):
        for record in source['global_probe']['two_arc_H881_endpoints']:
            rows=tuple(sources[i]['rows'])
            for a,b in record['arcs']:rows=flip(rows,a,b)
            need(marked(rows,sources[i]['root'])==tuple(record['state']),'two-arc endpoint')
            checked.append(tuple(record['state']))
    for record in three['H881_endpoints']:
        rows=tuple(sources[1]['rows'])
        for a,b in record['arcs']:rows=flip(rows,a,b)
        need(marked(rows,10)==tuple(record['state']),'three-arc endpoint')
        checked.append(tuple(record['state']))
    strong=[];good=[]
    for record in four['validated_H881_endpoints']:
        rows=tuple(sources[1]['rows'])
        for a,b in record['arcs']:rows=flip(rows,a,b)
        z=marked(rows,10)
        need(z==tuple(record['state']),'four-arc endpoint')
        isstrong=len(sccs(rows))==1
        need(isstrong==record['strong'],'four-arc strongness')
        if isstrong:strong.append(z)
        h,s,q,F=z;D=F+s-h-q;d=F-q
        hits=[]
        for k in range(1,401):
            t=d*2**k+q-D;b=3*t+2*D
            admissible=gcd(b,t)==1 and t%3==0 and b%3!=0
            wheel=record['wheel']
            need(admissible==(k%wheel['period'] in wheel['good']),'endpoint wheel control')
            if admissible:hits.append(k)
        if hits:good.append(z)
        checked.append(z)
    need(strong==[(881,97,97,941)],'unique filtered strong four-arc endpoint')
    need(len(good)==4,'four successful weak endpoints')
    chosen=four['selected_nine_vertex_core'];rows=tuple(sources[1]['rows']);names=list(range(11))
    for a,b in chosen['arcs']:rows=flip(rows,a,b)
    z=marked(rows,10);deletions=[]
    need(sccs(rows)==[tuple(range(8))+(10,),(8,),(9,)],'selected strong core followed by tail')
    for old in (9,8):
        remove=names.index(old);keep=[i for i in range(len(names)) if i!=remove]
        rows=relabel(rows,keep);names=[names[i] for i in keep];root=names.index(10)
        need(marked(rows,root)==z,'individual tail deletion preserves full state')
        deletions.append(dict(deleted_old_name=old,state=z,order=len(rows)))
    need(rows==(354,228,425,179,231,192,140,1,250),'nine-vertex final masks')
    need(z==(881,169,201,977) and len(sccs(rows))==1,'strong nine-vertex profile')
    ret=chosen['core_repair']['literal_return']
    need(marked(ret['rows'],ret['root'])==tuple(ret['state'])==(5195,1689,1689,5195),'literal nine-core return')
    need(all(h+F-s!=1691 for h,s,q,F in checked),'next1691 absent from supplied patch endpoint bank')
    result={'status':'FINITE-EXACT','optimized':not __debug__,
            'scope':'Independent endpoint recount; exhaustive patch enumeration completeness was reviewed in source code, not rerun.',
            'patch_endpoints_recounted':len(checked),'strong_filtered_width4':strong,
            'successful_weak_width4':len(good),'selected_deletion_steps':deletions,
            'selected_return':ret['state'],
            'reviewed_source_hashes':{name:hashlib.sha256((prep/name).read_bytes()).hexdigest()
                for name in ['NEUTRAL_881.md','neutral_881_four.cpp','neutral_881_four.py']}}
    (here/('preparation_audit.optimized.json' if not __debug__ else 'preparation_audit.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':'PASS','patch_endpoints_recounted':len(checked),'selected_core':z,'selected_return':ret['state']}))


if __name__=='__main__':main()
