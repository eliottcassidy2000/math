#!/usr/bin/env python3
"""Standard-library exact replay of the frozen, scoped mixed-(1,2), c=4 packet."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json
ROOT=Path(__file__).resolve().parent
FROZEN=ROOT
spec=importlib.util.spec_from_file_location('catalogue',ROOT/'verify_saturated_near_twin_lemma_20261019.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)
CHECK=B.check

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def pairs_mask(pairs,inc):
    mask=bits=0
    for a,b,w in pairs:
        CHECK(type(w) is int and w>0,'Certificate weight is not a positive integer')
        CHECK(a!=b and tuple(sorted((a,b))) in inc,'Certificate pair is not incomparable')
        i=inc.index(tuple(sorted((a,b))));bit=1<<i
        CHECK(not (mask&bit),'Repeated comparison in normalized certificate')
        mask|=bit
        if a<b:bits|=bit
    CHECK(mask!=0,'Empty certificate')
    return mask,bits

def generator_data(pred,words,pairs_lists):
    succ=B.successors(pred);sat={v for v in range(9) if 8-(pred[v]|succ[v]).bit_count()==7}
    totals=[sum(w for a,b,w in pairs) for pairs in pairs_lists]
    columns={};rankwords=[]
    for word in words:
        ranks=[0]*9
        for i,v in enumerate(word):ranks[v]=i
        rankwords.append((word,ranks))
        h=min(ranks[0],ranks[1]);k=max(ranks[0],ranks[1])+1
        jmax=min(h,min((ranks[v] for v in sat),default=9))
        tmin=max(k,max((ranks[v]+1 for v in sat),default=0))
        values=[1]+[3*sum(w for a,b,w in pairs if ranks[a]<ranks[b])-2*total
                    for pairs,total in zip(pairs_lists,totals)]
        for j in range(jmax+1):
            for t in range(tmin,10):
                col=columns.setdefault((word[:j],word[t:]),[0]*len(values))
                for i,value in enumerate(values):col[i]+=value
    keys=sorted(columns,key=lambda key:(len(key[0]),key[0],len(key[1]),key[1]))
    return keys,[columns[key] for key in keys],rankwords

def replay_core(record,pred,known):
    words=B.extensions(pred)
    inc=list((a,b) for a,b in combinations(range(9),2) if not(pred[a]>>b&1 or pred[b]>>a&1))
    certs=[dict(c) for c in record.get('distinct_certificates',[])]
    masks=[pairs_mask(c['weighted_forward_pairs'],inc) for c in certs]
    for c,(mask,bits) in zip(certs,masks):
        if 'mask' in c:CHECK(c['mask']==mask and c['bits']==bits,'Producer orientation mask mismatch')
    keys,columns,rankwords=generator_data(pred,words,[c['weighted_forward_pairs'] for c in certs])
    if 'core_extensions' in record:CHECK(record['core_extensions']==len(words),'Extension count mismatch')
    if 'retained_generators' in record:CHECK(record['retained_generators']==len(keys),'Generator count mismatch')
    for i,c in enumerate(certs):
        margins=[col[i+1] for col in columns]
        CHECK(max(margins)<=0,'Certificate has positive exact generator margin')
        CHECK(c['maximum_generator_margin']==max(margins) and c['minimum_generator_margin']==min(margins),
              'Producer exact margin receipt mismatch')
    known_orders={tuple(o) for o in known.get('orders',[])}
    count=supplements=0;stream=hashlib.sha256()
    for word,ranks in rankwords:
        span=ranks[1]-ranks[0]
        if span<=0 or any(abs(ranks[a]-ranks[b])>span for a,b in inc):continue
        count+=1;orientation=sum(1<<i for i,(a,b) in enumerate(inc) if ranks[a]<ranks[b])
        hit=next((i for i,(mask,bits) in enumerate(masks) if orientation&mask==bits),None)
        CHECK(hit is not None,'Uncovered admissible majority order')
        stream.update((json.dumps([word,hit],separators=(',',':'))+'\n').encode())
    CHECK(count==record['admissible_orders'],'Admissible-order count mismatch')
    normalized=[]
    for c in certs:
        normalized.append({key:c[key] for key in ['weighted_forward_pairs','maximum_generator_margin','minimum_generator_margin']})
    return {'core_index':record['core_index'],'predecessor_masks':list(pred),'core_extensions':len(words),
            'retained_generators':len(keys),'admissible_orders':count,'distinct_certificates':normalized,
            'prior_case_supplementary_certificates':supplements,
            'replay_case_certificate_stream_sha256':stream.hexdigest()}

def load_catalogue():
    coverage=json.loads(gzip.decompress((FROZEN/'c4_combined_analytic_coverage_20261019.json.gz').read_bytes()))
    kept=[];h=hashlib.sha256()
    for p in B.c4_cores():
        succ=B.successors(p);degree=[8-(p[v]|succ[v]).bit_count() for v in range(9)]
        if max(degree)<=7:kept.append(p);h.update((json.dumps(p,separators=(',',':'))+'\n').encode())
    CHECK(len(kept)==79121 and h.hexdigest()==coverage['degree_valid_catalogue_sha256'],'Catalogue mismatch')
    covered=set()
    for record in coverage['certificates']:
        i=record['core_index'];CHECK(type(i) is int and 0<=i<79121,'Invalid analytic core index');CHECK(i not in covered,'Repeated analytic core');covered.add(i)
        p=kept[i];CHECK(tuple(record['predecessor_masks'])==p,'Analytic core mismatch')
        U=B.successors(p);degree=[8-(p[v]|U[v]).bit_count() for v in range(9)]
        u,v=record['witness']['pair'];CHECK(not(p[v]>>u&1 or p[u]>>v&1),'Analytic pair is comparable')
        delta=lambda a,b:(p[b]&~p[a]).bit_count()+(U[a]&~U[b]).bit_count()
        if record['scope']=='AMBIENT_BALANCED_PAIR':
            e=lambda a:0 if a in (0,1) else 7-degree[a]
            CHECK(max(delta(u,v),delta(v,u))+e(u)+e(v)<=1,'Analytic balanced-pair inequality failed')
        else:
            CHECK(record['scope']=='MAXIMUM_SPAN_MAJORITY_CONTRADICTION','Unknown analytic scope')
            z=record['witness']['common_vertex'];CHECK(3<=z<7 and (u,v) in ((0,z),(z,1)),'Bad forced orientation')
            CHECK(delta(u,v)+7-degree[z]<=1,'Analytic directed inequality failed')
    residual=[i for i in range(79121) if i not in covered]
    CHECK(len(covered)==56630 and len(residual)==22491 and residual==coverage['residual_core_indices'],'Analytic partition mismatch')
    return kept,residual

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int,default=None,
        help='Optional prefix control only; omitted for the complete exact replay.')
    args=parser.parse_args()
    manifest_path=ROOT/'manifest.json'
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text())
        for name,item in manifest['files'].items():
            CHECK(digest(ROOT/name)==item['sha256'],'Manifest mismatch: '+name)
    kept,residual=load_catalogue()
    path=ROOT/'normalized_residual_certificates_20261019.jsonl.gz'
    records=[json.loads(line) for line in gzip.decompress(path.read_bytes()).splitlines()]
    if args.limit is None:CHECK(len(records)==22491,'Incomplete residual certificate packet')
    n=len(records) if args.limit is None else min(args.limit,len(records))
    CHECK(n>=0,'Negative control limit')
    counts=Counter();stream=hashlib.sha256()
    for pos,record in enumerate(records[:n]):
        index=residual[pos]
        CHECK(record['core_index']==index and record['residual_position']==pos,'Residual record order mismatch')
        result=replay_core(record,kept[index],{})
        CHECK(result['replay_case_certificate_stream_sha256']==record['replay_case_certificate_stream_sha256'],
              'Case/certificate stream mismatch')
        counts['residual_cores']+=1;counts['admissible_orders']+=result['admissible_orders']
        counts['weighted_certificates']+=len(result['distinct_certificates'])
        counts['core_extensions']+=result['core_extensions'];counts['retained_generators']+=result['retained_generators']
        stream.update((json.dumps([index,result['replay_case_certificate_stream_sha256']],separators=(',',':'))+'\n').encode())
        if (pos+1)%1000==0:print(json.dumps({'checked_residual_cores':pos+1}),flush=True)
    output={'status':'PASS_FULL_SCOPED_C4_REPLAY' if args.limit is None else 'PASS_PREFIX_CONTROL',
            'scope':'Seven-thin mixed-(1,2), c=4 maximum-span majority-root exclusion only.',
            'analytic_cores':56630,**dict(counts),'combined_covered_cores':56630+n,
            'normalized_certificates_sha256':digest(path),
            'analytic_coverage_sha256':digest(ROOT/'c4_combined_analytic_coverage_20261019.json.gz'),
            'aggregate_case_stream_sha256':stream.hexdigest(),'verifier_sha256':digest(Path(__file__))}
    if args.limit is None:
        (ROOT/'replayed_c4_controls_20261019.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
