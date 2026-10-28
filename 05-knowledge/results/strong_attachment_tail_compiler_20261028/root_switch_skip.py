"""Construct each root-switch block once; exact parent transport is inherited."""
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import argparse,contextlib,gzip,io,json
from native import oracle,need,ROOT

HERE=Path(__file__).resolve().parent

def generate(parents,N,*,skip_blocks=True):
    start=perf_counter();n=N-1
    auts=[oracle.isomorphisms(a,a,True) for a in parents]
    canonmasks=[tuple(min(oracle.move_mask(s,g) for g in group) for s in range(1<<n)) for group in auts]
    buckets=defaultdict(list)
    for k,a in enumerate(parents):buckets[oracle.colours(a)[1]].append(k)
    @lru_cache(None)
    def classify(a):
        for k in buckets[oracle.colours(a)[1]]:
            maps=oracle.isomorphisms(a,parents[k])
            if maps:
                p=maps[0]
                need(all(((a[i]>>j)&1)==((parents[k][p[i]]>>p[j])&1) for i in range(n) for j in range(n)),
                     'parent transport adjacency')
                return k,p
        raise RuntimeError('deletion parent missing')
    states=[(k,s) for k,masks in enumerate(canonmasks) for s in sorted(set(masks))]
    stateid={s:i for i,s in enumerate(states)}
    assigned=bytearray(len(states));children=[];aut_hist=Counter();blocks=Counter();trace=0;calls=0;computed=0
    for idx,(k,s) in enumerate(states):
        if assigned[idx] and skip_blocks:continue
        new_block=not assigned[idx];computed+=1
        child=oracle.extend(parents[k],s);row=Counter()
        for v in range(N):
            deletion,mask=oracle.delete(child,v);parent,p=classify(deletion);calls+=1
            key=parent,canonmasks[parent][oracle.move_mask(mask,p)]
            row[stateid[key]]+=1
        need(sum(row.values())==N and idx in row,'row complete and own root included')
        need(all(w%2 for w in row.values()),'odd vertex orbits')
        if not new_block:
            need(all(assigned[j] for j in row),'full-row control remains in assigned block')
            continue
        # A complete row lists every root orbit of this unrooted child.
        need(not any(assigned[j] for j in row),'different discovered blocks are disjoint')
        for j in row:assigned[j]=1
        children.append(child);trace+=sum(row.values());blocks[len(row)]+=1
        aut_hist[row[idx]*sum(oracle.move_mask(s,g)==s for g in auts[k])]+=1
    need(all(assigned),'all pointed states covered')
    need(trace==N*len(children),'root-switch trace')
    need(len(children)==oracle.burnside(N),'independent post-generation Burnside count')
    result=dict(order=N,pointed_states=len(states),classes=len(children),computed_rows=computed,
        omitted_identical_rows=len(states)-computed,parent_classification_calls=calls,
        parent_cache_size=classify.cache_info().currsize,trace=trace,
        automorphism_order_histogram=dict(sorted(aut_hist.items())),block_size_histogram=dict(sorted(blocks.items())),
        seconds=perf_counter()-start)
    return children,result

def main():
    p=argparse.ArgumentParser();p.add_argument('--max-order',type=int,default=9)
    p.add_argument('--output',required=True);p.add_argument('--compare-small',action='store_true');args=p.parse_args()
    need(1<=args.max_order<=9,'explicit complete universe through nine')
    parents=[()];results=[];small=[]
    for n in range(1,args.max_order+1):
        before=parents;parents,result=generate(before,n);results.append(result)
        print(json.dumps({k:v for k,v in result.items() if 'histogram' not in k}),flush=True)
        if args.compare_small and n<=8:
            with contextlib.redirect_stdout(io.StringIO()):reference,old=oracle.generate(before,n)
            need(reference==parents,'identical ordered child representatives')
            need(old['automorphism_order_histogram']==result['automorphism_order_histogram'],'same automorphism histogram')
            small.append(dict(order=n,identical_representatives=True,inherited_seconds_including_extra_child_audits=old['seconds']))
    comparison=None
    if args.max_order==9:
        path=ROOT/'05-knowledge/results/relation_tower_frobenius_ports_20261028/tournament_bank.json.gz'
        reference=json.loads(gzip.decompress(path.read_bytes()))
        need(reference['order']==9 and reference['rows']==[list(a) for a in parents],'entire saved order-nine bank')
        encoded=json.dumps({'order':9,'rows':parents},separators=(',',':'))+'\n'
        comparison=dict(all_ordered_representatives_identical=True,classes=len(parents),
                        canonical_json_sha256=sha256(encoded.encode()).hexdigest(),
                        source_bank='05-knowledge/results/relation_tower_frobenius_ports_20261028/tournament_bank.json.gz')
    receipt=dict(status='FINITE-EXACT',results=results,small_comparison=small,order_nine_bank=comparison,
                 note='Generation uses parent attachments and root switches; reference counts/bank are post-generation checks.')
    Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
