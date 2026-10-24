"""Complete native alphabet upper bound for an externally pinned pair section."""
from array import array
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,combinations_with_replacement
from pathlib import Path
import csv,gzip,json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
U,V,W=1965,2234,570;N=0
def need(ok,why):
    global N
    N+=1
    if not ok:raise RuntimeError(why)
def code(a):return(a[0]<<4)|sum(b<<j for j,b in enumerate(a[1:]))
def plus(a,b):return((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def main():
    nodespath=RESULTS/'metagraph_order8_20261023/merged_nodes.tsv';ownerpath=RESULTS/'metagraph_order8_20261023/owner_map.u16.gz'
    electricpath=RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    supplypath=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
    phasepath=RESULTS/'guarded_component_completion_20261024/sections/phase_sections_receipt.json'
    nodes=list(csv.DictReader(nodespath.open(),delimiter='\t'))
    merged={int(z[k]):int(z['merged_id'])for z in nodes for k in('class_a','class_b')}
    raw=array('H');raw.frombytes(gzip.decompress(ownerpath.read_bytes()))
    if sys.byteorder!='little':raw.byteswap()
    owner=array('H',(merged[x]for x in raw));full=len(owner)-1
    need(len(owner)==1<<21 and len(nodes)==3528,'entire order8 native normalized-word carrier')
    e=json.loads(electricpath.read_text())['order8'];X=e['higher_layers'][-1]['lifts'][0];Y=e['radical_lifts_mod4']
    labels={s:((X[s]%128)<<4)|sum((Y[j][s]%2)<<(j-1)for j in range(1,5))for s in range(len(nodes))}
    need([s for s in labels if labels[s]==labels[V]]==[W,U,V],'complete paired target-code fibre')
    support=defaultdict(int);wordcounts=Counter();bit={U:1,V:2,W:4}
    for m,s in enumerate(owner):
        t=owner[full^m]
        if t in bit:support[s]|=bit[t];wordcounts[t]+=1
    uv=sorted(s for s,b in support.items()if b&3==3)
    triple=sorted(s for s,b in support.items()if b&2 and b&5)
    allV=sorted(s for s,b in support.items()if b&2)
    need(len(uv)==100 and len(triple)==120 and len(allV)==381,'exact source supports including all loops and converse branches')
    need(dict(Counter(support.values()))=={1:283,2:261,4:165,3:95,5:23,6:20,7:5},'complete three-target support pattern histogram')
    supply=json.loads(supplypath.read_text());oldpath=RESULTS/supply['old_domain'];old=json.loads(oldpath.read_text())['sources']
    oldids=sorted(z['source']for z in old);newids=[z['source']for z in supply['all_new_source_labels']]
    need(uv==sorted(oldids+newids),'earlier72+28 alphabet is independently complete')
    for r in supply['helpers']:need(labels[r['source']]==code(r['label']),'original helper charge reconstructed')
    all54=sorted([z['source']for z in supply['helpers']]+[2125]);free51=[s for s in all54 if s not in(368,776,3139)]
    fixed=plus(labels[776],labels[3139]);need(fixed==1290,'two pinned V leaves exact fixed charge')
    def image(ids,repeats):
        pairs=combinations_with_replacement(ids,2)if repeats else combinations(ids,2)
        return{plus(labels[a],labels[b])for a,b in pairs}
    rows=[]
    for name,ids,expected in[
        ('existing_free51',free51,(972,991)),
        ('all_unused_free69',[s for s in newids if s not in(368,776,3139)],(1404,1415)),
        ('all100_UV',uv,(1853,1854)),
        ('all120_paired_target_switchable',triple,(1988,1988)),
        ('all381_V_incident_but_not_all_switchable',allV,(2048,2048))]:
        distinct=image(ids,False);repeated=image(ids,True)
        need((len(distinct),len(repeated))==expected,'complete pair-image sizes '+name+' got '+str((len(distinct),len(repeated))))
        missing=sorted(set(range(2048))-repeated)
        rows.append(dict(name=name,sources=ids,source_count=len(ids),distinct_pair_image=len(distinct),repeated_pair_image=len(repeated),
                         missing_even_with_repeats=missing,missing_total_helper_charges=[plus(fixed,q)for q in missing]))
    need(37 not in image(triple,True) and plus(fixed,37)==1327,'literal all120 hostile: residual37, total1327')
    # Positive comparison: at helper occupancy5 the same actual51 free slots suffice.
    counts=Counter();section={}
    for ss in combinations(free51,3):
        q=plus(fixed,plus(plus(labels[ss[0]],labels[ss[1]]),labels[ss[2]]));counts[q]+=1;section.setdefault(q,ss)
    need(set(counts)==set(range(2048))and sum(counts.values())==20825,'same actual supply complete at occupancy5')
    phase=json.loads(phasepath.read_text())
    need(phase['reset_clock_bound']==43 and phase['reset_support_bound']==6,'inherited actual fixed-phase costs')
    cost={z['source']:z['selected']['cost']for z in supply['helpers']};cost[2125]=8
    need(sum(sorted((cost[s]for s in free51),reverse=True)[:6])==43,'actual six-free-slot clock bound')
    # Core's cheap-grade comparison requested by the parent: no extra native ports.
    corecost={z['source']:sum(a!=b for a,b in zip(z['old_before'],z['old_after']))for z in old};best={};witness={}
    for ss in combinations(oldids,4):
        q=0
        for s in ss:q=plus(q,labels[s])
        weight=sum(corecost[s]for s in ss)
        if q not in best or weight<best[q]:best[q]=weight;witness[q]=ss
    need(len(best)==2048 and dict(Counter(best.values()))=={8:1899,9:135,10:14},'exact weighted old28 grade4 section')
    cheapest_pair=min((cost[a]+cost[b],(a,b))for a,b in combinations(free51,2))
    largest_free_pair=sum(sorted((cost[s]for s in free51),reverse=True)[:2])
    largest_core_four=sum(sorted(corecost.values(),reverse=True)[:4])
    need(cheapest_pair==(8,(1882,1952))and largest_free_pair==15 and largest_core_four==15,'literal endpoint and reference-pair weight bounds')
    grade4_bound=largest_core_four+max(best.values())+largest_free_pair+cheapest_pair[0]
    need(grade4_bound==48,'grade4 endpoint-to-fixed-helper hub bound')
    def digest(p):
        b=p.read_bytes()
        if p.suffix in('.json','.tsv','.md','.py'):b=b.replace(b'\r\n',b'\n')
        return sha256(b).hexdigest()
    out=dict(status='FINITE-EXACT full native pair-capacity obstruction and occupancy5 comparison',
             normalized_words=len(owner),merged_source_classes=len(nodes),target_word_counts=dict(wordcounts),
             target_support_pattern_histogram=dict(Counter(support.values())),target_code_fibre=[W,U,V],
             fixed_control=dict(U=[368],V=[776,3139],shared_owner=4,leaf_charge=fixed),pair_envelopes=rows,
             universal_singleunit_pair_hostile=dict(residual_charge=37,total_helper_charge=1327),
             occupancy5=dict(free_sources=free51,triples=20825,charge_count=2048,min_multiplicity=min(counts.values()),max_multiplicity=max(counts.values()),
                             section={q:section[q]for q in sorted(section)},helper_reset_support=6,helper_reset_clock=43),
             core_grade4=dict(minimum_weight_histogram=dict(Counter(best.values())),maximum=10,
                              hardest={q:witness[q]for q in sorted(best)if best[q]==10},
                              maximum_endpoint_selected_weight=largest_core_four,reference_free_helper_pair=cheapest_pair[1],
                              reference_pair_weight=cheapest_pair[0],maximum_free_pair_weight=largest_free_pair,
                              endpoint_to_reference_hub_clock_bound=grade4_bound),checks=N,
             scope='Fixed physical source rows, one selected paired unit each, no other variable incidences at those rows; positive occupancy5 comparison changes the target margin and is not a legal occupancy4 repair.',
             input_hashes={p.relative_to(RESULTS).as_posix():digest(p)for p in(nodespath,ownerpath,electricpath,supplypath,oldpath,phasepath)})
    out=json.loads(json.dumps(out));path=HERE/'receipt.json'
    if '--write'in sys.argv:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',newline='\n')
    else:need(out==json.loads(path.read_text()),'frozen pair-capacity receipt')
    print(f'PASS {N}; full120 paired-target alphabet misses60 pair charges; actual occupancy5 is complete')
if __name__=='__main__':main()
