"""Independent pair-count (not witness-bitset) audit of pair sections."""
from array import array
from collections import Counter
from itertools import combinations
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1];N=0
def need(v,m):
    global N
    N+=1
    if not v:raise RuntimeError(m)
def main():
    d=json.loads((RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json').read_text())
    r=json.loads((HERE/'pair_sections_receipt.json').read_text());ids=r['sources'];by={x['source']:x for x in d['helpers']}
    need(ids==next(b for b in d['banks']if len(b['source_ids'])==53)['source_ids'],'full authenticated53 bank')
    xs=[by[s]['label'][0]for s in ids];ys=[sum(b<<j for j,b in enumerate(by[s]['label'][1:]))for s in ids]
    pairs=list(combinations(range(53),2));index={p:i for i,p in enumerate(pairs)}
    totals=[0]*2048;singles=[array('H',[0])*53 for _ in range(2048)]
    counts=[array('H',[0])*1378 for _ in range(2048)]
    for ss in combinations(range(53),4):
        q=((sum(xs[i]for i in ss)%128)<<4)|(ys[ss[0]]^ys[ss[1]]^ys[ss[2]]^ys[ss[3]])
        totals[q]+=1
        for i in ss:singles[q][i]+=1
        for p in combinations(ss,2):counts[q][index[p]]+=1
    need(sum(totals)==292825 and min(totals)>0,'complete nonempty charge fibres')
    hist=Counter();safe=Counter();mixed=[]
    for k,(i,j)in enumerate(pairs):
        local=Counter();missing={m:[]for m in range(1,16)}
        for q in range(2048):
            both=counts[q][k];a=singles[q][i];b=singles[q][j]
            pattern=[totals[q]-a-b+both,b-both,a-both,both]
            need(min(pattern)>=0 and sum(pattern)==totals[q],'exact disjoint local pattern-count partition')
            available=sum(1<<j for j,n in enumerate(pattern)if n)
            need(available.bit_count()>=3,'at least three local patterns at every charge')
            hist[available]+=1;local[available]+=1
            for m in range(1,16):
                if not available&m:missing[m].append(q)
        good=[m for m in range(1,16)if not missing[m]]
        for m in good:safe[m]+=1
        for m in(2,4):
            if missing[m]:mixed.append(dict(pair=[ids[i],ids[j]],allowed_pattern='01'if m==2 else'10',missing=missing[m]))
        expected=dict(sources=[ids[i],ids[j]],available_pattern_histogram=dict(sorted(local.items())),
                      universally_safe_relation_masks=good,singleton11_covered=2048-len(missing[8]))
        need(json.loads(json.dumps(expected))==r['pair_rows'][k],'independent full pair row')
    need({str(k):v for k,v in hist.items()}==r['global_available_pattern_histogram'],'entire availability histogram')
    need({str(k):v for k,v in safe.items()}==r['safe_pair_counts_by_relation'],'all15 relation-mask counts')
    need(mixed==r['mixed_singleton_failures']and len(mixed)==32,'exact32 mixed singleton failures')
    # Hostile and repaired source-specific section, via a separate direct subset implementation.
    n=json.loads((HERE/'native_2125_receipt.json').read_text());alllabels={s:by[s]['label']for s in ids};alllabels[2125]=n['label']
    def charge(ss):
        low=0
        for s in ss:low^=sum(v<<j for j,v in enumerate(alllabels[s][1:]))
        return((sum(alllabels[s][0]for s in ss)%128)<<4)|low
    def accepts(ss):return((368 in ss)!=(776 in ss))and((368 in ss)!=(3139 in ss))
    images=[];repairs=[]
    for supply in(ids,ids+[2125]):
        image=set()
        for ss in combinations(supply,4):
            if accepts(ss):
                q=charge(ss);image.add(q)
                if len(supply)==54 and q==1509:repairs.append(list(ss))
        images.append(image)
    need(set(range(2048))-images[0]=={1509}and images[1]==set(range(2048)),'actual constrained53 hostile and54 repair')
    need(repairs==n['xor_tree']['repair_witnesses'],'both new1509 witnesses retained')
    # A branch-aware clock bound, with no free fixed-clock assumption for2125.
    costs={s:by[s]['selected']['cost']for s in ids};costs[2125]=8
    controls=[368,776,3139];free=sorted([v for s,v in costs.items()if s not in controls],reverse=True)
    clockbounds=[sum(costs[s]for s in controls)+sum(free[:5]),sum(free[:6]),sum(free[:4])]
    need([costs[s]for s in controls]==[6,5,6]and clockbounds==[53,43,29],'actual phase-specific53 helper clock bound')
    out=dict(status='FINITE-EXACT independent count-partition audit',checks=N,pairs=len(pairs),quadruples=sum(totals),
             xor_tree_missing=[1509],new_charge_count=len(images[1]),clock_cases=dict(cross_phase=clockbounds[0],center_V=clockbounds[1],center_U=clockbounds[2]),
             helper_support_bound=8,helper_clock_bound=max(clockbounds))
    out=json.loads(json.dumps(out));path=HERE/'check_receipt.json'
    if '--write'in sys.argv:path.write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    else:need(out==json.loads(path.read_text()),'independent frozen receipt')
    print(f'PASS {N:,}; all1378 pair predicates, exact XOR-tree repair, clock bound{max(clockbounds)}')
if __name__=='__main__':main()
