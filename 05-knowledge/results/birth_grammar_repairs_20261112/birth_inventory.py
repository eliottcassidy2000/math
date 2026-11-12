"""Unions and extension-ancestry statistics, with receipt-scoped identifiers."""
from common import *
from synchrony_probe import expand
from collections import Counter,defaultdict

def run():
    byorder=defaultdict(set);extensions=[]
    for r in read('birth_union_ten.json.gz',True)['records']:byorder[len(r['rows'])].add(key(r['rows'],r['selected']))
    paths=['birth_extensions_11.json.gz','birth_extensions_12.json.gz','birth_extensions_13.json.gz','birth_extensions_13_synchrony.json.gz','birth_extensions_14_rank5.json.gz']
    for name in paths:
        d=read(name);indegree=Counter();outdegree=Counter();edge_marks=Counter();edges=0
        for r in d['records']:
            byorder[len(r['rows'])].add(tuple(r['canonical_key']));parents={str(p[0]) for p in r['presentations']};indegree[len(parents)]+=1;edges+=len(parents)
            for p in parents:outdegree[p]+=1
            for p in {(str(p[0]),p[2]) for p in r['presentations']}:edge_marks[p[1]]+=1
        out_hist=Counter(outdegree.values());out_hist[0]=d['parent_births']-len(outdegree)
        extensions.append({'receipt':name,'parent_order':d['parent_order'],'parent_bank':d['parent_births'],'children':d['distinct_births'],'ordinary_parent_child_edges':edges,'child_indegree_histogram':dict(sorted(indegree.items())),'parent_outdegree_histogram':dict(sorted(out_hist.items())),'parent_child_mark_edges':dict(edge_marks),'rank_histogram':d['stable_port_histogram'],'scope':d['scope']})
    sync={}
    for m in [2,3]:
        keys=set()
        for r in read('synchrony_kernels.json.gz')['records']:
            rr,O,_=expand(r['rows'],r['selected'],dict.fromkeys(r['ports'][:2],m));keys.add(key(rr,O))
        n=8+2*m;sync[n]={'classes':len(keys),'overlap_with_extension_bank':len(keys&byorder[n])};byorder[n]|=keys
    a={tuple(r['canonical_key']) for r in read('birth_extensions_13.json.gz')['records']};b={tuple(r['canonical_key']) for r in read('birth_extensions_13_synchrony.json.gz')['records']}
    need(not a&b,'disjoint mixed and plain order13 layers');need(sum(map(len,byorder.values()))==4042,'declared retained through14 total')
    out={'scope':'Union of complete inherited births through10, the three full declared extension layers11..13, synchronized slices12 and14, all one-vertex children of synchronized12, and order14 children only of the eight retained five-port13 parents. Complete after10 only within those explicitly listed construction universes. Birth IDs are local to each receipt; cross-file identity always uses actual marked canonical keys.','birth_counts_by_order':{n:len(v) for n,v in sorted(byorder.items())},'retained_births_through14':sum(map(len,byorder.values())),'synchrony_slice_intersections':sync,'order13_plain_mixed_intersection':len(a&b),'extension_layers':extensions}
    save('birth_inventory.json.gz',out);print(json.dumps({k:v for k,v in out.items() if k!='extension_layers'}),flush=True)
if __name__=='__main__':run()
