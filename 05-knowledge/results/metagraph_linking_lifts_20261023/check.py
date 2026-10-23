"""Independent sparse consumer of the balanced-switch literal lift certificates.

No producer code is imported. The pinned earlier independent consumer supplies
highest-pivot nullspaces and sparse row multiplication, not correction search.
"""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import argparse,copy,csv,gzip,importlib.util,json

HERE=Path(__file__).resolve().parent
API=HERE.parent/'metagraph_allsizes_20261023/dyadic_saved_consumer.py'
spec=importlib.util.spec_from_file_location('consumer',API)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)

def read(path):
    raw=path.read_bytes() if path.exists() else gzip.decompress(Path(str(path)+'.gz').read_bytes())
    return raw.replace(b'\r\n',b'\n')

def all_layers(A,claim):
    V1,rankA,Q,rankQ,expected=c.independent_half_form(A)
    c.need(len(V1)==claim['critical_two_rank'] and rankQ==claim['half_form_rank'],'independent full first-layer agreement')
    layers=[{'level':1,'dimension':len(V1),'rank':rankQ,'next_dimension':len(expected)}]
    saved=[(2,claim['radical_lifts_mod4'],claim['quarter_form_rows']),
           (3,claim['eighth_radical_lifts_mod8'],claim['eighth_form_rows'])]
    saved += [(x['level'],x['lifts'],x['form_rows'])for x in claim['higher_layers']]
    c.need([x[0]for x in saved]==list(range(2,len(saved)+2)),'consecutive certificate layers')
    last_lift=None
    for level,lifts,form in saved:
        expected,entry=c.consume_layer(A,claim['vertices'],claim['root'],level,lifts,form,expected)
        layers.append(entry)
        if entry['rank']:last_lift=(level,lifts[0],entry['rank'],len(lifts))
    c.need(not expected and claim['complete_two_primary'],'termination without leftover modes')
    e,x,last_rank,last_dimension=last_lift
    c.need(last_rank==last_dimension==1,'highest constituent is one-dimensional in these four cases')
    image=c.sparse_image(A,x[:-1]);energy=sum(a*b for a,b in zip(x[:-1],image));modulus=1<<e
    c.need(all(v%modulus==0 for v in image) and any(v%2 for v in x),'primitive highest-order charge')
    c.need(energy%modulus==0 and (energy//modulus)%2==1,'nondegenerate self-linking unit')
    multiplicities={str(1<<row['level']):row['rank']for row in layers if row['rank']}
    return {'layers':layers,'critical_two_primary':multiplicities,
            'tree_valuation':sum(row['level']*row['rank']for row in layers),
            'highest_charge':{'order':modulus,'energy':energy,'linking_numerator':energy//modulus%modulus}},Q,V1

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    base=HERE.parent/'metagraph_order8_electrical_20261023';graph=HERE.parent/'metagraph_order8_20261023'
    raw=read(HERE/'certificates.json');cert=json.loads(raw);receipt=json.loads(read(base/'receipt.json'))
    for key,path in [('nodes',graph/'merged_nodes.tsv'),('edges',graph/'merged_edges.tsv'),('witness_receipt',graph/'receipt.json'),('native_electrical_receipt',base/'receipt.json'),('producer_code',base/'check.py')]:
        c.need(sha256(read(path)).hexdigest()==cert['input_sha256_lf'][key],'pinned input '+key)
    c.need(sha256(read(API)).hexdigest()=='4668a9e27a965fbc2901848da28f793b558138bf86cc7cd5309bd504599fcf03','independent consumer API pinned')
    A,n,root,inputs=c.read_graph(receipt)
    witness=json.loads(read(graph/'receipt.json'))['one_round_counterfeit'];edits=witness['edits']
    c.need(all(type(e['delta']) is int and e['delta'] in (-1,1) and type(e['native_count']) is int and e['native_count']>=0 for e in edits),'literal unit port edits')
    lo=max(-e['native_count']for e in edits if e['delta']==1);hi=min(e['native_count']for e in edits if e['delta']==-1)
    c.need([lo,hi]==[-1,2]==cert['feasible_integer_interval'],'exact port feasibility interval')
    rows=Counter();cols=Counter()
    for e in edits:
        a,b=e['source_endpoint_scores'];rows[e['source_class'],a,b]+=e['delta'];cols[e['target_class'],8-a,6-b]+=e['delta']
    c.need(not any(rows.values()) and not any(cols.values()),'original source and target endpoint capacities retained')
    c.need(all(e['native_count']+(lo-1)*e['delta']>=0 for e in edits) is False and all(e['native_count']+(hi+1)*e['delta']>=0 for e in edits) is False,'adjacent outside multiples infeasible')
    a,b,u,v=1238,2242,2246,2315
    delta={(a,u):2,(a,v):-2,(b,u):-2,(b,v):2}
    blue_vertices=set()
    for e in csv.DictReader(read(graph/'merged_edges.tsv').decode().splitlines(),delimiter='\t'):
        if int(e['blue']) and e['source']!=e['target']:blue_vertices.update((int(e['source']),int(e['target'])))
    c.need(not ({a,b,u,v}&blue_vertices),'four switch terminals blue-isolated')
    results=[];nativeQ=None;nativeV=None;calcs={x['switch_multiple']:x['calculation']for x in cert['cases']}
    expected_groups={-1:{'2':3360,'4':2,'32':1},0:{'2':3358,'4':4,'256':1},1:{'2':3360,'4':2,'64':1},2:{'2':3358,'4':4,'8':1}}
    for t in (0,-1,1,2):
        c.need(all(e['native_count']+t*e['delta']>=0 for e in edits),'scaled original port edits nonnegative')
        changed=copy.deepcopy(A)
        for (i,j),weight in delta.items():
            changed[i][j]=changed[i].get(j,0)+t*weight;changed[j][i]=changed[j].get(i,0)+t*weight
        c.need(all(sum(row.values())==sum(A[i].values()) and row.get(i,0)==A[i].get(i,0)for i,row in enumerate(changed)),'row sums and degrees retained')
        c.need(all(w<=0 for i,row in enumerate(changed)for j,w in row.items()if i!=j),'nonnegative graph conductances')
        adjacency=[set()for _ in range(n)]
        for i,row in enumerate(changed):
            for j,w in row.items():
                if i!=j and w<0:adjacency[i].add(j)
            if sum(row.values())>0:adjacency[i].add(root);adjacency[root].add(i)
        reached={root};stack=[root]
        while stack:
            for j in adjacency[stack.pop()]:
                if j not in reached:reached.add(j);stack.append(j)
        c.need(len(reached)==n,'each feasible weighted graph remains connected')
        out,Q,V=all_layers(changed,calcs[t]);out['switch_multiple']=t
        c.need(out['critical_two_primary']==expected_groups[t],'independently established primary group')
        if t==0:nativeQ=Q;nativeV=V
        else:
            c.need(V==nativeV,'same marked binary harmonic basis')
            expected_diff=[]
            for x in V:
                U=((x>>a)^(x>>b))&1;W=((x>>u)^(x>>v))&1
                expected_diff.append(sum(((t*(U*(((y>>u)^(y>>v))&1)+W*(((y>>a)^(y>>b))&1)))%2)<<j for j,y in enumerate(V)))
            c.need([x^y for x,y in zip(Q,nativeQ)]==expected_diff,'exact low-rank marked half-form switch')
            c.need(all(not((row>>i)&1)for i,row in enumerate(expected_diff)),'every diagonal half-energy unchanged')
            out['marked_half_form_changed']=Q!=nativeQ
            out['marked_half_form_difference_rank']=len(c.echelon(expected_diff))
        results.append(out)
    c.need(calcs[2]['radical_lifts_mod4']==calcs[0]['radical_lifts_mod4'],'same actual mod-four lifts for doubled switch')
    diff=[x^y for x,y in zip(calcs[2]['quarter_form_rows'],calcs[0]['quarter_form_rows'])]
    c.need(diff==[4,0,9,4,0] and len(c.echelon(diff))==2,'marked quarter-form changes while total rank stays four')
    # Hostile control: a diagonal-only comparison accepts all three wrong networks.
    c.need(((A[a].get(u,0)//2)%2)!=((A[a].get(u,0)+2)//2)%2,'direct two-terminal half-energy witness')
    result={'status':'FINITE-EXACT independent actual-port switch certificate consumer',
            'input_sha256_lf':inputs,'certificate_sha256_lf':sha256(raw).hexdigest(),
            'consumer_api_sha256_lf':sha256(read(API)).hexdigest(),'checks':c.CHECKS,
            'feasible_integer_interval':[lo,hi],'cases':sorted(results,key=lambda x:x['switch_multiple']),
            'doubled_switch_quarter_difference_rows':diff}
    output=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    if args.write:(HERE/'receipt.json').write_bytes(output)
    else:c.need(read(HERE/'receipt.json')==output,'frozen independent receipt')
    print(json.dumps({'checks':c.CHECKS,'cases':[{k:v for k,v in x.items()if k!='layers'}for x in result['cases']]},indent=2))

if __name__=='__main__':main()
