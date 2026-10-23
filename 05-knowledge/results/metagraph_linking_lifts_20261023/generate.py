"""Generate literal lifts for all feasible multiples of the native port switch.

The producer imports the frozen audited electrical engine. The separate
check.py only consumes its certificates, using an independent sparse engine.
"""
from pathlib import Path
from hashlib import sha256
import argparse,csv,gzip,importlib.util,json

HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'metagraph_order8_electrical_20261023'
spec=importlib.util.spec_from_file_location('producer',BASE/'check.py')
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    graph=HERE.parent/'metagraph_order8_20261023'
    nr=p.raw_read(graph/'merged_nodes.tsv');er=p.raw_read(graph/'merged_edges.tsv')
    wr=p.raw_read(graph/'receipt.json');br=p.raw_read(BASE/'receipt.json')
    witness=json.loads(wr)['one_round_counterfeit'];n=len(list(csv.DictReader(nr.decode().splitlines(),delimiter='\t')))
    L=[];B={}
    for row in csv.DictReader(er.decode().splitlines(),delimiter='\t'):
        u,v,blue,black=(int(row[k]) for k in ('source','target','blue','black'))
        p.need(black%2==0,'literal paired-black count')
        if blue:L.append((u,v,blue))
        if black:B[u,v]=black//2
    lower=max(-e['native_count'] for e in witness['edits'] if e['delta']==1)
    upper=min(e['native_count'] for e in witness['edits'] if e['delta']==-1)
    p.need((lower,upper)==(-1,2),'exact feasible integer port interval')
    cases=[]
    for t in range(lower,upper+1):
        p.need(all(e['native_count']+t*e['delta']>=0 for e in witness['edits']),'port capacities nonnegative')
        if t==0:out=json.loads(br)['order8']
        else:
            changed=B.copy()
            for (u,v),weight in witness['merged_black_delta']:changed[u,v]+=t*weight//2
            p.need(all(w>=0 for w in changed.values()),'graph conductances nonnegative')
            out=p.calculate(n,L,[(u,v,w) for (u,v),w in changed.items() if w])
        p.need(out['complete_two_primary'],'complete finite filtration')
        cases.append({'switch_multiple':t,'calculation':out})
        print(t,[(1,out['half_form_rank']),(2,out['quarter_form_rank']),(3,out['exact_valuation_three'])]+[(q['level'],q['rank'])for q in out['higher_layers']],flush=True)
    result={'status':'FINITE-EXACT literal lift certificates for actual-margin counterfeits',
      'input_sha256_lf':{k:sha256(raw).hexdigest() for k,raw in [('nodes',nr),('edges',er),('witness_receipt',wr),('native_electrical_receipt',br),('producer_code',p.raw_read(BASE/'check.py'))]},
      'feasible_integer_interval':[lower,upper],'cases':cases,'checks':p.CHECKS}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode()
    if args.write:(HERE/'certificates.json.gz').write_bytes(gzip.compress(raw,mtime=0))
    else:p.need(p.raw_read(HERE/'certificates.json')==raw,'frozen producer certificate')
    print('certificate_sha256',sha256(raw).hexdigest())

if __name__=='__main__':main()
