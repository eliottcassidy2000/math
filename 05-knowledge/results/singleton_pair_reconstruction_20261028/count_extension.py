"""Extend inherited exact scalar counts; never generate a large metagraph."""
import argparse, hashlib, importlib.util, json, sys
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def main():
    sys.set_int_max_str_digits(0)
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=350)
    ap.add_argument('--reference',action='store_true');args=ap.parse_args()
    start=perf_counter()
    if args.reference:
        module=load('tail','05-knowledge/results/strong_attachment_tail_compiler_20261028/counting.py')
        values,stats=module.contracted(args.bound,progress=True);name='reference'
    else:
        module=load('bucket','05-knowledge/results/session87_four_vertices_enumeration_20261028/fast_counts.py')
        values,stats=module.degree_bucket_dp(args.bound);name='bucket'
    elapsed=perf_counter()-start
    prior=json.loads((HERE.parent/'strong_attachment_tail_compiler_20261028/ordinary-300.json').read_text())['values']
    if values[:min(len(values),len(prior))]!=prior[:len(values)]:raise RuntimeError('old prefix disagrees')
    strong=[0]*len(values)
    for n in range(1,len(values)):
        strong[n]=values[n]-sum(strong[k]*values[n-k] for k in range(1,n))
    digest=hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest()
    result=dict(status='FINITE-EXACT',method=name,bound=args.bound,values=values,
                nonempty_strong=strong,seconds=elapsed,statistics=stats,values_sha256=digest)
    path=HERE/f'counts-{name}-{args.bound}.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(path=str(path),seconds=elapsed,final_digits=len(str(values[-1])),values_sha256=digest)),flush=True)

if __name__=='__main__':main()
