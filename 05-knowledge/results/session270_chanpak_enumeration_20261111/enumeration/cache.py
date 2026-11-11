"""Exact cached A000568 query; preprocessing and bounded support are explicit."""
import argparse,hashlib,json,subprocess,sys,time,types,statistics
from pathlib import Path
from counter import count
BASE='f2e03e29358d7a6e566f27f75794ec601d0449f5'
SOURCES=['05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py',
'05-knowledge/results/session239_sequence_engines_20261110/burnside/engines.py']
HERE=Path(__file__).resolve().parent
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)

def load():
    raw=(HERE/'prefix.json').read_bytes();meta=json.loads((HERE/'cache_receipt.json').read_text())
    if hashlib.sha256(raw).hexdigest()!=meta['file_sha256']:raise ValueError('cache checksum mismatch')
    return [int(x) for x in json.loads(raw)['values']]

def build(n):
    modules=[];hashes={}
    for path in SOURCES:
        raw=subprocess.check_output(['git','show',BASE+':'+path]);hashes[path]=hashlib.sha256(raw).hexdigest()
        mod=types.ModuleType(path);exec(compile(raw,path,'exec'),mod.__dict__);modules.append(mod)
    t=time.perf_counter();a,_=modules[0].local_dp(n);a=[int(x) for x in a];first=time.perf_counter()-t
    t=time.perf_counter();b,_=modules[1].integer_response(n);second=time.perf_counter()-t
    if a!=b:raise AssertionError('whole-prefix independent recurrence disagreement')
    raw=(json.dumps({'sequence':'A000568','max_n':n,'values':[str(x) for x in a]},indent=2)+'\n').encode()
    (HERE/'prefix.json').write_bytes(raw)
    receipt={'base_commit':BASE,'source_sha256':hashes,'max_n':n,'all_coefficients_compared':len(a),
      'degree_local_seconds':first,'independent_factorial_seconds':second,
      'file_sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'scope':'Precomputed exact table, not a faster general enumeration algorithm.'}
    (HERE/'cache_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    timings=[]
    for _ in range(5):
        t=time.perf_counter();arr=load();timings.append(time.perf_counter()-t)
    t=time.perf_counter()
    for i in range(100000): x=arr[i%len(arr)]
    receipt['cold_read_verify_parse_seconds']=timings
    receipt['warm_lookup_seconds_per_query']=(time.perf_counter()-t)/100000
    (HERE/'cache_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

def time_cache():
    receipt=json.loads((HERE/'cache_receipt.json').read_text());times=[]
    for _ in range(5):
        t=time.perf_counter()
        output=subprocess.check_output([sys.executable,str(Path(__file__).resolve()),str(receipt['max_n'])])
        times.append(time.perf_counter()-t)
    receipt['cold_cli_max_query_seconds']=times
    receipt['cold_cli_max_query_median']=statistics.median(times)
    receipt['max_query_output_digits']=len(output.strip())
    (HERE/'cache_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))

def main():
    p=argparse.ArgumentParser();p.add_argument('n',type=int,nargs='?');p.add_argument('--build',type=int);p.add_argument('--fresh',action='store_true');p.add_argument('--time-cache',action='store_true');args=p.parse_args()
    if args.time_cache:time_cache();return
    if args.build is not None:build(args.build);return
    if args.n is None or args.n<0:p.error('provide nonnegative n')
    if args.fresh:answer,_=count(args.n)
    else:
        values=load()
        if args.n<len(values):answer=values[args.n]
        else:answer,_=count(args.n)
    print(answer)
if __name__=='__main__':main()
