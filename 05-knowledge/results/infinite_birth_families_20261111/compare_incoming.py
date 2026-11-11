"""Compare independently developed concurrent arithmetic engines at our bounds."""
from pathlib import Path
import csv,hashlib,json,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'uniform_return_metagraph_20261111/arithmetic'
with tempfile.TemporaryDirectory(prefix='birth-arithmetic-crosscheck-') as td:
    td=Path(td)
    subprocess.run([sys.executable,'-B',str(OTHER/'maxima.py'),'--limit','10000',
                    '--table',str(td/'max.txt'),'--output',str(td/'max.json')],
                   check=True,stdout=subprocess.DEVNULL)
    subprocess.run([sys.executable,'-B',str(OTHER/'quadratic_primes.py'),'--root-bound','100000',
                    '--table',str(td/'prime.txt'),'--output',str(td/'prime.json')],
                   check=True,stdout=subprocess.DEVNULL)
    incoming_max=dict(tuple(map(int,line.split())) for line in (td/'max.txt').read_text().splitlines())
    own_max={int(r['n']):int(r['largest_predecessor']) for r in csv.DictReader(
             (HERE/'A135244_10000.tsv').open(encoding='utf-8'),delimiter='\t')}
    incoming_prime=dict(tuple(map(int,line.split())) for line in (td/'prime.txt').read_text().splitlines())
    own_prime={int(r['index']):int(r['prime']) for r in csv.DictReader(
               (HERE/'A002496_inputs_100000.tsv').open(encoding='utf-8'),delimiter='\t')}
    if incoming_max!=own_max or incoming_prime!=own_prime:
        raise RuntimeError('concurrent arithmetic engine mismatch')
result=dict(status='PASS',inverse_entries=len(own_max),quadratic_prime_entries=len(own_prime),
    source_packet='uniform_return_metagraph_20261111/arithmetic',
    source_hashes={p:hashlib.sha256((OTHER/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                   for p in ['maxima.py','quadratic_primes.py']})
(HERE/'incoming_comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))
