"""Costed99-ray recursion and exact small-prime wheel interfaces."""
from pathlib import Path
import hashlib,json,importlib.util
import sympy as s
import numpy as np
from census_nine import need,save,tables
from construct import tail
HERE=Path(__file__).resolve().parent

def run():
    source=HERE.parent/'path_spectrum_two_angles_metagraph_20261106/data/least_order_11.txt'
    least={}
    for line in source.read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split()
        if f[1].isdigit():least[int(f[0])]=int(f[1])
    wheelpath=HERE.parent/'endpoint_coverage_compiler_20261111/wheel.json';wheel=json.loads(wheelpath.read_text())
    primes=list(s.primerange(2,100));need(len(primes)==25 and primes==wheel['primes_below_100'],'first25 primes')
    results=[]
    for k in range(1,31):
        N=99*2**k+23;small=[];allf=[]
        for p in primes:
            if N%p==0:
                small.append({'p':p,'cofactor':N//p,'join_order':least[p]+least[N//p] if p in least and N//p in least else None})
        for d in s.divisors(N):
            d=int(d)
            if d==1 or d*d>N:continue
            if d in least and N//d in least:allf.append({'factors':[d,N//d],'join_order':least[d]+least[N//d]})
        allf.sort(key=lambda r:(r['join_order'],r['factors']))
        need(bool(small)==any(k%a['period']==a['residue'] for a in wheel['active']),'literal factor versus exponent-state wheel')
        if k<=9:
            rr=tail([40,17,3,6,13,30],[0,3],k);need(int(tables(np.array(rr,dtype=np.int64))[0])==N,'literal physical ray')
        results.append({'k':k,'H':N,'factorization':[[int(p),int(e)] for p,e in s.factorint(N).items()],'direct_tail_order':k+7,'least_order_in_complete11_spectrum':least.get(N),'small_prime_join_options':small,'best_known_two_factor_join':allf[0] if allf else None})
    need(results[0]['direct_tail_order']==8 and results[0]['best_known_two_factor_join']['join_order']==11,'factorization can cost more')
    need(not results[15]['small_prime_join_options'] and results[15]['best_known_two_factor_join']['join_order']==19,'small-prime wheel misses cheaper large-factor join')
    need(all((2-a['residue'])%a['period'] for a in wheel['active']),'wheel survivor419')
    out={'scope':'Exact physical99*2^k+23 family, k>=1. Cost comparison for k1..30 against the inherited complete least-order bank through11. Absence from the cost bank is unknown, not impossible. Prime wheel is only an algorithmic routing state, not a primality or Hamilton-realization criterion.',
         'source_lf_sha256':hashlib.sha256(source.read_text().encode()).hexdigest(),'wheel_lf_sha256':hashlib.sha256(wheelpath.read_text().encode()).hexdigest(),'wheel':wheel,'records':results}
    save('natural_recursion.json.gz',out);print(json.dumps(results[:17]),flush=True)
if __name__=='__main__':run()
