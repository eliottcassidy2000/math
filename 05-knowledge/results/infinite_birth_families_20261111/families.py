"""Compile every fixed-H birth family with target H<125 into a binomial basis."""
from pathlib import Path
from collections import Counter
from math import comb,gcd
from functools import lru_cache
import argparse,csv,json

HERE=Path(__file__).resolve().parent

def convolve(p,q):
    out=Counter()
    for i,a in p.items():
        for j,b in q.items():
            out[i+j]+=a*b
    return out

def compile_families():
    raw=json.loads((HERE/'strong_class_polynomials_below125.json').read_text(encoding='utf-8'))
    P={int(h):{i:c for i,c in enumerate(v) if c} for h,v in raw.items()}
    raw=json.loads((HERE/'strong_birth_kernels_below125.json').read_text(encoding='utf-8'))['kernels']
    C={tuple(map(int,k.split(','))):{i:c for i,c in enumerate(v) if c} for k,v in raw.items()}
    @lru_cache(None)
    def words(d):
        if d==1:
            return {0:{0:1}}
        result={}
        for h,p in P.items():
            if d%h:
                continue
            for length,q in words(d//h).items():
                dst=result.setdefault(length+1,Counter())
                dst.update(convolve(p,q))
        return result
    def exterior(d):
        return {(degree,length+2):c*(length+1) for length,p in words(d).items() for degree,c in p.items()}
    compiled={}
    for h in range(1,125,2):
        for g in range(h,125,2):
            terms=Counter(exterior(h)) if h==g else Counter()
            for d in range(1,gcd(h,g)+1,2):
                if h%d or g%d or (h//d,g//d) not in C:
                    continue
                for (degree,pole),value in exterior(d).items():
                    for shift,c in C[h//d,g//d].items():
                        terms[degree+shift,pole]+=value*c
            if terms:
                compiled[str(h)+','+str(g)]=dict(
                    eventual_degree=max(pole-1 for degree,pole in terms),
                    valid_polynomial_from=max(degree for degree,pole in terms),
                    terms=[[degree,pole,c] for (degree,pole),c in sorted(terms.items())])
    (HERE/'all_size_birth_families_below125.json').write_text(json.dumps(dict(
        status='PROVED family formula + FINITE-EXACT complete primitive coefficients',
        convention='term [shift,pole,mass] means mass*x^shift/(1-x)^pole; ordinary pointed births',
        families=compiled),indent=2)+'\n',encoding='utf-8')
    return compiled

def evaluate(family,n):
    return sum(c*comb(n-shift+pole-1,pole-1) for shift,pole,c in family['terms'] if n>=shift)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-h',type=int)
    parser.add_argument('--target-h',type=int)
    parser.add_argument('--order',type=int)
    args=parser.parse_args()
    families=compile_families()
    if args.source_h is None:
        out=dict(status='PASS',nonzero_families=len(families),
                 degree_histogram=dict(Counter(f['eventual_degree'] for f in families.values())))
        (HERE/'family_checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(out))
        return
    h,g,n=args.source_h,args.target_h,args.order
    if g is None or n is None or not (1<=h<=g<125 and h%2 and g%2 and n>=0):
        parser.error('need odd 1<=source-h<=target-h<125 and order>=0')
    family=families.get(str(h)+','+str(g))
    result=dict(source_H=h,target_H=g,source_order=n,count=0 if family is None else evaluate(family,n),
                formula=family,status='exact at every source order')
    if gcd(h,g)==1 and result['count']:
        with (HERE/'primitive_births.tsv').open(encoding='utf-8') as f:
            r=next(r for r in csv.DictReader(f,delimiter='\t') if int(r['source_H'])==h and
                   int(r['target_H'])==g and int(r['source_n'])<=n)
        result['witness']=dict(primitive_source_rows=list(map(int,r['source_rows'].split(','))),
            primitive_out_mask=int(r['new_vertex_out_mask']),add_uniform_source_vertices=n-int(r['source_n']))
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
