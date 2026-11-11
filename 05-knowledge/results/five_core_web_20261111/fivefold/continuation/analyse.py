"""Analyse the existing small-parent holes through their full H9 inverses."""
from collections import Counter, defaultdict
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[4]
source=ROOT/'05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('compiler',source)
c=importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
rows=json.loads((ROOT/'05-knowledge/results/chanpak_five_core_bridges_20261110/cover/signatures.json').read_text())


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def reverse(n,k):
    return int(f'{n:0{k}b}'[::-1],2)


def main():
    targets=json.loads((HERE/'targets.json').read_text())
    holes=set(targets['small_parent_holes'])
    witnesses=json.loads((HERE/'witnesses.json').read_text())
    by=defaultdict(list)
    high_pairs=[]
    internal_pairs=[]
    truncation_pairs=[]
    high_count=internal_count=truncation_count=0
    for h,ri,ci,a,b in witnesses:
        row=rows[ri]
        sig=row['signature']
        require(c.count(sig,a,b)==h,('inverse',h,ri,a,b))
        receipt={'row':ri,'core':row['core'],'cut':row['cut'],'A':a,'B':b,'parents':[]}
        for side,x in ((0,a),(1,b)):
            if x>=2:
                k=x.bit_length()-1
                n=x-(1<<k)
                aa,bb=(n,b) if side==0 else(a,n)
                ph=c.count(sig,aa,bb)
                high_count+=1
                item={'side':side,'k':k,'low_address':n,'H':ph,'jump':h-ph}
                receipt['parents'].append(item)
                if ph in holes:
                    high_pairs.append([h,ph,ri,a,b,side,k])
            for k in range(x.bit_length()):
                if x>>k&1:
                    aa,bb=(x-(1<<k),b) if side==0 else(a,x-(1<<k))
                    ph=c.count(sig,aa,bb)
                    internal_count+=1
                    if ph<h and ph in holes:
                        internal_pairs.append([h,ph,ri,a,b,side,k])
                aa,bb=(x%(1<<k),b) if side==0 else(a,x%(1<<k))
                ph=c.count(sig,aa,bb)
                truncation_count+=1
                if ph in holes:
                    truncation_pairs.append([h,ph,ri,a,b,side,k])
        by[h].append(receipt)
    require(set(by)==holes,'all inherited hole targets covered')
    short={h:min(min(w['A'],w['B']) for w in ws) for h,ws in by.items()}
    worst=max(short.values())
    # Independent all-bit derivatives and centered doubling controls.
    bit_checks=doubling_checks=0
    for k in range(1,11):
        for n in range(1<<k):
            D=(1<<k)-reverse(n,k)-(n>>(k-1))
            require(D>=c.g(n)-1,('doubling',k,n))
            equal=(n>>(k-1)==1 and f'{n:0{k}b}'.count('0')<=1)
            require((D==c.g(n)-1)==equal,('doubling equality',k,n))
            doubling_checks+=1
            for prefix in range(16):
                lo=(2*prefix)*(1<<k)+n
                hi=(2*prefix+1)*(1<<k)+n
                require(c.g(hi)-c.g(lo)==D-prefix-(prefix&1),('bit derivative',prefix,k,n))
                bit_checks+=1
    report={
      'status':'FINITE-EXACT correlation probes on inherited target universe; coverage OPEN',
      'target_count':len(holes),'largest_target':max(holes),'witness_count':len(witnesses),
      'inverse_universe':'all target-bounded H9 pairs from the44 inherited rows, modulo their flank-swap/signature quotient',
      'short_flank_minimum_histogram':dict(sorted(Counter(short.values()).items())),
      'worst_short_flank':worst,
      'worst_target_witnesses':{h:by[h] for h in sorted(by) if short[h]==worst},
      'rare_core_witnesses':{h:by[h] for h in targets['rare_S5a']+targets['rare_S5b']},
      'high_bit_comparisons':high_count,'high_bit_defect_pairs':high_pairs,
      'internal_bit_comparisons':internal_count,'internal_bit_defect_pairs':internal_pairs,
      'truncation_comparisons':truncation_count,'truncation_defect_incidence_count':len(truncation_pairs),
      'truncation_defect_examples':truncation_pairs[:12],
      'doubling_controls':doubling_checks,'internal_derivative_controls':bit_checks,
      'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
         for p in (source,ROOT/'05-knowledge/results/chanpak_five_core_bridges_20261110/cover/signatures.json',HERE/'targets.json',HERE/'witnesses.json')},
    }
    (HERE/'analysis.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in('worst_target_witnesses','rare_core_witnesses','source_sha256')},sort_keys=True))


if __name__=='__main__':
    main()
