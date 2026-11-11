"""Reproduce bounded-short-flank witnesses for the unresolved +4 dyadic ray.

These are target witnesses, not a path from a common source or a uniform law.
"""
from pathlib import Path
import importlib.util,json
import chart_inverse as inv

HERE=Path(__file__).resolve().parent


def main():
    source=inv.RESULTS/'session251_creative_completion_20261111/ternary/compiler.py'
    spec=importlib.util.spec_from_file_location('signed',source)
    oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)
    rows=[]
    for k in range(1,26):
        target=102*2**k+19
        stats={}
        witnesses=inv.bank_inverse(target,max_short=16,first=True,stats=stats)
        inv.require(bool(witnesses),'bounded probe lacks a witness')
        for w in witnesses:
            sig=oracle.signature(inv.CORES[w['core']],w['cut'])
            inv.require(oracle.count(sig,w['A'],w['B'])==target,'signed witness check')
            w['A_bits']=bin(w['A'])[2:]
            w['B_bits']=bin(w['B'])[2:]
        rows.append({'k':k,'target':target,'witnesses':witnesses,'work':stats})
    result={'status':'FINITE-EXACT target witnesses; no common source-return or uniform binary law claimed',
            'short_flank_cap':16,'rows':rows}
    text=json.dumps(result,indent=2)+'\n'
    (HERE/'ternary_plus_four.json').write_text(text,encoding='utf-8')
    print(json.dumps({'status':result['status'],'rows':len(rows),'largest_target':rows[-1]['target']}))


if __name__=='__main__':main()
