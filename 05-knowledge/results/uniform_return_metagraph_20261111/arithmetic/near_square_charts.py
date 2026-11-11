"""Fresh chart witnesses for the inherited first-three-core near-square support.

All-parameter compiler and inverse bounds are inherited. The new finite
statistic is minimum short flank in that normalized grammar.
"""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import csv,importlib.util,json,subprocess,tempfile

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
PRIOR=RESULTS/'chanpak_five_core_bridges_20261110'
COMPILER=RESULTS/'session251_creative_completion_20261111/ternary/compiler.py'
CORES=((0,),(2,4,1),(2,4,9,3))


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def main():
    records=json.loads((PRIOR/'cover/positive_raw_signatures.json').read_text())
    text=str(len(records))+'\n'+''.join(' '.join(map(str,[r['core_index'],r['cut'],*r['positive_observer']]))+'\n' for r in records)
    with tempfile.TemporaryDirectory(prefix='square-charts-') as td:
        t=Path(td)
        (t/'raw.txt').write_text(text)
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'square_charts.cpp'),'-o',str(t/'run.exe')],check=True,capture_output=True)
        result=subprocess.run([str(t/'run.exe'),str(t/'raw.txt'),'16777215',str(t/'table.tsv')],check=True,capture_output=True,text=True)
        data=(t/'table.tsv').read_bytes()
    table=[{k:int(v) for k,v in row.items()} for row in csv.DictReader(data.decode().splitlines(),delimiter='\t')]
    need(len(table)==2047 and all(r['core']>=0 for r in table),'all near-square targets')
    spec=importlib.util.spec_from_file_location('inherited_compiler',COMPILER)
    c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
    literal=0
    for r in table:
        sig=c.signature(CORES[r['core']],r['cut'])
        need(c.count(sig,r['A'],r['B'])==r['value'],'independent signed compiler evaluation')
        if r['order']<=13:
            graph=c.construct(CORES[r['core']],r['cut'],r['A'],r['B'])
            need(c.hp_subsets(graph)[-1]==r['value'],'literal graph DP')
            literal+=1
    (HERE/'near_square_charts.tsv').write_bytes(data)
    nonaxis=[r for r in table if r['min_flank']>0]
    old=json.loads((PRIOR/'cover/support_24.json').read_text())
    need(not set(old['small_parent_extra_holes'])&{r['value'] for r in table},'inherited support intersection')
    receipt={'status':'FINITE-EXACT; all-height near-square coverage OPEN',
        'target_bound':16777215,'root_domain':'positive even roots2..4094',
        'target_count':len(table),'min_flank_distribution':dict(Counter(r['min_flank'] for r in table)),
        'first_nonaxis':nonaxis[:10],'max_min_flank':max(r['min_flank'] for r in table),
        'literal_DP_checks_order_at_most13':literal,
        'signed_formula_checks':len(table),'search':result.stdout.strip(),
        'table_sha256':sha256(data.replace(b'\r\n',b'\n')).hexdigest(),
        'inherited_sources':{str(p.relative_to(RESULTS)):sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [COMPILER,PRIOR/'cover/positive_raw_signatures.json',PRIOR/'cover/support_24.json']}}
    out=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    (HERE/'near_square_charts.json').write_text(out,encoding='utf-8')
    print(out)


if __name__=='__main__':main()
