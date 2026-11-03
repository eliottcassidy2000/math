#!/usr/bin/env python3
"""Reproduce from pinned Git objects; compare normal/-O and C++ -O0/-O2."""
from pathlib import Path
from collections import defaultdict
import hashlib,json,shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def run(args,**kw):
    return subprocess.run([str(x) for x in args],check=True,**kw)

def main():
    for filename in ['SOURCES.json','INCOMING_SOURCES.json']:
        pin=json.loads((HERE/filename).read_text())
        for item in pin['sources']:
            raw=subprocess.check_output(['git','show',pin['public_commit']+':'+item['path']],cwd=HERE)
            need(hashlib.sha256(raw).hexdigest()==item['sha256'],'pinned source: '+item['path'])
    p=json.loads((HERE/'ports.json').read_text());sources=[];signs=defaultdict(set)
    for order in p['orders']:
        for v in order['strong_below125']:
            sources.append([order['n'],v['H'],*v['rows'],*v['I']])
            signs[v['H']].update(x for x in v['I'] if abs(x)==1)
    need(sources==[list(map(int,l.split())) for l in (HERE/'low_sources.txt').read_text().splitlines()],'low-source extraction')
    summary=dict(status='FINITE-EXACT',source_scope='All orders, by the proved SCC reduction and f(10)=125',strong_sources_below125=len(sources),
                 missing_counts_below125=[h for h in range(3,125,2) if h not in signs],
                 missing_signs=[dict(H=h,missing=[s for s in [-1,1] if s not in signs[h]]) for h in sorted(signs) if len(signs[h])<2],
                 no_binary_ports=[h for h in sorted(signs) if not signs[h]],
                 two_sided_class_counts={str(o['n']):o['both'] for o in p['orders']},
                 both_children_two_sided=sum(o['both_children_both'] for o in p['orders']))
    need(summary==json.loads((HERE/'BINARY_SUMMARY.json').read_text()),'summary derivation')
    with tempfile.TemporaryDirectory(prefix='metagraph-primegame-replay-') as temp:
        tmp=Path(temp);inputs=tmp/'inputs'
        run([sys.executable,HERE/'prepare_inputs.py',inputs])
        shutil.copyfile(HERE/'primegame_code.py',tmp/'primegame_code.py')
        for flags in [[],['-O']]:
            run([sys.executable,*flags,tmp/'primegame_code.py'])
            need((tmp/'primegame_code.json').read_bytes()==(HERE/'primegame_code.json').read_bytes(),'Python optimized parity')
        for flags in [[],['-O']]:
            run([sys.executable,*flags,HERE/'audit_induced.py',tmp/'induced.json'])
            need((tmp/'induced.json').read_bytes()==(HERE/'AUDIT_INDUCED.json').read_bytes(),'induced literal enumeration normal/-O')
        for opt in ['-O2','-O0']:
            print('Checking C++',opt,flush=True)
            for program,receipt in [('ports','ports.json'),('audit_ports','AUDIT_PORTS.json'),('word_search','WORD_SEARCH.json'),('audit_two_birth','AUDIT_TWO_BIRTH.json'),('induced_branches','INDUCED_BRANCHES.json'),('padded_ears','PADDED_EARS.json')]:
                binary=tmp/program
                run(['c++',opt,'-std=c++17','-Wall','-Wextra','-Werror',HERE/(program+'.cpp'),'-o',binary])
                if program=='ports':out=subprocess.check_output([str(binary),str(inputs),'9'])
                else:out=subprocess.check_output([str(binary)],input=(HERE/'low_sources.txt').read_bytes())
                need(out==(HERE/receipt).read_bytes(),program+' receipt '+opt)
                print('PASS',program,opt,flush=True)
    print('PASS all source pins, exact receipts, normal/-O and O0/O2')
if __name__=='__main__':main()
