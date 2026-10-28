"""Verify sealed receipts and pinned sources; replay exact computations on request."""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SCRIPTS=['fields','carriers','counting','transport','audit_fields']
PREFIX={'audit_fields':'audit-fields'}


def need(ok,message):
    if not ok: raise RuntimeError(message)


def digest(data): return sha256(data.replace(b'\r\n',b'\n')).hexdigest()
def read(name): return json.loads((HERE/name).read_text(encoding='utf-8'))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--replay',action='store_true')
    parser.add_argument('--replay-metagraph',action='store_true'); args=parser.parse_args()
    names=[]
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1)
        need(digest((HERE/name).read_bytes())==expected,'hash '+name); names.append(name)
    need(set(names)=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'},'inventory')
    for item in read('SOURCES.json')['repository_sources']:
        content=subprocess.check_output(['git','cat-file','blob',item['git_blob']],cwd=ROOT)
        need(digest(content)==item['sha256_lf'],'source '+item['path'])
    for p in HERE.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
            if '://' not in link and not link.startswith('#'):
                need((p.parent/link.split('#')[0]).resolve().exists(),'link '+link)
    for name in SCRIPTS:
        prefix=PREFIX.get(name,name)
        need(read(prefix+'-normal.json')==read(prefix+'-optimized.json'),'normal/O '+name)
        if args.replay:
            for flags in ([],['-O']):
                result=subprocess.run([sys.executable,'-B',*flags,str(HERE/(name+'.py'))],cwd=ROOT,
                                      capture_output=True,check=True,timeout=120)
                need(json.loads(result.stdout)==read(prefix+'-normal.json'),'replay '+name)
    counts=read('counting-normal.json')['ordinary_counts']
    external=[tuple(map(int,line.split())) for line in (HERE/'oeis-a000568-0-50.txt').read_text().splitlines() if line and not line.startswith('#')]
    need([n for n,v in external]==list(range(51)),'external universe')
    need([v for n,v in external]==counts,'OEIS comparison')
    old=ROOT/'05-knowledge/results/session67_metagraph_enumeration_information_20261026/session67-root/output/census.json'
    relocation=read('relocation.json')
    need(json.loads(old.read_text())==relocation,'inherited relocation receipt')
    for row in relocation:
        need(row['ordinary_classes']==counts[row['n']],'relocation versus Burnside')
        need(2*row['merged_classes']-row['self_converse']==row['ordinary_classes'],'converse lift')
    incoming=ROOT/'05-knowledge/results/relation_tower_frobenius_ports_20261028'
    root_switch=json.loads((incoming/'root_switch.json').read_text())
    for row in root_switch['results']:
        need(row['classes']==counts[row['order']],'incoming root-switch count')
        if row['order']<=8:
            need(row['tc_histogram']==relocation[row['order']-1]['tc_histogram'],'incoming whole presentation histogram')
    bank=json.loads(gzip.decompress((incoming/'tournament_bank.json.gz').read_bytes()))
    need(bank['order']==9 and len(bank['rows'])==counts[9],'incoming n9 bank size')
    need(len({tuple(row) for row in bank['rows']})==counts[9],'incoming no exact duplicate encodings')
    for rows in bank['rows']:
        need(len(rows)==9 and all(0<=x<512 for x in rows),'incoming row types')
        need(all(not(rows[i]>>i&1) for i in range(9)) and all(
            (rows[i]>>j&1)+(rows[j]>>i&1)==1 for i in range(9) for j in range(i+1,9)),
            'incoming tournament relation')
    if args.replay_metagraph:
        source=old.parent.parent/'enumerate.cpp'
        with tempfile.TemporaryDirectory(prefix='frobenius-metagraph-') as tmp:
            exe=Path(tmp)/('census.exe' if sys.platform=='win32' else 'census')
            subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-pedantic',str(source),'-o',str(exe)],
                           cwd=ROOT,check=True,capture_output=True,timeout=120)
            result=subprocess.run([str(exe),'8'],cwd=ROOT,check=True,capture_output=True,timeout=120)
            need(json.loads(result.stdout)==relocation,'compiled relocation replay')
    print(json.dumps(dict(sealed_files=len(names),source_pins='PASS',links='PASS',receipts='PASS',
                          external_counts='0..50 PASS',incoming_bank_types=counts[9],replay=args.replay,
                          metagraph_replay=args.replay_metagraph),sort_keys=True))


if __name__=='__main__':main()
