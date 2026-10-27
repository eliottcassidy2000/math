"""Verify sealed exact receipts; --replay repeats both Python and C++ paths."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def need(ok,message):
    if not ok: raise RuntimeError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def read(name):
    return json.loads((HERE/name).read_text(encoding='utf-8-sig'))

def run_json(args):
    result=subprocess.run(args,cwd=ROOT,check=True,capture_output=True)
    return json.loads(result.stdout.decode('utf-8-sig'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args()
    seals={}
    for line in (HERE/'SHA256SUMS').read_text().splitlines():
        h,name=line.split('  ',1)
        need(name not in seals and Path(name).name==name,'manifest path')
        seals[name]=h
    actual={p.name for p in HERE.iterdir() if p.is_file() and p.name!='SHA256SUMS'}
    need(actual==set(seals),'manifest coverage')
    for name,h in seals.items():need(digest(HERE/name)==h,'hash: '+name)
    for name in ['normal','hostile-audit-normal']:
        need(read(name+'.json')==read(name.replace('normal','optimized')+'.json'),'normal/-O: '+name)
    normal=read('normal.json');ranks=read('rank.json')
    need(ranks['modulus']==1000000007,'rank prime')
    for a,b in zip(normal['tournaments']['atlases'],ranks['atlases']):
        need(a['n']==b['n'] and a['tournament_types']==b['types'],'independent universe')
        need(a['counted_response_rank_from_constant']==b['rank_mod_prime']==a['converse_orbits'],'independent rank')
    r=ranks['atlases'][-1]
    need(r=={'n':7,'labelled_tournaments':2097152,'types':456,'converse_orbits':272,'rank_mod_prime':272,'last_pivot_power':271},'n7 certificate')
    hostile=read('hostile.json');audit=read('hostile-audit-normal.json')
    need(hostile[0]['L2_one']-hostile[1]['L2_one']==audit['L2_one_difference']==29097684,'hostile difference')
    for r in hostile:need(r['Aut']*r['sum_partner_H_over_presentations']==r['L2_one'],'path/presentation factor')
    for path in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
            if '://' not in target and not target.startswith('#'):
                need((path.parent/target.split('#')[0]).resolve().exists(),'local link: '+target)
    sources=read('SOURCES.json')
    need(sources['executable_inputs']==[],'standalone executable scope')
    print('Packet seals, links, normal/-O receipts and independent rank overlap agree.',flush=True)
    if args.replay:
        for script,receipt in [('check.py','normal.json'),('audit_hostile.py','hostile-audit-normal.json')]:
            for optimized in [False,True]:
                command=[sys.executable,'-B']+(['-O'] if optimized else [])+[str(HERE/script)]
                need(run_json(command)==read(receipt),'replay: '+script)
                print(script+(' optimized' if optimized else ' normal')+': reproduced',flush=True)
        compiler=shutil.which('g++');need(compiler,'g++ required for --replay')
        with tempfile.TemporaryDirectory(prefix='clock-even-') as directory:
            for stem in ['rank','hostile']:
                exe=str(Path(directory)/(stem+('.exe' if sys.platform=='win32' else '')))
                subprocess.run([compiler,'-std=c++17','-O3',str(HERE/(stem+'.cpp')),'-o',exe],check=True)
                need(run_json([exe])==read(stem+'.json'),'C++ replay: '+stem)
                print(stem+'.cpp: reproduced',flush=True)
    print('Finite claims verified; all-size conjecture and attachment minimality remain OPEN.')

if __name__=='__main__':main()
