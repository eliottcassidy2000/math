"""Compile/replay the exact menu and actual-port producer, then consume them.

Default is read-only comparison. --write intentionally refreshes owned ledgers.
Requires Python 3 and g++ with C++17; all executables live in a temporary folder.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,shutil,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def need(x,why):
    if not x:raise RuntimeError(why)
def source_bytes(path):
    return path.read_bytes() if path.suffix=='.gz' else path.read_text(encoding='utf-8-sig').replace('\r\n','\n').encode()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    inherited=[ROOT/'04-computation/experiments/tiling_line_metagraph_census_20261020_classes_n8.tsv',ROOT/'05-knowledge/results/blue_quotient_20261020b/n9_classes.tsv.gz',ROOT/'05-knowledge/results/blue_quotient_20261020b/census.cpp',ROOT/'04-computation/experiments/sierpinski_ocf_20261022/drt23_classes.json']
    with tempfile.TemporaryDirectory(prefix='rooted_five_ports_') as td:
        tmp=Path(td);exe=tmp/('menu.exe' if sys.platform=='win32' else 'menu')
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'menu.cpp'),'-o',str(exe)],check=True)
        for n,path in zip((8,9),inherited[:2]):
            with (gzip.open(path,'rt') if n==9 else path.open()) as f:rows=list(csv.DictReader(f,delimiter='\t'))
            inp=tmp/f'n{n}.tsv';inp.write_text(''.join(f"{r['canon'] if n==8 else r['code']}\t{r['Aut']}\n" for r in rows))
            output=tmp/f'menu_n{n}.json';subprocess.run([str(exe),str(n),str(inp),str(output)],check=True)
            if args.write:shutil.copyfile(output,HERE/output.name)
            else:need(json.loads(output.read_text())==json.loads((HERE/output.name).read_text()),f'n{n} menu differs')
        subprocess.run([sys.executable,str(HERE/'probe.py')]+(['--write'] if args.write else []),check=True)
        probe=json.loads((HERE/'probe.json').read_text());targets=tmp/'targets.tsv';targets.write_text(''.join(f"{r['code']}\n" for r in probe['orders']['9']['records']))
        pexe=tmp/('ports.exe' if sys.platform=='win32' else 'ports')
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'ports.cpp'),'-o',str(pexe)],check=True)
        subprocess.run([str(pexe),str(tmp/'n9.tsv'),str(targets),str(tmp)],check=True)
        for name in ('focused_targets.tsv','port_producer.json'):
            if args.write:shutil.copyfile(tmp/name,HERE/name)
            else:need((tmp/name).read_text()==(HERE/name).read_text(),f'{name} differs')
        binary=(tmp/'focused_ports.bin').read_bytes()
        if args.write:(HERE/'focused_ports.bin.gz').write_bytes(gzip.compress(binary,compresslevel=9,mtime=0))
        else:need(gzip.decompress((HERE/'focused_ports.bin.gz').read_bytes())==binary,'literal port ledger differs')
    subprocess.run([sys.executable,str(HERE/'verify.py')]+(['--write'] if args.write else []),check=True)
    owned=('check.py','menu.cpp','ports.cpp','probe.py','verify.py','menu_n8.json','menu_n9.json','probe.json','focused_targets.tsv','focused_ports.bin.gz','port_producer.json','verification.json')
    manifest={'status':'FINITE-EXACT rooted-five incidence menu and complete focused n9 port fibres',
              'inherited_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(source_bytes(p)).hexdigest() for p in inherited},
              'owned_sha256':{name:hashlib.sha256(source_bytes(HERE/name)).hexdigest() for name in owned},
              'hash_convention':'gzip files hashed as committed bytes; text canonical LF UTF-8 without BOM',
              'producer':json.loads((HERE/'port_producer.json').read_text())}
    if args.write:(HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    else:need(json.loads((HERE/'manifest.json').read_text())==manifest,'manifest differs')
    print('rooted-five menu, literal ports, and independent consumer:', 'written' if args.write else 'PASS')
if __name__=='__main__':main()
