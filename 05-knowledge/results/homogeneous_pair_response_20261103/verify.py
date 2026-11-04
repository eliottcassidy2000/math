"""Frozen packet verifier; replay writes only to a temporary directory."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path, normalize=False):
    data = path.read_bytes()
    if normalize:
        data = data.replace(b'\r\n',b'\n')
    return sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    args = parser.parse_args()
    manifest = json.loads((HERE/'manifest.json').read_text(encoding='utf8'))
    inventory = {p.name for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json'}
    need(inventory==set(manifest['files']), 'inventory')
    for name,value in manifest['files'].items():
        need(digest(HERE/name)==value,'packet hash: '+name)
    for name,value in manifest['dependencies_lf_sha256'].items():
        need(digest(ROOT/name,True)==value,'dependency hash: '+name)
    for stem,script in [('check','check.py'),('response','response.py')]:
        receipt = HERE/(stem+'.json')
        need(receipt.read_bytes()==(HERE/(stem+'.optimized.json')).read_bytes(),'optimized: '+stem)
        pins = json.loads(receipt.read_text(encoding='utf8'))['source_sha256']
        if isinstance(pins,str):
            pins = {script:pins}
        for name,value in pins.items():
            need(digest(HERE/name)==value,'source pin: '+name)
    need((HERE/'response_tables.json').read_bytes()==
         (HERE/'response_tables.optimized.json').read_bytes(),'complete tables normal/optimized')
    response = json.loads((HERE/'response.json').read_text(encoding='utf8'))
    need(digest(HERE/'response_tables.json')==response['tables_sha256'],'table receipt')
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='pair-replay-') as folder:
            for optimized in (False,True):
                prefix = [sys.executable]+(['-O'] if optimized else [])
                for script,stem in [('check.py','check'),('response.py','response')]:
                    out = Path(folder)/(stem+'.json')
                    command = prefix+[str(HERE/script),'--output',str(out)]
                    if stem=='response':
                        tables = Path(folder)/'tables.json'
                        command += ['--tables',str(tables)]
                    subprocess.run(command,check=True)
                    need(out.read_bytes()==(HERE/(stem+'.json')).read_bytes(),'replay: '+stem)
                    if stem=='response':
                        need(tables.read_bytes()==(HERE/'response_tables.json').read_bytes(),'table replay')
    print(json.dumps(dict(status='PASS',files=len(inventory),
                          dependencies=len(manifest['dependencies_lf_sha256']),
                          rooted_expansions=5405,compiled_queries=42668,full_named_cores=6,
                          receipt_pairs=2,replayed=args.replay,
                          scope='Pair branch excluded; prime branch and universal reconstruction OPEN'),
                     indent=2,sort_keys=True))


if __name__=='__main__':
    main()

