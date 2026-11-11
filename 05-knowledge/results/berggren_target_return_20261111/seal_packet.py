"""Seal/check research artifacts, inherited sources and local Markdown links."""
from pathlib import Path
from urllib.parse import unquote
import argparse, hashlib, json, re

BASE=Path(__file__).resolve().parent
REPO=BASE.parents[2]
RESULTS=REPO/'05-knowledge/results'
OUTPUTS={'MANIFEST.json','SOURCE_PINS.json','LINK_CHECK.json'}
CODE_TREES=[
 '05-knowledge/results/uniform_return_metagraph_20261111/reset',
 '05-knowledge/results/chanpak_five_core_bridges_20261110/alternatives',
 '05-knowledge/results/polarized_path_calculus_20261103',
 '05-knowledge/results/two_entry_funnels_20261103',
 '04-computation/nonzero_preparation_20261103',
]
EXPLICIT=[
 '05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py',
 '05-knowledge/results/uniform_return_metagraph_20261111/arithmetic/tail_squares.json',
 '05-knowledge/results/uniform_return_metagraph_20261111/reset/ATLAS_TWO_BIRTH.md',
 '05-knowledge/results/uniform_return_metagraph_20261111/reset/READY_REBASE.md',
 '05-knowledge/results/uniform_return_metagraph_20261111/arithmetic/DYADIC_DEFECT.md',
 '05-knowledge/results/session256_chan_pak_triples_20261111/spinor/PROOF.md',
]

def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def formatted(value):
    return json.dumps(value,indent=2,sort_keys=True)+'\n'

def collect():
    records=[];inherited=set()
    for file in sorted(BASE.rglob('*.md')):
        for raw in re.findall(r'(?<!!)\[[^\]]+\]\(([^\n]+?)\)',file.read_text(encoding='utf-8')):
            target=raw.strip().strip('<>')
            if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',target) or target.startswith('#'):
                continue
            target=unquote(target.split('#')[0])
            path=(file.parent/target).resolve()
            if not path.exists():
                raise RuntimeError('Missing local link: '+str(file.relative_to(BASE))+' -> '+target)
            if not path.is_relative_to(REPO):
                raise RuntimeError('Local reference escapes repository: '+str(path))
            records.append(dict(source=file.relative_to(BASE).as_posix(),
                                target=path.relative_to(REPO).as_posix()))
            if path.is_file() and not path.is_relative_to(BASE):inherited.add(path)
    for name in CODE_TREES:
        folder=REPO/name
        if not folder.is_dir():raise RuntimeError('Missing dependency tree '+name)
        inherited.update(p for p in folder.rglob('*') if p.is_file() and p.suffix in{'.py','.cpp','.h'})
    for name in EXPLICIT:
        path=REPO/name
        if not path.is_file():raise RuntimeError('Missing explicit dependency '+name)
        inherited.add(path)
    pins=dict(hash_basis='SHA256 of LF-normalized bytes',
              scope='All external local Markdown references, explicit data dependencies and source code in declared inherited compiler trees; not a whole-repository proof audit.',
              sources={p.relative_to(REPO).as_posix():sha(p) for p in sorted(inherited)})
    links=dict(status='PASS',scope='Every relative Markdown file link in this packet; anchors and external URL availability are not checked.',
               checked=len(records),links=records)
    artifacts={p.relative_to(BASE).as_posix():sha(p) for p in sorted(BASE.rglob('*'))
               if p.is_file() and p.name not in OUTPUTS and '__pycache__' not in p.parts}
    manifest=dict(hash_basis='SHA256 of LF-normalized bytes',artifacts=artifacts,
                  source_pins_sha256=hashlib.sha256(formatted(pins).encode()).hexdigest(),
                  link_check_sha256=hashlib.sha256(formatted(links).encode()).hexdigest())
    return {'SOURCE_PINS.json':pins,'LINK_CHECK.json':links,'MANIFEST.json':manifest}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
    reports=collect()
    for name,value in reports.items():
        target=BASE/name
        if args.check:
            if json.loads(target.read_text(encoding='utf-8'))!=value:
                raise RuntimeError('Publication/dependency mismatch: '+name)
        else:target.write_text(formatted(value),encoding='utf-8')
    print('PASS:',len(reports['MANIFEST.json']['artifacts']),'artifacts;',
          len(reports['SOURCE_PINS.json']['sources']),'inherited sources;',
          reports['LINK_CHECK.json']['checked'],'local links')

if __name__=='__main__':main()
