"""Pin imported repository sources, check local links, and hash the research packet."""
from pathlib import Path
import sys,importlib,hashlib,json,re
from urllib.parse import unquote
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from reproduce import JOBS
PROOFS=[
'01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md',
'01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md',
'05-knowledge/results/polarized_path_calculus_20261103/FUNNELS.md',
'05-knowledge/results/two_entry_funnels_20261103/READY.md',
'05-knowledge/results/ready_return_pentagon_20261103/OBSTRUCTION.md',
'05-knowledge/results/nonzero_preparation_20261103/surgery/ORDER_OBSTRUCTION.md',
'05-knowledge/results/berggren_target_return_20261111/reset/TARGET_RETURN.md',
'05-knowledge/results/cycle_opening_six_return_20261112/compiler/MODULE_RETURN.md',
'05-knowledge/results/cycle_opening_six_return_20261112/compiler/SIX_STATE_MINIMALITY.md',
'05-knowledge/results/session295_construction_atlas_20261113/bridges.md',
'05-knowledge/results/session296_multicore_metagraph_20261113/cores/PROOF.md',
'05-knowledge/results/session296_multicore_metagraph_20261113/cores/FINITE_ATLAS.md',
'05-knowledge/results/session296_multicore_metagraph_20261113/metagraph/PROOF.md',
'05-knowledge/results/mrb_seven_cutoffs_20261113/module/PROOF.md',
'05-knowledge/results/session298_exact_target_geometry_20261113/README.md',
'05-knowledge/results/session298_exact_target_geometry_20261113/TWO_WIDTH.md']

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run():
    for name,args in JOBS:importlib.import_module(name)
    paths={Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None)}
    paths={p for p in paths if p.is_relative_to(ROOT) and p.suffix=='.py'}
    paths.update(ROOT/p for p in PROOFS)
    pins=[dict(path=p.relative_to(ROOT).as_posix(),sha256=digest(p)) for p in sorted(paths)]
    (HERE/'SOURCE_PINS.json').write_text(json.dumps(dict(scope='Imported repository Python modules and named proof dependencies at seal time',files=pins),indent=2)+'\n',encoding='utf8')
    links=0
    for path in HERE.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text(encoding='utf8')):
            if target.startswith(('http:','https:','#')):continue
            local=unquote(target.split('#')[0])
            if not (path.parent/local).resolve().exists():raise RuntimeError('broken local link '+str(path)+': '+target)
            links+=1
    files=[dict(path=p.relative_to(HERE).as_posix(),sha256=digest(p),bytes=p.stat().st_size) for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.json']
    (HERE/'MANIFEST.json').write_text(json.dumps(dict(status='PASS',source_pins=len(pins),local_links=links,files=files),indent=2)+'\n',encoding='utf8')
    print('SEALED',len(pins),'source pins',links,'local links',len(files),'files')
if __name__=='__main__':run()
