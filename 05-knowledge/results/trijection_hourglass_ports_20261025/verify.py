"""Read-only publication checks; --refresh explicitly regenerates file pins."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MANIFEST=HERE/'manifest.json'
PAIRS=[('audit.normal.json','audit.optimized.json'),
       ('exception.normal.json','exception.optimized.json'),
       ('endpoint/normal.json','endpoint/optimized.json'),
       ('endpoint/tail-normal.json','endpoint/tail-optimized.json'),
       ('endpoint/general-normal.json','endpoint/general-optimized.json'),
       ('trijection/receipt.json','trijection/receipt-optimized.json'),
       ('codes/normal.json','codes/optimized.json'),
       ('incidence/normal.json','incidence/optimized.json'),
       ('incidence/common_normal.json','incidence/common_optimized.json')]
INHERITED=[
 '05-knowledge/results/ports_owners_f1_20261025/endpoint.md',
 '05-knowledge/results/endpoint_four_frontier_20261021/README.md',
 '05-knowledge/results/endpoint_three_order10_20261021/README.md',
 '05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md',
 '05-knowledge/results/reflection_bridge_20261025/geometric_carrier.md',
 '05-knowledge/results/section_defect_two_error_geometry_20261019/README.md',
 '05-knowledge/results/session32_square_monster_20261023/SYLOW_FIXED_POINT_PARTITION.md',
 '05-knowledge/results/session60_six_seven_geometry_20261025/seven-lines/RESULT.md',
 '05-knowledge/results/session60_six_seven_geometry_20261025/hourglass-symmetry/RESULT.md',
]


def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return path.read_text(encoding='utf-8-sig')


def main():
    files=sorted(p for p in HERE.rglob('*') if p.is_file() and
        '__pycache__' not in p.parts and p.suffix in {'.py','.md','.json'} and
        p.name not in {'manifest.json','verification.json','verification-optimized.json'})
    if '--refresh' in sys.argv:
        record={'packet':{str(p.relative_to(HERE)).replace('\\','/'):digest(p) for p in files},
                'inherited':{name:digest(ROOT/name) for name in INHERITED}}
        MANIFEST.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    record=json.loads(read(MANIFEST))
    need(set(record['packet'])=={str(p.relative_to(HERE)).replace('\\','/') for p in files},'manifest file set changed')
    for name,expected in record['packet'].items():
        need(digest(HERE/name)==expected,('packet hash mismatch',name))
    for name,expected in record['inherited'].items():
        need(digest(ROOT/name)==expected,('inherited hash mismatch',name))
    for a,b in PAIRS:
        need(json.loads(read(HERE/a))==json.loads(read(HERE/b)),('normal/optimized mismatch',a,b))
    local_links=0; parsed=0
    for path in files:
        if path.suffix=='.py':
            ast.parse(read(path),filename=str(path)); parsed+=1
        if path.suffix!='.md':
            continue
        for url in re.findall(r'\]\(([^)]+)\)',read(path)):
            if '://' in url or url.startswith('#'):
                continue
            target=url.split('#',1)[0].strip('<>')
            need((path.parent/target).exists(),('broken local link',path.name,url))
            local_links+=1
    print(json.dumps({'status':'PASS','packet_files':len(files),
        'inherited_files':len(INHERITED),'receipt_pairs':len(PAIRS),
        'parsed_python_files':parsed,'local_links':local_links},indent=2))


if __name__=='__main__':
    main()
