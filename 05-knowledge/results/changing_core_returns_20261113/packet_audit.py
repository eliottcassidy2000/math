"""Audit packet links, inherited source pins, and normal/optimized receipts."""
from pathlib import Path
import hashlib,json,re,subprocess
from urllib.parse import unquote
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(data):return hashlib.sha256(data).hexdigest()
def raw(path):return path.read_bytes()
def normalized(path):return raw(path).replace(b'\r\n',b'\n')
def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items() if k not in ('python_optimized','optimized')}
    if isinstance(value,list):return [clean(x) for x in value]
    return value
links=0
for path in HERE.rglob('*.md'):
    text=path.read_text(encoding='utf-8');need('\ufffd' not in text,('replacement character',str(path)))
    for target in re.findall(r'\[[^\]\n]*\]\(([^()\n]+)\)',text):
        if target.startswith(('https:','http:','mailto:','#')):continue
        file=unquote(target.split('#')[0]).strip('<>')
        # Algebra such as [r=0](p-h) is not a document link.
        if not file or ('.' not in file and '/' not in file):continue
        resolved=(path.parent/file).resolve();need(resolved.exists(),('broken local link',str(path),target));links+=1
pairs=[]
for opt in HERE.rglob('*.optimized.json'):
    normal=opt.with_name(opt.name.replace('.optimized.json','.json'))
    need(normal.exists(),('missing normal receipt',str(opt)))
    need(clean(json.loads(normal.read_text()))==clean(json.loads(opt.read_text())),('receipt mismatch',str(opt)))
    pairs.append(str(normal.relative_to(HERE)))
pins=0
geo=json.loads((HERE/'geometry/SOURCE_PINS.json').read_text())
for item in geo['sources']:
    data=subprocess.check_output(['git','show',geo['source_commit']+':'+item['path']],cwd=ROOT)
    need(digest(data)==item['sha256'],('historical source pin',item['path']));pins+=1
reviewed=0
for key,items in geo.items():
    if not key.startswith('audited_'):continue
    for item in items:
        if not isinstance(item,dict) or 'path' not in item or 'sha256' not in item:continue
        path=ROOT/item['path']
        need(item['sha256'] in (digest(raw(path)),digest(normalized(path))),('final reviewed source pin',item['path']))
        reviewed+=1
prep=json.loads((HERE/'preparation/source_manifest.json').read_text())
for file,sha in prep['source_sha256'].items():
    path=ROOT/file
    need(sha in (digest(raw(path)),digest(normalized(path))),('preparation current dependency pin',file));pins+=1
for name,item in prep['receipt_pairs'].items():
    path=HERE/'preparation'/(name+'.json')
    need(item['sha256'] in (digest(raw(path)),digest(normalized(path))),('preparation receipt hash',name))
for file,sha in json.loads((HERE/'returns/source_pins.json').read_text()).items():
    path=ROOT/file;need(sha in (digest(raw(path)),digest(normalized(path))),('return source pin',file));pins+=1
# The source census is independently pinned inside the finite-head receipt.
receipt=json.loads((HERE/'returns/birth_budget_receipt.json').read_text())
census=ROOT/receipt['source']['path']
need(receipt['source']['sha256'] in (digest(raw(census)),digest(normalized(census))),'inherited census hash');pins+=1
report=dict(status='PASS',local_markdown_links=links,normal_optimized_pairs=pairs,inherited_pins_checked=pins,final_reviewed_sources=reviewed,
    scope='Does not regenerate inherited isomorphism census or replace the mathematical proof reviews; every optimized receipt is compared with its normal partner.')
(HERE/'packet_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(report,indent=2))
