"""Pin external inputs, check prose links, and seal the complete research packet.

Run after all producer files and both replays are stable. --verify is read-only.
The manifest excludes only itself and Python bytecode caches. Text uses
LF-normalized SHA256; gzip and other binary inputs use their actual bytes.
"""
from pathlib import Path
import argparse
import json
import re
from urllib.parse import unquote
import reproduce_all as replay

HERE=replay.HERE
ROOT=replay.ROOT
native=replay.native
digest=replay.digest
read_json=replay.read_json
write_json=replay.write_json
require=replay.require

# Explicit entry points supplement runtime capture. Every recursive repository
# import and data read is collected separately by reproduce_all's audited child.
EXPLICIT_INPUTS=[
 '05-knowledge/results/metagraph_aspect_sequences_20261111/records.jsonl.gz',
 '05-knowledge/results/chanpak_five_core_bridges_20261110/alternatives/check.py',
 '05-knowledge/results/berggren_target_return_20261111/reset/funnel_six_state.py',
]

def prose_only(text):
    lines=[];fence=None
    for line in text.splitlines():
        marker=re.match(r'^\s*(`{3,}|~{3,})',line)
        if fence is not None:
            if marker and marker.group(1)[0]==fence[0] and len(marker.group(1))>=fence[1]: fence=None
            continue
        if marker:
            fence=(marker.group(1)[0],len(marker.group(1)));continue
        lines.append(line)
    return re.sub(r'(`+)(.*?)\1','', '\n'.join(lines), flags=re.S)

def inventory():
    return [p for p in sorted(HERE.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo') and p.name!='MANIFEST.json']

def check_links():
    linked=set();links=[]
    for note in sorted(HERE.rglob('*.md')):
        text=native(note).read_text(encoding='utf-8-sig')
        for target in re.findall(r'\]\((<[^>]+>|[^)\n]+)\)',prose_only(text)):
            target=target.strip()
            if target.startswith('<'):target=target[1:target.index('>')]
            else:target=re.sub(r'\s+[\"\'][^\"\']*[\"\']\s*$','',target)
            if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:',target) or target.startswith('#'):continue
            target=unquote(target.split('#',1)[0])
            if not target:continue
            path=(note.parent/target).resolve()
            require(native(path).exists(),f'Broken prose link {note.relative_to(HERE)} -> {target}')
            links.append({'note':note.relative_to(HERE).as_posix(),'target':target})
            if native(path).is_file() and not path.is_relative_to(HERE):
                require(path.is_relative_to(ROOT),'External local link lies outside the shared repository: '+str(path))
                linked.add(path.relative_to(ROOT).as_posix())
    return linked,{'status':'PASS','scope':'Explicit local Markdown prose links; inline and fenced code removed; URL fragments excluded; URLs not fetched','local_links':len(links),'links':links}

def collect_sources(linked):
    sources=set(linked)|set(EXPLICIT_INPUTS)
    reasons={name:['explicit entry point'] for name in EXPLICIT_INPUTS}
    for name in linked:reasons.setdefault(name,[]).append('linked external prose source')
    for capture in sorted((HERE/'_replay').rglob('*.inputs.json')):
        data=read_json(capture)
        for name,expected in data['external_repository_inputs'].items():
            require(digest(ROOT/name)==expected,'Captured source changed since engine replay: '+name)
            sources.add(name)
            reasons.setdefault(name,[]).append('runtime '+data['script'])
    for name in sources:require(native(ROOT/name).is_file(),'Missing pinned input '+name)
    return {name:digest(ROOT/name) for name in sorted(sources)}, {name:sorted(set(reasons.get(name,[]))) for name in sorted(sources)}

def verify_replays():
    declared={e['script'] for e in replay.ENGINES}|{'reproduce_all.py','seal_packet.py'}
    present={p.relative_to(HERE).as_posix() for p in HERE.rglob('*.py') if '__pycache__' not in p.parts}
    require(present<=declared,'Unregistered packet Python files: '+str(sorted(present-declared)))
    ids=[e['id'] for e in replay.ENGINES]
    reports={mode:read_json(HERE/('REPLAY_'+mode.upper()+'.json')) for mode in ('normal','optimized')}
    tables={}
    for mode,data in reports.items():
        require(data.get('complete') is True,'Incomplete '+mode+' replay')
        require(data.get('declared_engine_ids')==ids,'Outdated '+mode+' replay inventory')
        rows={row['id']:row for row in data['engines']}
        require(set(rows)==set(ids),'Missing/extra '+mode+' engine records')
        tables[mode]=rows
        for spec in replay.ENGINES:
            row=rows[spec['id']]
            require(row['exit_code']==0,'Failed '+mode+' engine')
            require(row['script_lf_sha256']==digest(HERE/spec['script']),'Engine source changed after '+mode+' replay: '+spec['script'])
            for item in row['outputs']:
                require(item['sha256']==digest(HERE/item['receipt']),'Receipt changed after '+mode+' replay: '+item['receipt'])
    for name in ids:
        first=tables['normal'][name];second=tables['optimized'][name]
        require(first['outputs']==second['outputs'],'Normal/optimized output mismatch: '+name)
        a=read_json(HERE/first['input_capture']);b=read_json(HERE/second['input_capture'])
        require(a['external_repository_inputs']==b['external_repository_inputs'],'Normal/optimized source-input mismatch: '+name)
    return {'status':'PASS','engines':len(ids),'normal_optimized_receipts_match':True,'hash_basis':'LF-normalized text bytes; binary bytes unchanged'}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify','--check',dest='verify',action='store_true',help='Read-only verification of the existing manifest and pins')
    ap.add_argument('--links-only',action='store_true',help='Check current prose links without writing a seal')
    args=ap.parse_args()
    if args.verify:
        pins=read_json(HERE/'SOURCE_PINS.json')
        for name,expected in pins.items():require(digest(ROOT/name)==expected,'External source pin differs: '+name)
        expected=read_json(HERE/'MANIFEST.json')
        current={p.relative_to(HERE).as_posix():digest(p) for p in inventory()}
        require(current==expected,'Packet manifest differs; inspect files and reseal intentionally')
        verify_replays()
        check_links()
        print(f'VERIFIED {len(current)} packet files and {len(pins)} source pins')
        return
    for path in inventory():
        if path.suffix.lower() in replay.TEXT_SUFFIXES:
            native(path).read_text(encoding='utf-8-sig')
    linked,link_report=check_links()
    if args.links_only:
        print(json.dumps({'status':'PASS','local_links':link_report['local_links'],'linked_external_files':len(linked)}));return
    replay_report=verify_replays()
    pins,reasons=collect_sources(linked)
    write_json(HERE/'SOURCE_PINS.json',pins)
    write_json(HERE/'SOURCE_PIN_DETAILS.json',{'status':'PASS; complete normal and optimized replay inputs verified','hash_basis':'SHA256 of LF-normalized text; binary gzip untouched','capture':'Recursive repository imports and audited opens, plus exact external prose links','reasons':reasons})
    write_json(HERE/'LINK_CHECK.json',link_report)
    write_json(HERE/'REPLAY_COMPARISON.json',replay_report)
    manifest={p.relative_to(HERE).as_posix():digest(p) for p in inventory()}
    write_json(HERE/'MANIFEST.json',manifest)
    print(f'SEALED {len(manifest)} files, {len(pins)} external source pins, {link_report["local_links"]} prose links')

if __name__=='__main__':main()
