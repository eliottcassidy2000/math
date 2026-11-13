#!/usr/bin/env python3
"""Reproduce a bounded, source-pinned recovery catalog; not a new corpus scan."""
import hashlib
import json
import posixpath
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = '68d0d86bef480b6a36efc3d9864fd01ff3244bf5'
PREFIX = '05-knowledge/results/'
INDEX = PREFIX + 'natural_number_roles_20261004.md'
ATLAS = PREFIX + 'six_roots_hidden_bits_20261103/NUMBER_ATLAS.md'
MERGE = PREFIX + 'merge_18_19_11_22_collatz_hn_20261018.md'
CHAIN = PREFIX + 'joint_metagraph_dyadic_20261113/sporadic/BRIDGES.md'
CONIC = PREFIX + 'pythagorean_form_conic_graphs_rp_margin_20261018.md'
CARD = PREFIX + 'cardinality_bridges_20261018m.md'
ZERO = PREFIX + 'odd_tables_signed_carries_20261111/history/NATIVE_ZERO.md'
EIGHT = PREFIX + 'eight_nine_twentyseven_20261112/history/HISTORY.md'
ELEVEN = PREFIX + 'eleven_wedge_honeycomb_20261102/RECOVERY.md'
F64 = PREFIX + 'sextic_subfield_clock_20261004.md'
F8 = PREFIX + 'session54_reframed_numerology_20261025/sixtythree-bridge/README.md'
PALETTE = PREFIX + 'session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md'
RESPONSE = PREFIX + 'hn_metagraph_information_20261026/information/RESULT.md'
SOURCES = {
 INDEX: ('curated small-number index', '20261004 difference-families-oct04'),
 PREFIX+'constants_atlas_20260926_recurrent_numbers.md': ('lexical census; deliberately filters small integers', '20260926 collatz-exponent-atlas'),
 ATLAS: ('six-numeral lexical census and selective proof atlas', 'six_roots_hidden_bits_20261103'),
 MERGE: ('explicit owner-requested four-number survey', '20261018 sixcolour second brief'),
 CHAIN: ('audited typed bridges; not a numeral occurrence census', 'joint_metagraph_dyadic_20261113 sporadic'),
 CONIC: ('typed conic dictionary with finite controls', '20261018 pythagorean/form/conic'),
 CARD: ('targeted cardinality atlas; not an occurrence census', '20261018 frontier_petersen_oct18'),
 ZERO: ('native zero-digit construction; synthesized lookup for numeral zero', 'odd_tables_signed_carries_20261111'),
 EIGHT: ('targeted 8/9/27 recovery and new construction', 'eight_nine_twentyseven_20261112'),
 ELEVEN: ('targeted 20 and 7/15/24/26 recovery', 'eleven_wedge_honeycomb_20261102'),
 F64: ('proved finite-field norm/tensor clock', '20261004 sextic-subfield-clock'),
 F8: ('proved bridge conditional on independently reproduced finite blocker census', 'session54_reframed_numerology_20261025'),
 PALETTE: ('proved pair-palette, Fano-completion and holonomy transfer; finite controls', 'session61_hadwiger_six_seven_20261025'),
 RESPONSE: ('proved full labelled Hamilton-path flip-response reconstruction', 'hn_metagraph_information_20261026'),
}

def read(path):
    return subprocess.check_output(['git','show',f'{BASE}:{path}'], cwd=ROOT).decode()


def route(path, needle=None):
    s=read(path)
    line=next((i for i,v in enumerate(s.splitlines(),1) if needle and needle in v),1)
    return {'path':path, 'line':line, 'source_kind':SOURCES[path][0], 'session_route':SOURCES[path][1]}

rows = {}
for line_no, line in enumerate(read(INDEX).splitlines(),1):
    m=re.match(r'^\| \*\*(\d+)\*\* \| (.*?) \| (.*?) \|$',line)
    if not m or int(m[1])>32:
        continue
    n=int(m[1])
    rows[n]={
      'number':n, 'old_role_index_entry':True, 'recovery_status':'RECOVERED_INDEX_ENTRY',
      'claim_status':'INHERITED_SCOPED: exact roles and echoes must retain linked-source scope',
      'roles':m[2], 'relation_and_source_verbatim':m[3],
      'source':dict(route(INDEX), line=line_no),
      'linked_source_paths':[posixpath.normpath(str((Path(INDEX).parent / p).as_posix())) for p in re.findall(r'\]\(([^)]+)\)',m[3])],
      'request_provenance':'User reports prior request; individual archived request not independently located',
      'additional_routes':[]}

supplements={
 0:(ZERO,'native zero-digit return; a preserved two-birth interface multiplies H by3','PROVED in inherited source; route is not an old numeral census','## 2.'),
 22:(MERGE,'explicit four-number survey; 22 is a merge vertex congruent4 mod18, with separate HN/Jacobian roles','INHERITED_SCOPED; no cross-problem theorem inferred','| 22 |'),
 23:(ATLAS,'dedicated lexical census plus residue slots, square-sum order23 and Paley/Golay routes','INHERITED_SCOPED; lexical coverage is FINITE-EXACT in its stated source revision','## 23:'),
 25:(ATLAS,'dedicated lexical census; 25-coordinate score/Lorentzian route retains scores but loses H','INHERITED_SCOPED; score-hostile boundary retained','## 25:'),
 26:(CHAIN,'rank26 = Leech rank24 plus hyperbolic plane; signature(25,1)','INHERITED PROVED/CITED; 26 is total rank, not eta level','### The recovered 23/24/25/26'),
 28:(CONIC,'28 exterior conic points at p=7 carry the Coxeter graph','INHERITED FINITE-EXACT/classical identification; not a global number census','**p = 7**'),
 30:(CARD,'30 nonmatching Lagrangians in binary symplectic dimension6: labelled affine/Fano structures','INHERITED PROVED/FINITE-EXACT; inherited classical interpretations are CITED','At eight labels')}
for n,(path,meaning,status,needle) in supplements.items():
    rows[n]={'number':n,'old_role_index_entry':False,'recovery_status':'SUPPLEMENTAL_TYPED_ROUTE',
      'claim_status':status,'roles':meaning,'source':route(path,needle),
      'linked_source_paths':[], 'request_provenance':'User reports prior request; individual archived request not independently located',
      'additional_routes':[]}
for n in (1,23,25):
    rows[n]['additional_routes'].append(route(ATLAS,f'## {n}:'))
    rows[n]['request_provenance']='Inherited atlas explicitly identifies the annotation {1,23,25,63,65,67}; no individual archived chat opened'
for n in (11,18,19,22):
    rows[n]['additional_routes'].append(route(MERGE,f'| {n} |'))
    rows[n]['request_provenance']='Inherited note quotes owner request to recover 18,19,11,22; not an independently opened archived chat'
for n in (8,9,27): rows[n]['additional_routes'].append(route(EIGHT))
for n in (7,11,15,20,24,26): rows[n]['additional_routes'].append(route(ELEVEN))
for n in (23,24,25,26): rows[n]['additional_routes'].append(route(CHAIN,'### The recovered 23/24/25/26'))
for n in (4,6,8,14,15,16,21,28,30): rows[n]['additional_routes'].append(route(CARD))
for n in (2,3,4,6,7,8,9,21,32): rows[n]['additional_routes'].append(route(F64))
for n in (1,2,4,6,7,8,21): rows[n]['additional_routes'].append(route(F8))
for n in (2,3,6,7,8,15,30): rows[n]['additional_routes'].append(route(PALETTE))
for n in (1,2): rows[n]['additional_routes'].append(route(RESPONSE))
if set(rows)!=set(range(33)): raise RuntimeError('catalog must have exactly 0 through 32')
manifest=[]
for path,(kind,session) in SOURCES.items():
    data=read(path).encode()
    manifest.append({'path':path,'source_kind':kind,'session_route':session,'sha256':hashlib.sha256(data).hexdigest(),
                     'git_blob':subprocess.check_output(['git','rev-parse',f'{BASE}:{path}'],cwd=ROOT,text=True).strip()})
obj={'status':'SOURCE-PINNED RECOVERY CATALOG; not an exhaustive occurrence or request-history census',
     'source_commit':BASE,'universe':list(range(33)), 'literal_entries_in_old_role_index':26,
     'absent_from_old_role_index':sorted(supplements),'sources':manifest,
     'scope_limits':['bounded recovery of named atlases and targeted source notes',
       'no complete archived chat or deleted-only history scan',
       'a supplemental typed route does not establish a dedicated old census for that integer',
       'inherited mathematical statuses are not newly promoted by inclusion'],
     'rows':[rows[n] for n in range(33)]}
(HERE/'catalog.json').write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
# Compact readable table; source paths are relative to the present file.
lines=['# Recovered routes for every integer 0 through 32','',
 'Status: **source-pinned recovery index**, not an exhaustive numeral or request-history census.',
 'All inherited sources are pinned to `'+BASE+'` in [catalog.json](catalog.json).',
 'The original role index has26 rows in this range; † marks seven supplemental lookups.',
 'Every row preserves the source scope; inclusion does not prove a cross-domain identification.','',
 '| Integer | Source kind | Recovered role | Route |','|---:|---|---|---|']
for n in range(33):
    row=rows[n]; source=row['source']; p=source['path'].removeprefix(PREFIX)
    role=row['roles'].replace('|','\\|')
    kind='old role-index entry' if row['old_role_index_entry'] else source['source_kind']
    lines.append(f"| {n}{'' if row['old_role_index_entry'] else '†'} | {kind} | {role} | [{source['session_route']}](../../{p}) (line{source['line']}) |")
lines+=['','The JSON also preserves each original row verbatim, additional targeted routes, explicit',
 'owner-annotation evidence where found, source hashes, and the limits of request attribution.','']
(HERE/'TABLE.md').write_text('\n'.join(lines))
print(json.dumps({'rows':len(rows),'old_index_rows':26,'supplemental':sorted(supplements),'sources':len(manifest)},sort_keys=True))
