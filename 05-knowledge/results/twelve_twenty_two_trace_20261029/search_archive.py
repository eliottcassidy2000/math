"""Rebuild a mechanical 12/22 archive index; this is not a semantic audit.

Run explicitly with --write to refresh the search snapshot. Mathematical
verifiers do not rescan the advancing repository.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ROOTS = ['01-canon/theorems','05-knowledge/results','07-reflections']
PATTERN = r'(?<!\d)(12|22)(?!\d)'


def read(path):
    absolute = str(ROOT/path)
    if os.name == 'nt':
        absolute = '\\\\?\\'+absolute
    return Path(absolute).read_text(encoding='utf-8',errors='replace')


def run():
    result = subprocess.run(['rg','-l','--pcre2',PATTERN,*ROOTS,'-g','*.md'],
                            cwd=ROOT,capture_output=True,text=True,encoding='utf-8')
    if result.returncode not in (0,1):
        raise RuntimeError(result.stderr)
    files = sorted(path.replace('\\','/') for path in result.stdout.splitlines()
                   if 'twelve_twenty_two_trace_20261029' not in path)
    both, counts, failures = [], Counter(), []
    for path in files:
        try:
            text = read(path)
        except OSError as error:
            failures.append(dict(path=path,error=str(error)))
            continue
        counts[path.split('/')[0]] += 1
        if re.search(r'(?<!\d)12(?!\d)',text) and re.search(r'(?<!\d)22(?!\d)',text):
            both.append(path)
    return dict(cached_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                roots=ROOTS,pattern=PATTERN,excluded_packet=HERE.name,
                candidate_count=len(files),readable_by_root=dict(counts),both_count=len(both),
                both_paths=both,read_failures=failures,
                scope='Mechanical digit-block hits, including references and equation numbers; not claims of mathematical relevance or exhaustive reading.')


if __name__ == '__main__':
    if '--write' not in sys.argv:
        raise SystemExit('Use --write to explicitly refresh this archive-search snapshot.')
    data=run()
    (HERE/'archive_search.json').write_text(json.dumps(data,separators=(',',':'))+'\n',encoding='utf-8')
    print('Indexed',data['candidate_count'],'candidate documents;',data['both_count'],'contain both; unreadable',len(data['read_failures']))
