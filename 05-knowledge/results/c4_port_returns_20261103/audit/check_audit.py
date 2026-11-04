"""Read-only independent whole-attachment audit; writes only this directory."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def main():
    with tempfile.TemporaryDirectory(prefix='c4-independent-attachment-') as directory:
        exe=Path(directory)/'audit.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'whole_attachment.cpp'),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],check=True,capture_output=True,text=True)
    need(not result.stderr,'unexpected compiler audit stderr')
    pattern=r'order=(\d+) graphs=(\d+) all_portsets=(\d+) strong_portsets=(\d+) f63=(\d+) hits22_24_42_44=(\d+),(\d+),(\d+),(\d+),'
    found=[list(map(int,m)) for m in re.findall(pattern,result.stdout)]
    expected=[[1,1,2,1,0,0,0,0,0], [2,1,4,2,0,0,0,0,0],
              [3,2,16,11,0,0,0,0,0], [4,4,64,45,0,0,0,0,0],
              [5,12,384,308,10,0,0,0,0], [6,56,3584,3133,22,3,0,3,0],
              [7,3584,458752,426587,153,7,0,0,0]]
    need(found==expected,'independent audit universe or result changed')
    source_names=['regeneration/small_bank.cpp','regeneration/full_small.json',
                  'regeneration/complete7_small.json','native/PROOF.md',
                  'native/check_native.py','prefix/PROOF.md','prefix/probe.cpp','prefix/check.py']
    pins={name:hashlib.sha256((PACKET/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for name in source_names}
    records=[dict(order=row[0],graphs=row[1],all_portsets=row[2],strong_portsets=row[3],
                  fixed_f63_candidates=row[4],hits_by_digit_word=dict(zip(('22','24','42','44'),row[5:])))
             for row in found]
    report=dict(status='INDEPENDENT_FINITE_EXACT_WHOLE_ATTACHMENT_AUDIT',
        class_generation='All labeled codes through order6, marking their complete permutation orbits; all56 representatives and64 birth cuts at order7.',
        order7_scope='3584 extensions cover every isomorphism class but retain duplicates; candidate and hit multiplicities are not distinct rooted-class counts.',
        cover_method='Enumerate all vertex permutations and all valid two-slot splits; no endpoint convolution.',
        strongness_method='Construct each actual full clock1 attachment; forward and reverse reachability.',
        independent_full_counts='For every one of185 fixed-f63 candidates, directly count full Hamilton paths at clocks1 and2; extract slope and intercept.',
        all_order_dependency=dict(status='CITED',source='Arthur H. Busch, EJC13(2006) N3, Theorem1/Corollary1, pp3-4',
            url='https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/',
            statement='Strong tournaments of order at least10 have at least125 Hamilton paths; a fixed-chart clock1 target91/97 therefore has suffix order at most7.'),
        mathematical_review='REVIEW.md separately audits enumeration completeness, native padded guards, and fixed-prefix coefficient logic.',
        source_pins=pins,orders=records,total_fixed_f63_candidates=185,
        total_hits=dict(zip(('22','24','42','44'),(10,0,6,0))),
        scope='No fixed-suffix strong common-clock attachment has(f,c)=(63,-35)or(63,-29); not an unrestricted Hamilton-count or shifted-clock exclusion.',
        reproduction=['python 05-knowledge/results/c4_port_returns_20261103/audit/check_audit.py',
                      'python -O 05-knowledge/results/c4_port_returns_20261103/audit/check_audit.py'])
    (HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(result.stdout.strip())


if __name__=='__main__':
    main()
