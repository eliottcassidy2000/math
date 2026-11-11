"""Full independent support-byte comparison, with fresh literal observers.

No inherited Python implementation is imported. Seven actual context counts
determine each raw-cut positive observer by the proved seven-coordinate law.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time

HERE=Path(__file__).resolve().parent
CORES=((0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11))
CONTEXTS=((0,0),(1,0),(3,0),(2,0),(0,1),(0,3),(1,1))

def need(p,msg):
    if not p:raise RuntimeError(msg)

def hamilton(rows):
    # First-vertex recurrence, independently implemented from the inherited
    # last-vertex path enumerator and residual-path signature extraction.
    n=len(rows);p=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for v in range(n):
            if not(mask>>v&1):continue
            rest=mask^(1<<v)
            if not rest:p[mask][v]=1
            else:p[mask][v]=sum(p[rest][w] for w in range(n) if rest>>w&1 and rows[v]>>w&1)
    return sum(p[-1]) if rows else 1

def child(core,cut,A,B):
    left=[int(b) for b in bin(A)[2:]] if A else []
    right=[1-int(b) for b in reversed(bin(B)[2:])] if B else []
    l=len(left);n=len(core);m=l+n+len(right)
    rows=[0]*(m+1)
    # Start from the ordinal-sum transitive edges; replace the core block.
    for u in range(m):
        for v in range(u+1,m):
            forward=not(l<=u<v<l+n) or bool(core[u-l]>>(v-l)&1)
            winner,loser=(u,v) if forward else (v,u)
            rows[winner]|=1<<loser
    out=left+[(cut>>v)&1 for v in range(n)]+right
    for v,bit in enumerate(out):rows[m if bit else v]|=1<<(v if bit else m)
    return tuple(rows)

def raw_observers():
    records=[]
    for ci,core in enumerate(CORES):
        for cut in range(1<<len(core)):
            f=[hamilton(child(core,cut,A,B)) for A,B in CONTEXTS]
            h,p10,p30,p20,p01,p03,p11=f
            need((p30-p10)%2==0 and (p03-p01)%2==0,'integral slope recovery')
            i=(p30-p10)//2;o=(p03-p01)//2
            s=p10-h-i;e=p01-h-o
            need((p20-h-2*i)%2==0,'integral core recovery')
            d=(p20-h-2*i)//2;L=p11-p10-p01+h
            need(d==hamilton(core) and d>0 and L>0,'recovered core controls')
            need(min(i,o,s,e)>=0,'nonnegative response')
            need((i-s)%2==0 and (o-e)%2==0 and L%2==0,'parity controls')
            records.append({'core_index':ci,'core_rows':core,'cut':cut,
                            'literal_context_counts':f,'positive_observer':[h,d,L,i,o,s,e]})
    need(len(records)==90,'raw-cut universe')
    return records

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=16777215)
    ap.add_argument('--fail',action='store_true');args=ap.parse_args()
    if args.fail:need(False,'intentional full-support audit failure')
    N=args.bound;need(N>0 and N%2==1,'positive odd target bound')
    records=raw_observers()
    observer_text='90\n'+''.join(' '.join(map(str,[r['core_index'],r['cut'],*r['positive_observer']]))+'\n' for r in records)
    started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='chanpak-full-support-') as td:
        tmp=Path(td);(tmp/'raw90.txt').write_text(observer_text)
        for stem in ('support','support_positive'):
            subprocess.run(['g++','-O3','-std=c++17',str(HERE/(stem+'.cpp')),'-o',str(tmp/(stem+'.exe'))],check=True,capture_output=True)
        old_json=tmp/'old.json';new_json=tmp/'positive.json'
        subprocess.run([str(tmp/'support.exe'),str(HERE/'signatures.txt'),str(N),str(old_json),str(tmp/'old.bin')],check=True,capture_output=True)
        subprocess.run([str(tmp/'support_positive.exe'),str(tmp/'raw90.txt'),str(N),str(new_json),str(tmp/'positive.bin')],check=True,capture_output=True)
        old=(tmp/'old.bin').read_bytes();new=(tmp/'positive.bin').read_bytes()
        need(old==new,'complete per-target per-core support bytes differ')
        need(len(old)==(N+1)//2,'mask byte count')
        old_report=json.loads(old_json.read_text());new_report=json.loads(new_json.read_text())
        need(old_report['support_mask_counts']==new_report['support_mask_counts'],'support distribution')
        retained=HERE/f'support_{N.bit_length()}.json'
        if retained.exists():need(old_report==json.loads(retained.read_text()),'original retained full receipt')
    # Keep constant-size raw observer provenance, never bulky support binaries.
    (HERE/'positive_raw_signatures.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
    return {'status':'PASS','bound':N,'raw_cut_count':len(records),
            'literal_observer_context_counts':len(records)*len(CONTEXTS),
            'support_bytes_compared':len(old),'all_support_bytes_equal':True,
            'support_sha256':hashlib.sha256(old).hexdigest(),
            'old_parameter_pairs':old_report['evaluated_parameter_pairs'],
            'positive_parameter_pairs':new_report['evaluated_parameter_pairs'],
            'old_signature_count':old_report['signature_count'],
            'direct_g_formula':'one plus binary-prefix/adjacent-one contributions at zero bits; no g-table recurrence',
            'observer_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'old_cpp_sha256':hashlib.sha256((HERE/'support.cpp').read_bytes()).hexdigest(),
            'positive_cpp_sha256':hashlib.sha256((HERE/'support_positive.cpp').read_bytes()).hexdigest(),
            'raw_observer_sha256':hashlib.sha256(observer_text.encode()).hexdigest(),
            'scope':'Complete value-support comparison at this bound: old signed44-row A-first inverse versus independent literal90-row positive B-first inverse. All-order coverage remains open.'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
