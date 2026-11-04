#!/usr/bin/env python3
"""Exact order-eight conference characterization; no third-party packages."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib,json,struct,subprocess
ROOT=Path(__file__).resolve().parent
BANK=ROOT.parent/'eight_pair_quotient'/'blind8.bin'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def rows(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return a

def seidel(a):return [[0 if i==j else 1 if a[i]>>j&1 else -1 for j in range(len(a))] for i in range(len(a))]

def main():
    files=['PROOF.md','verify.py','census.cpp'];pins={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files};raw=BANK.read_bytes();bank_pin=hashlib.sha256(raw).hexdigest();need(bank_pin=='2fb59e122ab3625ba1c5a3dd311ff3abe39ef7ceb1390083a6fda4a7367457ba','complete reviewed h-zero dependency')
    with TemporaryDirectory() as td:
        outputs=[]
        for optimization in ['-O2','-O0']:
            exe=Path(td)/('census'+optimization);subprocess.run(['g++','-std=c++17',optimization,'-Wall','-Wextra','-Werror',str(ROOT/'census.cpp'),'-o',str(exe)],check=True);outputs.append(subprocess.check_output([str(exe)],cwd=ROOT))
        need(outputs[0]==outputs[1],'optimized and unoptimized C++ exact replay');(ROOT/'census.json').write_bytes(outputs[0])
    data=json.loads(outputs[0]);base=data['normalized_seven_codes'];need(len(base)==240 and base==sorted(set(base)),'all normalized codes retained')
    for c in base:
        a=rows(7,c);S=seidel(a)
        need(all(sum(r)==0 for r in S),'normalized row sums')
        need(all(sum(S[i][k]*S[j][k] for k in range(7))==(6 if i==j else -1) for i in range(7) for j in range(7)),'all normalized integer Gram matrices')
    pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
    masks=[sum((((switches>>(i-1))&1 if i else 0)^((switches>>(j-1))&1))<<k for k,(i,j) in enumerate(pairs)) for switches in range(128)]
    generated=sorted(((c<<7)|127)^m for c in base for m in masks);expected=list(struct.unpack('<30720I',raw))
    need(len(set(generated))==30720 and generated==expected,'independent full switching set equality')
    control=seidel(rows(4,63));bad=sum(control[0][k]*control[1][k] for k in range(4));need(bad==2,'order-four nonconference boundary')
    need(pins=={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files} and BANK.read_bytes()==raw,'frozen source/dependency')
    return dict(status='PASS',all_seven_tournaments=2097152,regular_seven_tournaments=data['regular_seven_tournaments'],normalized_blocks=240,switchings_per_normalized_block=128,labelled_eight_conference_tournaments=30720,complete_h_zero_code_equality=True,independent_integer_gram_entries=240*49,nonconference_order_four_row_inner_product=bad,source_sha256=pins,dependency_sha256={'eight_pair_quotient/blind8.bin':bank_pin},data_sha256={'census.json':hashlib.sha256(outputs[0]).hexdigest()},scope='Order-eight equivalence of the complete h-zero kernel and skew-conference tournaments, not a general-order characterization.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
