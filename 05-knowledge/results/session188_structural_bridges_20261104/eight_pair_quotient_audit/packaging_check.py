#!/usr/bin/env python3
"""Independent two-bank lossless packaging and prevalidation review."""
from pathlib import Path
import base64,gzip,hashlib,json,os,shutil,struct,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent.parent/'eight_pair_quotient'

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    meta=json.loads((ROOT/'REASSEMBLY.json').read_text());need(len(meta['groups'])==2,'two complete bank groups');parts=[];banks={};streams=[]
    for group in meta['groups']:
        need(len(group['outputs'])==1,'one output per bank');chunks=[]
        for item in group['parts']:
            raw=(ROOT/item['file']).read_bytes();need(len(raw)==item['size'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'every chunk size/hash');parts.append(item['file']);chunks.append(b''.join(raw.split()))
        compressed=base64.b64decode(b''.join(chunks),validate=True);raw=gzip.decompress(compressed);item=group['outputs'][0]
        need(raw==(ROOT/item['file']).read_bytes(),'every original bank byte equal')
        need(len(raw)==item['size'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'complete output integrity');banks[item['file']]=raw
        streams.append(dict(file=item['file'],raw_bytes=len(raw),raw_sha256=hashlib.sha256(raw).hexdigest(),gzip_bytes=len(compressed),gzip_sha256=hashlib.sha256(compressed).hexdigest(),parts=len(group['parts'])))
    table=banks['t7.bin'];raw=banks['blind8.bin'];need(len(table)==1<<21 and all(v<128 and v.bit_count()%2==0 for v in table),'all seven-bit vector coordinate conventions');need(table.count(0)==50688,'full zero vector count')
    codes=list(struct.unpack('<30720I',raw));need(codes==sorted(set(codes)) and all(c<1<<28 for c in codes),'all unsigned little-endian labelled words')
    pairs=[(i,j) for i in range(8) for j in range(i+1,8)]
    keep=[[k for k,(i,j) in enumerate(pairs) if deleted not in (i,j)] for deleted in range(8)]
    for c in codes:
        for positions in keep:
            sub=sum(((c>>old)&1)<<new for new,old in enumerate(positions));need(table[sub]==0,'all physically named deletion indices')
    scenarios=['roundtrip','corrupt_first','corrupt_late','conflict_first','conflict_late','correct_first_conflicting_late','bad_late_hash','nonlocal_late_output','duplicate_output','nonlocal_late_chunk']
    controls=[];files=['REASSEMBLY.md','REASSEMBLY.json','reassemble.py']+parts;pins={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}
    for optimized in (False,True):
        for scenario in scenarios:
            with tempfile.TemporaryDirectory(prefix='eight_pair_package_') as td:
                dst=Path(td)
                for name in ['REASSEMBLY.json','reassemble.py']+parts:shutil.copyfile(ROOT/name,dst/name)
                m=json.loads(json.dumps(meta));sentinel=b'preserve this existing output';existing={};mtime=None
                if scenario in ('corrupt_first','corrupt_late'):
                    group=m['groups'][0 if scenario=='corrupt_first' else 1];name=group['parts'][0 if scenario=='corrupt_first' else -1]['file'];b=bytearray((dst/name).read_bytes());b[0]^=1;(dst/name).write_bytes(b)
                elif scenario in ('conflict_first','conflict_late','correct_first_conflicting_late'):
                    name='t7.bin' if scenario=='conflict_first' else 'blind8.bin';(dst/name).write_bytes(sentinel);existing[name]=sentinel
                    if scenario=='correct_first_conflicting_late':
                        (dst/'t7.bin').write_bytes(table);os.utime(dst/'t7.bin',ns=(1000000,1000000));mtime=(dst/'t7.bin').stat().st_mtime_ns;existing['t7.bin']=table
                elif scenario=='bad_late_hash':m['groups'][1]['outputs'][0]['sha256']='0'*64
                elif scenario=='nonlocal_late_output':m['groups'][1]['outputs'][0]['file']='../escape.bin'
                elif scenario=='duplicate_output':m['groups'][1]['outputs'][0]['file']='t7.bin'
                elif scenario=='nonlocal_late_chunk':m['groups'][1]['parts'][0]['file']='../escape.bin'
                (dst/'REASSEMBLY.json').write_text(json.dumps(m))
                command=[sys.executable]+(['-O'] if optimized else [])+[str(dst/'reassemble.py')];r=subprocess.run(command,capture_output=True,text=True)
                if scenario=='roundtrip':
                    need(r.returncode==0 and json.loads(r.stdout)['status']=='PASS','isolated complete reconstruction')
                    for name,raw in banks.items():need((dst/name).read_bytes()==raw,'complete isolated bytes')
                    again=subprocess.run(command,capture_output=True,text=True);need(again.returncode==0 and again.stdout==r.stdout,'idempotent paired outputs')
                else:
                    need(r.returncode!=0,'required refusal')
                    for name in banks:
                        if name in existing:need((dst/name).read_bytes()==existing[name],'existing output preserved')
                        else:need(not(dst/name).exists(),'no early output before late validation')
                    if mtime is not None:need((dst/'t7.bin').stat().st_mtime_ns==mtime,'valid first output was not rewritten before late refusal')
                controls.append(dict(optimized=optimized,scenario=scenario,status='PASS'))
    need(pins=={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},'packaging frozen during review')
    return dict(status='PASS',total_original_bytes_compared=sum(map(len,banks.values())),complete_seven_vectors_checked=len(table),complete_eight_words_checked=len(codes),physical_deletion_indices_checked=len(codes)*8,streams=streams,controls=controls,source_sha256=pins,scope='Lossless complete finite-certificate banks; mathematical completeness is separately audited, not inferred from packaging hashes.')
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
