#!/usr/bin/env python3
"""Independent lossless package roundtrip and refusal controls."""
import base64,gzip,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent/'covered_laminar'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    meta=json.loads((ROOT/'SIX_REASSEMBLY.json').read_text());need(len(meta['groups'])==1,'one bank group');group=meta['groups'][0];need(len(group['outputs'])==1,'one output')
    chunks=[];parts=[]
    for item in group['parts']:
        raw=(ROOT/item['file']).read_bytes();need(len(raw)==item['size'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'part hash and length');parts.append(item['file']);chunks.append(raw.strip())
    compressed=base64.b64decode(b''.join(chunks),validate=True);raw=gzip.decompress(compressed);original=(ROOT/'mod4_six_vectors.bin').read_bytes();item=group['outputs'][0]
    need(raw==original,'complete original byte equality');need(len(raw)==item['size'] and hashlib.sha256(raw).hexdigest()==item['sha256'],'output manifest integrity')
    controls=[]
    for optimized in [False,True]:
        for scenario in ['roundtrip','corrupt_first_chunk','corrupt_last_chunk','conflicting_output','late_conflict','bad_output_hash','nonlocal_output']:
            with tempfile.TemporaryDirectory(prefix='six_bank_packaging_') as td:
                dst=Path(td)
                for name in parts+['SIX_REASSEMBLY.json','reassemble_six.py']:shutil.copyfile(ROOT/name,dst/name)
                outfile=dst/'mod4_six_vectors.bin';sentinel=b'keep this existing output'
                if scenario.startswith('corrupt_'):
                    target=dst/(parts[0] if 'first' in scenario else parts[-1]);data=bytearray(target.read_bytes());data[0]^=1;target.write_bytes(data)
                elif scenario=='conflicting_output':outfile.write_bytes(sentinel)
                elif scenario=='late_conflict':
                    m=json.loads(json.dumps(meta));second=json.loads(json.dumps(m['groups'][0]));second['outputs'][0]['file']='late.bin';m['groups'].append(second);(dst/'SIX_REASSEMBLY.json').write_text(json.dumps(m));(dst/'late.bin').write_bytes(sentinel)
                elif scenario=='bad_output_hash':
                    m=json.loads(json.dumps(meta));m['groups'][0]['outputs'][0]['sha256']='0'*64;(dst/'SIX_REASSEMBLY.json').write_text(json.dumps(m))
                elif scenario=='nonlocal_output':
                    m=json.loads(json.dumps(meta));m['groups'][0]['outputs'][0]['file']='../escape.bin';(dst/'SIX_REASSEMBLY.json').write_text(json.dumps(m))
                cmd=[sys.executable]+(['-O'] if optimized else [])+[str(dst/'reassemble_six.py')];r=subprocess.run(cmd,capture_output=True,text=True)
                if scenario=='roundtrip':
                    need(r.returncode==0 and outfile.read_bytes()==raw,'isolated reconstruction');need(json.loads(r.stdout)['status']=='PASS','roundtrip receipt')
                    # Repeating with the correct existing output must remain byte-identical.
                    r2=subprocess.run(cmd,capture_output=True,text=True);need(r2.returncode==0 and outfile.read_bytes()==raw and r2.stdout==r.stdout,'idempotent exact output')
                else:
                    need(r.returncode!=0,'refusal exit code')
                    if scenario=='conflicting_output':need(outfile.read_bytes()==sentinel,'conflicting existing output unchanged')
                    else:need(not outfile.exists(),'no partial write before all validation')
                    if scenario=='late_conflict':need((dst/'late.bin').read_bytes()==sentinel,'late conflict unchanged')
                controls.append(dict(optimized=optimized,scenario=scenario,status='PASS'))
    files=['SIX_REASSEMBLY.md','SIX_REASSEMBLY.json','reassemble_six.py']+parts
    return dict(status='PASS',raw_bytes=len(raw),raw_sha256=hashlib.sha256(raw).hexdigest(),gzip_bytes=len(compressed),gzip_sha256=hashlib.sha256(compressed).hexdigest(),parts=len(parts),all_original_bytes_compared=len(raw),controls=controls,source_pins=[dict(file=name,size=(ROOT/name).stat().st_size,sha256=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()) for name in files])
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
