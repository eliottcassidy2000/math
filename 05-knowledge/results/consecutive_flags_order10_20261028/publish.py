"""Build/replay the complete census, or pack an already completed run."""
import argparse,gzip,hashlib,json,shutil,subprocess,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(b,message):
    if not b: raise RuntimeError(message)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pack-from',type=Path);args=ap.parse_args()
    workspace=args.pack_from or Path(tempfile.mkdtemp(prefix='consecutive-order10-'))
    workspace.mkdir(parents=True,exist_ok=True)
    if not args.pack_from:
        for name,source in (('generate','generate.cpp'),('verify','verify_bank.cpp')):
            subprocess.run(['g++','-O3','-std=c++17',str(HERE/source),'-o',str(workspace/(name+'.exe'))],check=True)
        subprocess.run([str(workspace/'generate.exe'),'10',str(workspace)],check=True)
        subprocess.run([str(workspace/'verify.exe'),'10',str(workspace),str(workspace/'verification10.json')],check=True)
    census=json.loads((workspace/'generation.json').read_text(encoding='utf-8'))
    audit=json.loads((workspace/'verification10.json').read_text(encoding='utf-8'))
    last=census['results'][-1]
    for key in ('classes','labelled_mass','automorphism_order_histogram'):
        need(last[key]==audit[key],'independent full-bank '+key)
    need(audit['full_class_uniqueness'] and audit['status']=='PASS','full-bank audit')
    need([r['classes'] for r in census['results']]==[1,1,2,4,12,56,456,6880,191536,9733056],
         'classical count reference is audit-only')
    manifest={'status':'PASS','artifacts':{},'decoded_tournament_order':10,
              'flag_parent_order':9,'address_bits':27,
              'parent_index_bits':18,'attachment_mask_bits':9,
              'decoded_dense_bank_bytes':8*last['classes']}
    for name in ('parents9.u64','flags10.varint'):
        raw=(workspace/name).read_bytes();zipped=gzip.compress(raw,mtime=0)
        (HERE/(name+'.gz')).write_bytes(zipped)
        manifest['artifacts'][name+'.gz']={'bytes':len(zipped),'uncompressed_bytes':len(raw),
             'sha256':hashlib.sha256(zipped).hexdigest(),'uncompressed_sha256':hashlib.sha256(raw).hexdigest()}
    for name in ('generation.json','verification10.json'):
        shutil.copyfile(workspace/name,HERE/name)
    manifest['combined_compressed_bytes']=sum(v['bytes'] for v in manifest['artifacts'].values())
    (HERE/'bank_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(manifest))

if __name__=='__main__':main()
