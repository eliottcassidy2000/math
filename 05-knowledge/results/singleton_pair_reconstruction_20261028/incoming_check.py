"""Integrate concurrent low bits and independently replay the exported ten-bank."""
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import argparse,gzip,json,subprocess,sys,tempfile
import reconstruct as r

HERE=Path(__file__).resolve().parent
INCOMING=HERE.parent/'consecutive_flags_order10_20261028'

def low_bits(values,layer,color):
    counts=Counter()
    for n,value in enumerate(values):
        if n==0:r.need(value==1,'empty channel');counts['empty']+=1;continue
        if n==2:r.need(value==1<<(layer+2*color-1),'order-two channel');counts['two']+=1;continue
        base=layer*(n-1)//2+color if n%2 else layer*n//2+2*color
        parity=1 if n%2 else (sum(gcd(i,n)==1 for i in range(1,n))//2)%2
        r.need(value%(1<<(base+1))==parity<<base,'incoming independent low-bit theorem')
        counts['odd_exact' if n%2 else 'even_equal' if parity else 'even_strict']+=1
    return dict(counts)

def main(full_bank=False):
    sys.set_int_max_str_digits(0)
    ordinary=json.loads((HERE/'counts-bucket-350.json').read_text())['values']
    paired=json.loads((HERE/'pair-square-receipt.json').read_text())['paired_counts_through_100']
    result=dict(status='PASS',ordinary_low_bits=low_bits(ordinary,1,0),paired_low_bits=low_bits(paired,2,1))
    a=(6,4,8,3);b=r.core.complement(a)
    r.need(b==(10,12,9,0),'literal native partner')
    r.need(bool(r.oracle.isomorphisms(a,r.core.converse(a))) and not r.oracle.isomorphisms(b,r.core.converse(b)),
           'self-converse family not native-complement closed')
    result['native_self_converse_hostile']=[list(a),list(b)]
    manifest=json.loads((INCOMING/'bank_manifest.json').read_text());decoded={}
    for name,info in manifest['artifacts'].items():
        raw=(INCOMING/name).read_bytes();data=gzip.decompress(raw)
        r.need(sha256(raw).hexdigest()==info['sha256'] and sha256(data).hexdigest()==info['uncompressed_sha256'],'incoming bank bytes')
        decoded[name[:-3]]=data
    result['bank_compressed_bytes']=manifest['combined_compressed_bytes']
    if full_bank:
        # No incoming artifact is rewritten; all executable/decoded files stay
        # in a new bounded temporary directory. No cleanup/delete is required.
        work=Path(tempfile.mkdtemp(prefix='singleton-ten-audit-'))
        for name,data in decoded.items():(work/name).write_bytes(data)
        exe=work/'verify.exe';receipt=work/'audit.json'
        subprocess.run(['g++','-O3','-std=c++17',str(INCOMING/'verify_bank.cpp'),'-o',str(exe)],check=True)
        subprocess.run([str(exe),'10',str(work),str(receipt)],check=True)
        bank=json.loads(receipt.read_text());expected=json.loads((INCOMING/'verification10.json').read_text())
        for key in ('classes','labelled_mass','full_class_uniqueness','automorphism_order_histogram'):
            r.need(bank[key]==expected[key],'independent replay of incoming full ten-bank')
        r.need(bank['classes']==ordinary[10]==9733056 and bank['labelled_mass']==1<<45,'incoming scalar/constructive agreement')
        result['full_bank_audit']=bank
        (HERE/'incoming-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--full-bank',action='store_true');args=ap.parse_args()
    main(args.full_bank)
