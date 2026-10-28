"""Replay saved banks and exact finite quotients; --full rebuilds every C++ producer."""
from array import array
from collections import Counter
from pathlib import Path
import argparse, gzip, hashlib, importlib.util, json, subprocess, sys, tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from closed_flags import intern, canonical_probability

def need(b,message):
    if not b:
        raise RuntimeError(message)

def read(name):
    return json.loads((HERE/name).read_text(encoding='utf-8'))

def without_times(value):
    if isinstance(value,dict):
        return {k:without_times(v) for k,v in value.items() if k!='seconds'}
    if isinstance(value,list):
        return [without_times(v) for v in value]
    return value

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');args=ap.parse_args()
    pins=read('source_pins.json')
    for path,expected in pins['lf_normalized_sha256'].items():
        data=(ROOT/path).read_bytes().replace(b'\r\n',b'\n')
        need(hashlib.sha256(data).hexdigest()==expected,'source pin '+path)
    manifest=read('bank_manifest.json');raw={}
    for name,meta in manifest['artifacts'].items():
        data=(HERE/name).read_bytes();plain=gzip.decompress(data)
        need(len(data)==meta['bytes'] and len(plain)==meta['uncompressed_bytes'],'bank lengths')
        need(hashlib.sha256(data).hexdigest()==meta['sha256'],'compressed bank checksum')
        need(hashlib.sha256(plain).hexdigest()==meta['uncompressed_sha256'],'decoded bank checksum')
        raw[name[:-3]]=plain
    need(len(raw['parents9.u64'])==191536*8,'parent dictionary complete length')
    address=0;delta=0;shift=0;count=0
    for b in raw['flags10.varint']:
        delta |= (b&127)<<shift
        if b&128:
            shift+=7;need(shift<63,'bounded varint');continue
        need(count==0 or delta>0,'strict accepted address order')
        address+=delta;need(address>>9<191536,'decoded parent address')
        count+=1;shift=delta=0
    need(shift==0 and count==9733056,'complete flag stream')
    generation=read('generation.json');audit=read('verification10.json')
    for key in ('classes','automorphism_order_histogram','labelled_mass'):
        need(generation['results'][-1][key]==audit[key],'full census agreement '+key)
    need(audit['status']=='PASS' and audit['full_class_uniqueness'],'saved full audit')
    carrier_counts=[]
    for n in range(2,8):
        payload=json.loads(gzip.decompress((HERE/f'closed_flags{n}.json.gz').read_bytes()))
        q=payload['pointed_id_by_state'];meta=payload['pointed_types'];groups=[[] for _ in meta]
        for s,t in enumerate(q):groups[t].append(s)
        need(all(groups),'every pointed type occurs')
        colors,k=intern((meta[t][2],meta[t][0]) for t in q)
        total=1<<((n-1)*(n-2)//2)
        need(len(q)==n*total,'full rooted path universe')
        tau=[((s//n)^(total-1))*n+s%n for s in range(len(q))]
        while True:
            profiles,_=intern(canonical_probability(Counter(colors[s] for s in group)) for group in groups)
            new,nk=intern((colors[s],colors[tau[s]],profiles[q[s]]) for s in range(len(q)))
            if nk==k:
                need(new==colors,'stable exact partition');break
            colors,k=new,nk
        need(colors==payload['closed_color_by_state'],'complete minimum classifier replay')
        carrier_counts.append(k)
    expected=read('structure_audit.json')
    result=subprocess.run([sys.executable,'-B',str(HERE/'check_structure.py')],check=True,capture_output=True,text=True)
    need(json.loads(result.stdout)==expected,'literal structural audit replay')
    if args.full:
        with tempfile.TemporaryDirectory(prefix='consecutive-full-') as tmp:
            tmp=Path(tmp)
            for name,source in (('generate','generate.cpp'),('verify_bank','verify_bank.cpp'),
                                ('flag_states','flag_states.cpp'),('pair_bridge','pair_bridge.cpp')):
                subprocess.run(['g++','-O3','-std=c++17',str(HERE/source),'-o',str(tmp/(name+'.exe'))],check=True)
            subprocess.run([str(tmp/'generate.exe'),'10',str(tmp)],check=True)
            need(without_times(json.loads((tmp/'generation.json').read_text()))==without_times(generation),
                 'full generation replay')
            for name,data in raw.items():
                need((tmp/name).read_bytes()==data,'generated bank byte agreement '+name)
            subprocess.run([str(tmp/'verify_bank.exe'),'10',str(tmp),str(tmp/'verification10.json')],check=True)
            need(without_times(json.loads((tmp/'verification10.json').read_text()))==without_times(audit),
                 'full independent canonical audit replay')
            subprocess.run([str(tmp/'flag_states.exe'),str(tmp)],check=True)
            for n in range(2,8):
                a=array('I');a.frombytes((tmp/f'flag_states{n}.u32').read_bytes())
                if sys.byteorder!='little':a.byteswap()
                payload=json.loads(gzip.decompress((HERE/f'closed_flags{n}.json.gz').read_bytes()))
                need(list(a)==payload['pointed_id_by_state'],'full pointed classifier replay')
                need(json.loads((tmp/f'flag_meta{n}.json').read_text())==payload['pointed_types'],'pointed metadata replay')
            subprocess.run([str(tmp/'pair_bridge.exe'),str(tmp/'pair_bridge.json')],check=True)
            need(json.loads((tmp/'pair_bridge.json').read_text())==read('pair_bridge.json'),'full doubled correspondence replay')
    print(json.dumps({'status':'PASS','full_rebuild':args.full,'bank_classes':count,
                      'minimum_carrier_counts':carrier_counts,'source_pins':len(pins['lf_normalized_sha256'])}))

if __name__=='__main__':main()
