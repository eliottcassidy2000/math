"""Package all VT21 literal words, with the full-target stream commitment."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import zlib

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('raw',type=Path)
    args=parser.parse_args()
    raw=args.raw.read_bytes()
    if len(raw)!=8+2048*7*48 or raw[:8]!=b'VT21KEY1':
        raise ValueError('complete VT21 producer stream')
    words=b''.join(raw[i:i+21] for i in range(8,len(raw),48))
    for i in range(0,len(words),21):
        if set(words[i:i+21])!=set(range(65,86)):
            raise ValueError('literal21 permutation')
    result={'schema':'VT21-two-groups-v1','frames':2048,'word_bytes':len(words),
            'word_sha256':hashlib.sha256(words).hexdigest(),
            'full_stream_sha256':hashlib.sha256(raw).hexdigest(),
            'word_base85_zlib':base64.b85encode(zlib.compress(words,9)).decode('ascii'),
            'discovery_only':{'total_nodes':169513,'maximum_nodes':6147,'maximum_completed_words':13,'node_cap_per_frame':10000000}}
    (HERE/'inputs21.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    print('Packed',len(words),'word bytes.')


if __name__=='__main__':
    main()
