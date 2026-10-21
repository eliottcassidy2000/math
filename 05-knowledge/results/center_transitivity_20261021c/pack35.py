"""Losslessly package the separate discovery output; no exclusion logic."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import zlib

def main():
    p = argparse.ArgumentParser()
    p.add_argument('raw')
    args = p.parse_args()
    words = {}
    packed = bytearray()
    stream = hashlib.sha256()
    nodes = max_nodes = max_words = 0
    lines = Path(args.raw).read_text(encoding='ascii').splitlines()
    if len(lines) != 1 << 17:
        raise ValueError('complete frame count')
    for mask, line in enumerate(lines):
        fields = line.split()
        if len(fields) != 11 or fields[:2] != ['0', str(mask)]:
            raise ValueError('ordered literal frame schema')
        for word in fields[2:9]:
            if len(word) != 35 or set(word) != {chr(65+i) for i in range(35)}:
                raise ValueError('literal permutation')
            if word not in words:
                words[word] = len(words)
            packed.extend(words[word].to_bytes(2, 'little'))
            stream.update(word.encode('ascii'))
        a, b = map(int, fields[9:])
        nodes += a
        max_nodes = max(max_nodes, a)
        max_words = max(max_words, b)
    out = {'schema': 'C35-seven-black-v1', 'frames': len(lines),
           'words': list(words), 'index_bytes': len(packed),
           'index_sha256': hashlib.sha256(packed).hexdigest(),
           'expanded_word_sha256': stream.hexdigest(),
           'index_base85_zlib': base64.b85encode(zlib.compress(packed, 9)).decode('ascii'),
           'discovery_only': {'nodes': nodes, 'max_nodes': max_nodes, 'max_completed_words': max_words}}
    target = Path(__file__).resolve().parent / 'inputs35.json'
    target.write_text(json.dumps(out, sort_keys=True, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(len(words), len(packed), target.stat().st_size)

if __name__ == '__main__':
    main()
