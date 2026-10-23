"""Losslessly package discovered words; commit to every literal full target."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import struct
import zlib

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('raw', type=Path)
    args = parser.parse_args()
    raw = args.raw.read_bytes()
    if raw[:8] != b'D15Q302!':
        raise ValueError('producer version')
    at = 8
    index = bytearray(b'D15IDX01')
    words, ids = [], {}
    counts = [0, 0]
    for mask in range(1 << 18):
        tag = raw[at]
        at += 1
        if tag not in (0, 1):
            raise ValueError('literal tag')
        index.append(tag)
        counts[tag] += 1
        if tag:
            for j in range(7):
                value = int.from_bytes(raw[at:at+8], 'little')
                text = ''.join(format(value >> (4*i) & 15, 'x') for i in range(15))
                if len(text) != 15 or set(text) != set('0123456789abcde'):
                    raise ValueError('literal Hamilton permutation')
                at += 22  # 8 word bytes and 14 full-target bytes.
                if text not in ids:
                    ids[text] = len(words)
                    words.append(text)
                index.extend(struct.pack('<H', ids[text]))
    if at != len(raw) or counts != [91392, 170752]:
        raise ValueError('complete declared universe')
    result = {
        'schema': 'D15-q3-full-key-v1', 'frames': 1 << 18,
        'score_rejected': counts[0], 'seven_black': counts[1],
        'words': words, 'index_bytes': len(index),
        'index_sha256': hashlib.sha256(index).hexdigest(),
        'expanded_bytes': len(raw),
        'expanded_sha256': hashlib.sha256(raw).hexdigest(),
        'index_base85_zlib': base64.b85encode(zlib.compress(index, 9)).decode('ascii'),
        'discovery_only': {'total_nodes': 9438508, 'maximum_nodes': 1085,
                           'maximum_completed_words': 42, 'node_cap_per_frame': 10000000}}
    (HERE/'inputs_q3.json').write_text(json.dumps(result, sort_keys=True, indent=2)+'\n', encoding='utf-8')
    print('Packed', len(raw), 'bytes;', len(words), 'distinct words;', len(index), 'index bytes.')


if __name__ == '__main__':
    main()
