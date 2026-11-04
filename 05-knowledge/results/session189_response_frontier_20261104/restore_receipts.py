#!/usr/bin/env python3
"""Verify and restore losslessly compressed Session189 receipts (stdlib only).

Run: python3 restore_receipts.py ARCHIVES_CP02.json [--verify-only]
Paths in each manifest are relative to that manifest's directory.
Existing files are accepted only when they already have the exact expected bytes.
"""
import argparse
import base64
import gzip
import hashlib
import json
import lzma
from pathlib import Path


def check(condition, message):
    if not condition:
        raise ValueError(message)


def safe_path(root, name):
    path = Path(name)
    check(not path.is_absolute() and '..' not in path.parts,
          'Unsafe relative path: ' + name)
    result = (root / path).resolve()
    check(result.is_relative_to(root.resolve()), 'Path escapes root: ' + name)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    check(manifest['format'] in ('deterministic-gzip-base64-parts-v1', 'lossless-compressed-base64-parts-v2'),
          'Unsupported receipt format')
    root = args.manifest.resolve().parent
    verified = []
    for record in manifest['archives']:
        text = ''.join(safe_path(root, part).read_text().strip()
                       for part in record['parts'])
        compressed = base64.b64decode(text, validate=True)
        check(len(compressed) == record.get('compressed_bytes', record.get('gzip_bytes')), 'Compressed byte mismatch')
        check(hashlib.sha256(compressed).hexdigest() == record.get('compressed_sha256', record.get('gzip_sha256')),
              'Compressed SHA-256 mismatch')
        codec = record.get('compression', 'gzip')
        check(codec in ('gzip', 'xz'), 'Unsupported compression')
        data = lzma.decompress(compressed) if codec == 'xz' else gzip.decompress(compressed)
        check(len(data) == record['bytes'], 'Original byte mismatch')
        check(hashlib.sha256(data).hexdigest() == record['sha256'],
              'Original SHA-256 mismatch')
        blob = b'blob ' + str(len(data)).encode() + b'\0' + data
        check(hashlib.sha1(blob).hexdigest() == record['git_blob_sha'],
              'Original Git blob hash mismatch')
        output = safe_path(root, record['path'])
        if output.exists():
            check(output.read_bytes() == data,
                  'Existing output differs: ' + str(output))
        elif not args.verify_only:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(data)
        verified.append({'path': record['path'], 'bytes': len(data),
                         'sha256': record['sha256']})
    print(json.dumps({'status': 'passed', 'verify_only': args.verify_only,
                      'verified_archives': verified}, indent=2))


if __name__ == '__main__':
    main()
