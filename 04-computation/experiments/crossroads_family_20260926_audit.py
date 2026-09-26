"""Reproduce four family lanes. --record intentionally refreshes stdout files."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    for lane in ('flow', 'orbits', 'squares', 'fractal'):
        stem = 'crossroads_family_20260926_'+lane
        script = root/'04-computation'/'experiments'/(stem+'.py')
        target = root/'05-knowledge'/'results'/(stem+'.out')
        normal = subprocess.run([sys.executable, '-B', str(script)], cwd=root,
                                capture_output=True, check=True).stdout
        optimized = subprocess.run([sys.executable, '-B', '-O', str(script)], cwd=root,
                                   capture_output=True, check=True).stdout
        if normal != optimized:
            raise RuntimeError(lane+' normal/optimized mismatch')
        if args.record:
            target.write_bytes(normal)
        if normal.replace(b'\r\n', b'\n') != target.read_bytes().replace(b'\r\n', b'\n'):
            raise RuntimeError(lane+' retained-output mismatch')
        print(lane+': normal = optimized = retained stdout')
    print('PASS: four exact lane replays; claims remain scoped in their notes.')


if __name__ == '__main__':
    main()
