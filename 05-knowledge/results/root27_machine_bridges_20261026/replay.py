"""Replay all exact controls; only --refresh writes canonical receipts."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOBS = [
    ('check_incoming_s65.py', 'incoming-normal.json', 'incoming-optimized.json'),
    ('anchor/check.py', 'anchor/normal.json', 'anchor/optimized.json'),
    ('anchor/review_check.py', 'anchor/review-normal.json', 'anchor/review-optimized.json'),
    ('anchor/LEGENDRE_SIDECAR.md', 'anchor/legendre-normal.json', 'anchor/legendre-optimized.json'),
    ('atlas/check.py', 'atlas/normal.json', 'atlas/optimized.json'),
    ('atlas/INDEPENDENT_REVIEW.md', 'atlas/review-normal.json', 'atlas/review-optimized.json'),
    ('machines/check.py', 'machines/normal.json', 'machines/optimized.json'),
    ('machines/INDEPENDENT_REVIEW.md', 'machines/review-normal.json', 'machines/review-optimized.json'),
    ('compiler/check.py', 'compiler/normal.json', 'compiler/optimized.json'),
    ('compiler/INDEPENDENT_REVIEW.md', 'compiler/review-normal.json', 'compiler/review-optimized.json'),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run_job(job, refresh):
    source, normal, optimized = job
    path = HERE / source
    if path.suffix == '.md':
        blocks = re.findall(r'^```python\n(.*?)^```', path.read_text(encoding='utf-8'), re.M | re.S)
        need(len(blocks) == 1, (source, 'expected one Python control'))
        args = ['-c', blocks[0]]
    else:
        args = [str(path)]
    outputs = []
    for flags, name in (([], normal), (['-O'], optimized)):
        proc = subprocess.run([sys.executable, '-B', *flags, *args], cwd=ROOT,
                              capture_output=True, text=True, timeout=300)
        need(proc.returncode == 0, (source, proc.returncode, proc.stderr[-4000:]))
        data = json.loads(proc.stdout)
        need(data.get('status') == 'PASS', (source, 'non-PASS receipt'))
        output = json.dumps(data, indent=2, sort_keys=True) + '\n'
        outputs.append(output)
        target = HERE / name
        if refresh:
            target.write_bytes(output.encode('utf-8'))
        else:
            need(target.read_text(encoding='utf-8') == output, (source, 'receipt changed'))
    need(outputs[0] == outputs[1], (source, 'normal/optimized mismatch'))
    return source


def main():
    refresh = '--refresh' in sys.argv
    with ThreadPoolExecutor(max_workers=3) as pool:
        for source in pool.map(lambda job: run_job(job, refresh), JOBS):
            print('PASS', source, flush=True)
    print(json.dumps({'status':'PASS', 'paired_replays':len(JOBS), 'refresh':refresh}))


if __name__ == '__main__':
    main()
