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
    ('order12/check.py', 'order12/normal.json', 'order12/optimized.json'),
    ('order12/INDEPENDENT_REVIEW.md', 'order12/review-normal.json', 'order12/review-optimized.json'),
    ('order14/check.py', 'order14/normal.json', 'order14/optimized.json'),
    ('order14/INDEPENDENT_REVIEW.md', 'order14/review-normal.json', 'order14/review-optimized.json'),
    ('counting/check.py', 'counting/normal.json', 'counting/optimized.json'),
    ('counting/INDEPENDENT_REVIEW.md', 'counting/review-normal.json', 'counting/review-optimized.json'),
    ('arithmetic/check.py', 'arithmetic/normal.json', 'arithmetic/optimized.json'),
    ('arithmetic/INDEPENDENT_REVIEW.md', 'arithmetic/review-normal.json', 'arithmetic/review-optimized.json'),
    ('arithmetic/CONGRUENCE_REVIEW.md', 'arithmetic/congruence-normal.json', 'arithmetic/congruence-optimized.json'),
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
