#!/usr/bin/env python3
"""Independent source jobs; process parallelism only, no agents or shared writes."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import scout

HERE=Path(__file__).resolve().parent


def one(job):
    i,row=job
    result=scout.search(tuple(row['target']),[row['clone']])
    (HERE/f'next_{i}.json').write_text(json.dumps(result,indent=2)+'\n')
    return i,result['statistics']


if __name__=='__main__':
    data=json.loads((HERE/'scout31.json').read_text())['programmes']
    # 7 is marked-isomorphic to0, and11 to10; those identifications are
    # logged separately and never treated as free physical operations.
    jobs=[(i,row) for i,row in enumerate(data) if i not in (0,1,7,11)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        tasks=[pool.submit(one,job) for job in jobs]
        for task in as_completed(tasks):
            print(json.dumps(task.result()),flush=True)
