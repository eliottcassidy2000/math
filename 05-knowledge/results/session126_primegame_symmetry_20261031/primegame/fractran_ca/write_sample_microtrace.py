#!/usr/bin/env python3
"""Regenerate the compact exact first130-step CA sample, standard-library only."""
from pathlib import Path
import json
from primegame_ca import encode,ca_step,observe,ready,HEAD_STATES
ROOT=Path(__file__).resolve().parent
configuration=encode(2)
with (ROOT/'seed2_first130_microsteps.jsonl').open('w') as stream:
    for t in range(131):
        info=observe(configuration,False)
        row={'CA_step':t,'head_position':info['head_position'],
             'head_state':HEAD_STATES[info['head_state']-1],
             'READY':ready(configuration),'cells':sorted(configuration.items())}
        if row['READY']:row['observed_integer']=info['integer']
        stream.write(json.dumps(row,separators=(',',':'))+'\n')
        if t<130:configuration=ca_step(configuration)
if observe(configuration)['integer']!=4:
    raise RuntimeError('sample endpoint is not the first prime output4')
print(json.dumps({'status':'PASS','CA_updates':130,'final_integer':4}))
