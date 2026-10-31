#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from core import Core,DEFAULT_BASE
p=argparse.ArgumentParser(description='Construct and stream-check the certified cycle; O(N) memory/time.')
p.add_argument('N',type=int);p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);a=p.parse_args()
print(json.dumps(Core(a.base_dir).cycle(a.N),indent=2))
