#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from core42 import create_core,DEFAULT_CORE,DEFAULT_BASE
p=argparse.ArgumentParser();p.add_argument('N',type=int);p.add_argument('--core-dir',type=Path,default=DEFAULT_CORE);p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);a=p.parse_args()
print(json.dumps(create_core(a.core_dir,a.base_dir).cycle(a.N),indent=2))
