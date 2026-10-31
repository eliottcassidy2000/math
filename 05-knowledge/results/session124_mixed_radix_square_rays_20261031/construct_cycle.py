#!/usr/bin/env python3
import argparse,json
from pathlib import Path
from core import Constructor,DEFAULT_BASE,DEFAULT_SMALL,DEFAULT_Q49,threshold
p=argparse.ArgumentParser();p.add_argument('word');p.add_argument('N',type=int,nargs='?');p.add_argument('--base-dir',type=Path,default=DEFAULT_BASE);p.add_argument('--small-base-file',type=Path,default=DEFAULT_SMALL);p.add_argument('--q49-dir',type=Path,default=DEFAULT_Q49);a=p.parse_args()
print(json.dumps(Constructor(a.base_dir,a.small_base_file,a.q49_dir).cycle(a.N if a.N is not None else threshold(a.word),a.word),indent=2))
