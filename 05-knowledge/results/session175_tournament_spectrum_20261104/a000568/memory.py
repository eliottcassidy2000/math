#!/usr/bin/env python3
"""Traced Python allocation peaks; not RSS and not native allocator memory."""
import gc,hashlib,json,platform,sys,tracemalloc
from pathlib import Path
from projective import packed_count,need
here=Path(__file__).resolve().parent
rows=[]
for mode in [True,False]:
 gc.collect();tracemalloc.start()
 answer,stats=packed_count(150,mode)
 current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
 rows.append({'projective':mode,'order':150,'traced_peak_bytes':peak,'answer':str(answer)})
need(rows[0]['answer']==rows[1]['answer'],'memory comparison same answer')
r={'python':sys.version,'platform':platform.platform(),'source_sha256':hashlib.sha256((here/'projective.py').read_bytes()).hexdigest(),'rows':rows,'scope':'tracemalloc Python allocation peak, not RSS; collection before each run'}
(here/'memory-receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r))
