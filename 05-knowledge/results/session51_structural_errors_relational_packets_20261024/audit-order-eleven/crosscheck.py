#!/usr/bin/env python3
"""Check frozen producer certificates using the independent audit algorithms."""
import collections
import hashlib
import json
from pathlib import Path
import audit
HERE=Path(__file__).resolve().parent

def main():
 pin=json.loads((HERE/'PRODUCER_PIN.json').read_text())
 audit.require(hashlib.sha256((HERE/'producer/MANIFEST.json').read_bytes()).hexdigest()==pin['manifest_sha256'],'frozen producer manifest pin')
 produced=json.loads((HERE/'producer/normal.json').read_text());independent=json.loads((HERE/'normal.json').read_text())
 records=[]
 for sector in produced['uniform_substitutions'].values():records.extend(sector['classes'])
 for sector in produced['five_cells']['nonuniform_sectors']:records.extend(sector['classes'])
 uniform=dict(produced['five_cells']['uniform_zero_control']);uniform['error']=0;records.append(uniform)
 records.extend(produced['triangle_cells']['classes'])
 audit.require(len(records)==1232,'producer certificate has1232 actual representatives')
 histogram=collections.Counter();pattern_histogram=collections.Counter()
 for record in records:
  t=tuple(record['rows']);audit.valid(t);audit.require(audit.strong(t),'producer representative independently strong')
  group=audit.automorphisms(t);e,sizes=audit.defect(t,group)
  audit.require(len(group)==record['automorphisms'] and e==record['error'],'producer full Aut and error independently recomputed')
  pattern=tuple(sorted(map(len,audit.refine(t))))
  audit.require(pattern!=(11,),'no regular eleven entered frozen producer residual')
  histogram[e]+=1;pattern_histogram[(pattern,e)]+=1
 inherited_hist=collections.Counter((tuple(sorted(map(len,audit.refine(tuple(r['rows']))))),r['error']) for r in independent['representatives'])
 audit.require(pattern_histogram==inherited_hist,'producer and independent sector distributions match')
 histogram.pop(0,None)
 audit.require(sorted(histogram.items())==[tuple(x) for x in independent['positive_nonregular_histogram']],'producer independent positive distribution exact match')
 audit.require(sum(e*n for e,n in histogram.items())==9166,'producer independent residual9166')
 # Every orbit size of an odd-order group is odd. Enumerate only integer
 # partitions to verify the permitted error support, not regular tournaments.
 support=sorted({11*(len(p)-1) for p in audit.oddparts(11) if len(p)>1})
 audit.require(support==[22,44,66,88,110],'regular eleven odd-orbit support')
 audit.require(produced['unknown']['positive_regular_eleven_error_support']==support,'producer unknown support agrees with independent parity deduction')
 audit.require(19480%22==10,'known total error residue modulo twenty two')
 return {'status':'PASS','producer_manifest_sha256':pin['manifest_sha256'],'independently_checked_producer_representatives':len(records),'nonregular_positive_distribution':sorted(histogram.items()),'regular_eleven_allowed_positive_support':support,'total_error_modulo_twenty_two':10,'total_eleven_error':'UNKNOWN regular-eleven contribution +19480','checks':dict(sorted(audit.CHECKS.items()))}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
