#!/usr/bin/env python3
"""Larger-order comparison with inherited compressed Fraction recurrence."""
import importlib.util,time,json,sys
from pathlib import Path
from hashlib import sha256
from engines import integer_response,require
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('inherited',HERE.parents[1]/'tournament_future_boundary_20261020/verify.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
rows=[]
for n in [120,160,200]:
 t=time.perf_counter();a,stats=integer_response(n);ti=time.perf_counter()-t
 t=time.perf_counter();b,_=old.profile_dp(n);tf=time.perf_counter()-t
 require(a==b,'large prefix exact equality')
 t=time.perf_counter();r,rs=integer_response(n,binary_precision=16);tr=time.perf_counter()-t
 require(r==[x%65536 for x in a],'large binary precision equality')
 row={'n':n,'integer_seconds':ti,'inherited_fraction_seconds':tf,'speedup':tf/ti,'binary16_seconds':tr,'states':stats,'binary16_states':rs,'prefix_sha256':sha256(','.join(map(str,a)).encode()).hexdigest()}
 rows.append(row);print(json.dumps(row),flush=True)
HERE.joinpath('scale.json').write_text(json.dumps(rows,indent=2)+'\n')
