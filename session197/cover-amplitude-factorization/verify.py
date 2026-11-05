#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,subprocess,tempfile
P=Path(__file__).resolve().parent
expected={'complete_unlabeled_through7':{'graphs':532,'roots':3613,'minimum_covers':6919,'rank_histogram':[209,1053,1873,478,0,0,0]},'complete_labelled_through5':{'graphs':1099,'roots':5405,'minimum_covers':8747,'rank_histogram':[1215,3350,840,0,0,0,0]},'seeded_larger_controls':{'graphs':72,'roots':756,'minimum_covers':1336,'rank_histogram':[5,31,115,211,210,151,33]},'status':'PASS'}
with tempfile.TemporaryDirectory(prefix='cover-factor-',dir=P)as td:
 exe=Path(td)/'check';subprocess.run(['g++','-O3','-std=c++17',str(P/'verify.cpp'),'-o',str(exe)],check=True)
 out=subprocess.check_output([str(exe)],text=True)
 if json.loads(out)!=expected:raise RuntimeError('factorization census mismatch')
 print(out,end='')
