#!/usr/bin/env python3
"""Complete finite eight-quotient certificate; standard library only."""
from pathlib import Path
from tempfile import TemporaryDirectory
import hashlib,json,struct,subprocess,sys
ROOT=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def main():
    sources=['PROOF.md','DATA.md','verify.py','census.cpp','check_endpoint.cpp','literal.cpp','orbits.py'];pins={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources}
    outputs={}
    with TemporaryDirectory() as td:
        for name in ['census','check_endpoint','literal']:
            binary=Path(td)/name
            subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/(name+'.cpp')),'-o',str(binary)],check=True)
            raw=subprocess.check_output([str(binary)],cwd=ROOT);(ROOT/(name+'.json')).write_bytes(raw);outputs[name]=json.loads(raw)
    table=(ROOT/'t7.bin').read_bytes();need(len(table)==1<<21 and all(x<128 for x in table),'complete seven-vector bank')
    need(table.count(0)==50688,'all seven zero signatures')
    raw=(ROOT/'blind8.bin').read_bytes();need(len(raw)==30720*4,'complete eight residual bytes');codes=list(struct.unpack('<30720I',raw));need(codes==sorted(set(codes)) and all(x<1<<28 for x in codes),'exact complete residual codes')
    need(all(table[c>>7]==0 for c in codes),'deletion-zero prerequisite')
    census=outputs['census'];need(census['seven_zero_vectors']==50688 and census['extension_candidates']==50688*128 and census['eight_shared_cubic_blind']==30720 and census['prime_blind']==26880,'complete census counts')
    check=outputs['check_endpoint'];need(check=={'verified_blind_quotients':30720,'endpoint_xor_histogram':[[255,30720]]},'independent endpoint residual')
    lit=outputs['literal'];need(lit['small_literal_queries']==570 and lit['source_mod8']==5 and lit['target_mod8']==7 and lit['quotient_start_parity']==155 and lit['quotient_end_parity']==100,'literal physical controls')
    orbit_raw=subprocess.check_output([sys.executable,str(ROOT/'orbits.py')],cwd=ROOT);(ROOT/'orbits.json').write_bytes(orbit_raw);orbits=json.loads(orbit_raw)
    need([x['code'] for x in orbits]==[600626,600653,602967,1472106],'complete orbit representatives')
    need([x['orbit_size'] for x in orbits]==[1920,13440,13440,1920] and [x['prime'] for x in orbits]==[False,True,True,False] and [x['converse_orbit'] for x in orbits]==[3,2,1,0],'exact orbit controls')
    need(pins=={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in sources},'source frozen')
    data=['t7.bin','blind8.bin','census.json','check_endpoint.json','literal.json','orbits.json']
    return dict(status='PASS',all_eight_quotients=1<<28,all_seven_quotients=1<<21,seven_zero_vectors=50688,candidate_extensions=6488064,complete_shared_cubic_blind=30720,prime_blind=26880,endpoint_balanced_blind=0,blind_relabelling_orbits=4,literal_physical_queries=572,source_sha256=pins,data_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in data},scope='Finite-exact eight-label TT2 quotient-converse separation, with arbitrary fixed ambient lift; not a theorem for all even orders or unrelated quotient reconstruction.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
