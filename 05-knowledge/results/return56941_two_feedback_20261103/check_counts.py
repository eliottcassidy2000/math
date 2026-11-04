"""Fast deterministic replay, distinct from workload-dependent benchmarks."""
from pathlib import Path
from hashlib import sha256
from importlib.util import spec_from_file_location,module_from_spec
import json,sys
from count_prefix import prefix,need
from bucket_prefix import prefix as bucket
HERE=Path(__file__).resolve().parent;sys.set_int_max_str_digits(0)
def main():
    spec=spec_from_file_location('session87',HERE.parent/'session87_four_vertices_enumeration_20261028/fast_counts.py');old=module_from_spec(spec);spec.loader.exec_module(old)
    for n in range(31):
        expected=old.degree_bucket_dp(n)[0]
        need(prefix(n,audit=True)[0]==prefix(n,'loop','python',True)[0]==expected,'common denominator controls')
        for tail in ('fold','convolution'):
            for backend in ('flint','python'):need(bucket(n,tail,backend)[0]==expected,'degree bucket controls')
    a=list(map(int,json.loads((HERE/'a000568_256.json').read_text())));b=list(map(int,json.loads((HERE/'a000568_200.json').read_text())))
    need(a[:201]==b and len(a)==257,'overlapping exact prefixes')
    # Independently published incoming receipt, commit83d4954e15.
    incoming=sha256(('\n'.join(hex(x) for x in a)+'\n').encode()).hexdigest()
    need(incoming=='f32a38d60933e6525ffa0c8c08eb8a8ca786a2b1f2a6bf85990ed05cf00f579f','incoming all-degree hex-line hash')
    n=a[-1];raw=n.to_bytes((n.bit_length()+7)//8,'big');h=sha256(raw).hexdigest()
    need(h=='2b415fefb721acef1bd9e740d4d53b57b4d6aa42db0921b33a0a43d0f99a8ef3','incoming full a256 value')
    need(n.bit_length()==30957 and (n&-n).bit_length()-1==131,'incoming exact binary length and valuation')
    report=dict(status='FINITE_EXACT_NORMALIZED_OUTPUT_CONTROLS',small_prefix_bound=30,small_variants=7,full_prefix_bound=256,
        incoming_commit='83d4954e15',incoming_source='orbit_spectrum_session_20261103/counting/full256.json',
        prefix_sha256_hex_lines=incoming,a256_sha256_big_endian=h,a256_bits=n.bit_length(),a256_v2=131,
        raw_decimal_file_hashes={f'a000568_{N}.json':sha256((HERE/f'a000568_{N}.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest() for N in (200,256)})
    (HERE/'count_controls.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS full incoming256 identity and seven small-prefix variants')
if __name__=='__main__':main()
