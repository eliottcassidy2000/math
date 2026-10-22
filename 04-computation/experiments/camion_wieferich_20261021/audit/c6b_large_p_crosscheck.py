# c6b_large_p_crosscheck.py -- q(p) = D(p)/p mod p at large p by the audit's FLINT mod-p^2 formula s(p) = (p-1)! [x^(p-1)] T'/T^2
# (c5_cw_lists_flint.py), against the values printed by the session-code replay in c6_validate_cw_flint_np.out.
import flint, re, os
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'c5_cw_lists_flint.py')).read()
ns = {}; exec('import flint\n' + 'def primes' + src.split('# --- self-test')[0].split('def primes', 1)[1], ns)
sess = {int(a): int(b) for a, b in re.findall(r'p=(\d+): .*q\(p\) = (\d+)', open(os.path.join(here, 'c6_validate_cw_flint_np.out')).read())}
for p, qs in sorted(sess.items()):
    M = p * p; s = ns['s_mod_p2'](2, p); D = (s - pow(2, (p - 1) // 2, M)) % M
    print(f"p={p}: D(p) mod p = {D % p} (Euler: 0); q(p) by T'/T^2 mod p^2 = {D // p}; session code path {qs}: {'MATCH' if D // p == qs else 'MISMATCH'}")
