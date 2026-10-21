# AUDIT a7 (independent auditor, 2026-10-21): the auditor's lemma  s(p^2) - c(p^2) = eps * D(p)  (mod p^2),
# eps = (2/p), which makes "p^2 is a Camion pseudoprime" EQUIVALENT to "D(p) = 0 mod p^2" (THM-4746 sec. 5 asserts the
# equivalence parenthetically; HYP-9358(a) treats it as finite evidence).  Checked by DIRECT computation of s(p^2) mod p^2
# (a2_smod single p^2) for every odd prime p <= 97, and D(p) mod p^2 (a2_smod single p p^2).
import subprocess
def run(*a):
    out = subprocess.run(["./a2_smod", "single", *map(str, a)], capture_output=True, text=True).stdout
    return int(out.split("s(n) mod M = ")[1].split()[0]), int(out.split("c(n) mod M = ")[1].split()[0])
ok = True
for p in [q for q in range(3, 98) if all(q % d for d in range(2, int(q**.5) + 1))]:
    M = p * p
    s2, c2 = run(M)            # s(p^2), c(p^2) mod p^2
    s1, _ = run(p, M)          # s(p) mod p^2
    eps = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    Dp = (s1 - pow(2, (p - 1) // 2, M)) % M
    lhs = (s2 - c2) % M; rhs = (eps * Dp) % M
    ok &= (lhs == rhs)
    print(f"p={p:3d} eps={eps:+d}  s(p^2)-c(p^2) mod p^2 = {lhs:6d}   eps*D(p) mod p^2 = {rhs:6d}   equal={lhs == rhs}   p^2 passes={lhs == 0}")
print("IDENTITY HOLDS FOR ALL p <= 97" if ok else "IDENTITY FAILS")
