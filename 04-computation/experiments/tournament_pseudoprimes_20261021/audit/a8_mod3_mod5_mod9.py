# AUDIT a8 (independent auditor, 2026-10-21): THM-4746 statement 3 (mod 3, mod 5) for m <= 200 and the FINITE-EXACT
# mod-9 pattern of sec. 5 / HYP-9358(d) for odd m <= 1401, from a2_seq_mod45.txt (s(m) mod 45, m = 1..1401).
s = {}
for line in open("a2_seq_mod45.txt"):
    m, v = map(int, line.split()); s[m] = v
c = lambda m, M: (1 % M if m == 1 else (pow(2, (m - 1) // 2, M) if m % 2 else 0))
bad3 = [m for m in range(1, 201) if (s[m] - ((-1) ** ((m - 1) // 2) if m % 2 else 0)) % 3]
bad3b = [m for m in range(1, 201) if m % 2 and (s[m] - pow(2, (m - 1) // 2, 3)) % 3]
print("mod 3: m <= 200 violating s(m) = sin-coefficient (odd: (-1)^((m-1)/2), even: 0):", bad3, "; odd m with s != 2^((m-1)/2):", bad3b)
bad5 = [m for m in range(1, 201, 2) if (s[m] - pow(2, (m - 1) // 2, 5)) % 5]
print("mod 5: odd m <= 200 with s(m) != 2^((m-1)/2):", bad5)
ev2 = [m for m in range(2, 201, 4) if s[m] % 5]
ev0 = sorted(set((s[m] % 5) for m in range(4, 201, 4)))
print("mod 5: m = 2 mod 4 with s(m) != 0:", ev2, "; residues of s(m) mod 5 for 4|m:", ev0, "(+-24 = +-4 mod 5)")
bad9 = [m for m in range(1, 1402, 2) if (s[m] - c(m, 9)) % 9]
print("mod 9: odd m <= 1401 with D(m) != 0 mod 9: count", len(bad9))
print("   all = 1 mod 6 and >= 7:", all(m % 6 == 1 and m >= 7 for m in bad9),
      "; every m = 1 mod 6, 7 <= m <= 1401 is bad:", set(range(7, 1402, 6)) == set(bad9))
print("   residue of D(m) mod 9 on the bad class, by m mod 18 / mod 54:",
      sorted(set(((m % 54), (s[m] - c(m, 9)) % 9) for m in bad9))[:40])
# even m mod 9 (not claimed; for the record)
print("mod 9: even m <= 60, s(m) mod 9:", [s[m] % 9 for m in range(2, 61, 2)])
