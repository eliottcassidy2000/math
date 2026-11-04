"""Independent finite controls for the scoped transitive-chart theorem.

All guards are explicit and remain active under python -O. No old receipts
are rewritten. The only output file is transitive_chart.json beside this file.
"""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def joint_pairs(f):
    ports = [i for i in range(f.bit_length()) if f >> i & 1]
    return sum((1 << (k-i-1)) + (k == i+1)
               for index, i in enumerate(ports) for k in ports[index+1:])


def joint_bits(f):
    g = j = e = 0
    for i in range(f.bit_length()):
        b = (f >> i) & 1
        j, g, e = j+b*(g+e), 2*g+b, b
    return j


def ordered_subset_signature(f):
    """Literal ordered two-slot covers: each induced TT has one path."""
    n = f.bit_length()
    full = (1 << n)-1
    fc = gc = jc = 0
    for first in range(full+1):
        second = full ^ first
        end_port = bool(first and f >> (first.bit_length()-1) & 1)
        start_port = bool(second and f >> ((second & -second).bit_length()-1) & 1)
        fc += end_port
        gc += start_port
        jc += end_port and start_port
    return (1, 1 << n, fc, gc, jc, 1, f & 1)


def attachment_rows(f, clock):
    """Independent literal adjacency for prefix, suffix, and root."""
    m = f.bit_length()
    n = clock+m+1
    root = n-1
    rows = [0]*n
    for u in range(n):
        for v in range(u+1, n):
            if v != root:
                # All prefix and suffix vertices lie in one ordered TT.
                rows[u] |= 1 << v
            elif u < clock or not (f >> (u-clock) & 1):
                rows[root] |= 1 << u
            else:
                rows[u] |= 1 << root
    return rows


def strong(rows):
    full = (1 << len(rows))-1
    for source in range(len(rows)):
        reached = 1 << source
        frontier = reached
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length()-1
            fresh = rows[vertex] & ~reached
            reached |= fresh
            frontier |= fresh
        if reached != full:
            return False
    return True


def hamilton_count(rows):
    """Subset DP on literal directed adjacency; no port-count identities."""
    n = len(rows)
    full = (1 << n)-1
    dp = [[0]*n for _ in range(full+1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full+1):
        for end, count in enumerate(dp[mask]):
            if not count:
                continue
            extensions = rows[end] & ~mask
            while extensions:
                bit = extensions & -extensions
                extensions ^= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[full])


def value(f, t):
    need(t >= 1, 'illegal prefix clock')
    return 1-joint_pairs(f)+(1 << t)*f


def shifted_mask(f, k):
    if k <= 0:
        return (3*f) << (-k)
    divisor = 1 << k
    return 3*f//divisor if 3*f % divisor == 0 else None


def main():
    subset_checks = shift_checks = attachment_checks = 0
    for f in range(1, 256):
        sig = ordered_subset_signature(f)
        need(sig[2] == f, ('literal marginal', f, sig))
        need(sig[3] == int(bin(f)[2:][::-1], 2), ('literal reverse marginal', f, sig))
        need(sig[4] == joint_pairs(f) == joint_bits(f), ('literal joint', f, sig))
        need(sig[4] % 2 == 0, ('joint parity', f))
        subset_checks += 1
        d = 3*joint_pairs(f)-joint_pairs(3*f)-2
        for a in range(5):
            need(joint_pairs(f << a) == sig[4], ('shift invariance', f, a))
        for k in range(-4, 5):
            target = shifted_mask(f, k)
            if target is None:
                continue
            need(joint_pairs(target) == joint_pairs(3*f), ('target correction', f, k))
            clock = max(1, 1-k)
            for t in (clock, clock+1):
                need(value(target, t+k) == 3*value(f, t)+d,
                     ('fixed clock shift', f, k, t))
            shift_checks += 1

    for f in range(1, 32):
        for t in (1, 2, 3):
            rows = attachment_rows(f, t)
            need(strong(rows), ('literal attachment strongness', f, t))
            need(hamilton_count(rows) == value(f, t), ('literal full count', f, t))
            attachment_checks += 1

    controls = []
    for f, expected in ((1, -4), (3, 0), (5, -10), (7, 4), (11, 6), (21, -28)):
        j, target_j = joint_pairs(f), joint_pairs(3*f)
        d = 3*j-target_j-2
        need(d == expected, ('named control', f, d, expected))
        controls.append(dict(f=f, j=j, target_f=3*f, target_j=target_j, digit=d))

    digit_two = []
    digit_four = []
    odd_scan = 0
    for f in range(1, 1 << 16, 2):
        j, target_j = joint_pairs(f), joint_pairs(3*f)
        need(j == joint_bits(f), ('scan source evaluators', f))
        need(target_j == joint_bits(3*f), ('scan target evaluators', f))
        d = 3*j-target_j-2
        need(d % 2 == 0, ('digit parity', f))
        if d == 2:
            digit_two.append(f)
        if d == 4:
            digit_four.append(f)
        odd_scan += 1
    need(digit_two == [], ('bounded digit two classification', digit_two))
    need(digit_four == [7], ('bounded digit four classification', digit_four))
    need((subset_checks, shift_checks, attachment_checks, odd_scan) == (255, 1511, 93, 32768),
         'declared finite universe changed')

    source_path = HERE.parent/'session182_physical_returns_20261104'/'induction'/'PROOF.md'
    source_pin = hashlib.sha256(source_path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    report = dict(
        status='FINITE_EXACT_INDEPENDENT_TRANSITIVE_CHART_CONTROLS',
        scope='All-order uniqueness is proved in TRANSITIVE_CHART.md; no global digit classification, two-digit closure, or native legality claimed.',
        source=dict(path='../session182_physical_returns_20261104/induction/PROOF.md',
                    read_commit='b486271d24', normalized_sha256=source_pin),
        ordered_subset_masks=dict(first=1, last=255, count=subset_checks),
        fixed_clock_shifts=dict(first=-4, last=4, integral_cases=shift_checks,
                                clocks_per_case=2),
        literal_attachments=dict(first_mask=1, last_mask=31, prefix_lengths=[1, 2, 3],
                                 count=attachment_checks, independent_strongness=True),
        controls=controls,
        finite_odd_scan=dict(lower_inclusive=1, upper_exclusive=1 << 16,
                             count=odd_scan, digit_two=digit_two, digit_four=digit_four,
                             all_order_classification_claimed=False),
        reproduction=[
            'python 05-knowledge/results/multiport_defect_grammar_20261103/check_transitive_chart.py',
            'python -O 05-knowledge/results/multiport_defect_grammar_20261103/check_transitive_chart.py'],
        guards='Explicit RuntimeError checks; no Python assertions.')
    (HERE/'transitive_chart.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('PASS: 255 literal TT subset signatures; 1511 shifted-clock cases; 93 literal attachments; 32768 odd masks')
    print('FINITE ONLY: no digit-2 odd mask below 65536; unique digit-4 odd mask below 65536 is 7')


if __name__ == '__main__':
    main()
