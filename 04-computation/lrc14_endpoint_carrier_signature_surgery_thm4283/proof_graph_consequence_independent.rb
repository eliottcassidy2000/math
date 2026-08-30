#!/usr/bin/env ruby
# Independent Ruby replay of the THM-4283 typed proof-graph union.

require 'digest'
require 'set'

FNV_OFFSET = 0xcbf29ce484222325
FNV_PRIME = 0x100000001b3
MASK64 = (1 << 64) - 1

def fail!(message)
  raise message
end

def read_pairs(path)
  rows = File.readlines(path, chomp: true).map.with_index do |line, index|
    fail!("blank row at #{path}:#{index + 1}") if line.empty?
    fields = line.split(',')
    fail!("malformed pair at #{path}:#{index + 1}") unless fields.length == 2
    pair = fields.map { |value| Integer(value, 10) }
    fail!("invalid pair at #{path}:#{index + 1}") unless pair[0] > 0 && pair[0] < pair[1]
    pair.freeze
  end
  fail!("unordered pair ledger #{path}") unless rows.each_cons(2).all? { |a, b| (a <=> b) == -1 }
  rows
end

def fnv(rows)
  state = FNV_OFFSET
  rows.flatten.each do |word|
    8.times do |byte|
      state ^= (word >> (8 * byte)) & 0xff
      state = (state * FNV_PRIME) & MASK64
    end
  end
  state
end

def sha(rows)
  Digest::SHA256.hexdigest(rows.map { |q, r| "#{q},#{r}\n" }.join)
end

fail!('usage: independent.rb RESIDUAL COMMON_UNION CARRIER_SCAN BOUNDARY_WITNESS') unless ARGV.length == 4
residual = read_pairs(ARGV[0])
common = read_pairs(ARGV[1])
scan_lines = File.readlines(ARGV[2], chomp: true)
boundary_lines = File.readlines(ARGV[3], chomp: true)

fail!('residual identity changed') unless residual.length == 23_373 &&
  fnv(residual) == 0xc6ab0ae49ee32273 &&
  sha(residual) == 'c3e5bf37887aa57af79cb166fce4a6e933e5daffc26dd8032fdfc52ce31240f3'
fail!('common identity changed') unless common.length == 640 &&
  fnv(common) == 0x45e9ecdf240e6417 &&
  sha(common) == '3246d76e82e9e19d07e5851810da3107ad8fe98a1dfbd087edb2b9c5d8b27fa1'

pair_rows = scan_lines.grep(/^PAIR /).map do |line|
  fields = line.split
  q, r = fields[1].split(',').map { |value| Integer(value, 10) }
  nested = Integer(fields[fields.index('NESTED_FAILURES') + 1], 10)
  [q, r, nested]
end
expected_scan = residual.select { |_, r| r >= 638 }.sort_by { |q, r| [-r, q] }
fail!('scan pair set changed') unless pair_rows.map { |q, r, _| [q, r] } == expected_scan
fail!('nested failure inside claimed band') unless pair_rows.all? { |_, r, n| r < 639 || n.zero? }
fail!('radius-638 hostile control changed') unless pair_rows.select { |_, r, _| r == 638 }
  .sum { |_, _, n| n } == 40
fail!('radius-638 witness replay changed') unless
  boundary_lines[0] == 'THM4283_ENDPOINT638_EXACT_RESPONSE_WITNESS_V1' &&
  boundary_lines[1].include?('FAILURES 40 FAILURE_FNV 917d107c4536efc9') &&
  boundary_lines[1].end_with?('EXACT_MINIMUM 9') &&
  boundary_lines[2] == 'WITNESSES 9 FNV 2b936529030e4bc MASKS 02203226 081e1084 08a89440 180a8281 18261042 18a0d040 1a82a200 202a9440 280a0a88' &&
  boundary_lines[3] == 'REPAIRED9006_FNV fdc1c57ae4dc1bb6 REPLAY_FAILURES 0' &&
  boundary_lines[-1] == 'VERDICT PASS EXACT_MINIMUM_NINE_AND_EXPLICIT_WITNESS'

carrier = residual.select { |_, r| r >= 638 }
residual_set = residual.to_set
common_set = common.to_set
carrier_set = carrier.to_set
fail!('typed common/carrier subset changed') unless common_set.subset?(residual_set) &&
  carrier_set.subset?(residual_set)
overlap = (common_set & carrier_set).to_a.sort
union = (common_set | carrier_set).to_a.sort
remaining = (residual_set - union.to_set).to_a.sort
max_r = remaining.map(&:last).max
top = remaining.select { |_, r| r == max_r }

expected = {
  'CARRIER638_644' => [carrier, 64, 0xc230e22462f7f3ab,
    '36baf6505a470f0bd63e306b5bafc895fb37841187d6a407cb67f7d5e5a2c2a3'],
  'OVERLAP' => [overlap, 13, 0xc9171a79d21e375d,
    '4660a98fcfc14a2f4319df70a3fec52f8848771312121b970bd017eceaa67f67'],
  'UNION' => [union, 691, 0x4b299b49d107a139,
    'c5646e81b3815bdef5168e36bcd76174065ed21339a5d8853d9efddc8fa3efae'],
  'FINAL_RESIDUAL' => [remaining, 22_682, 0xf7563445f15efebf,
    '7d0044bc4c32f08b9d09dca420171a05666d26e03f38fbc48a9baa1fcb027102']
}
expected.each do |label, (rows, count, expected_fnv, expected_sha)|
  fail!("#{label} identity changed") unless rows.length == count &&
    fnv(rows) == expected_fnv && sha(rows) == expected_sha
end
fail!('final boundary changed') unless max_r == 637 && top == [
  [100, 637], [294, 637], [520, 637]
]

puts 'THM4283_PROOF_GRAPH_INDEPENDENT_RUBY_V1'
expected.each do |label, (rows, _, _, _)|
  puts format('%s count=%d fnv=%016x sha256=%s', label, rows.length, fnv(rows), sha(rows))
end
puts "FINAL_BOUNDARY max_r=#{max_r} top=#{top.map { |q, r| "#{q},#{r}" }.join(';')}"
puts 'VERDICT PASS INDEPENDENT_RUBY_SET_REPLAY'
