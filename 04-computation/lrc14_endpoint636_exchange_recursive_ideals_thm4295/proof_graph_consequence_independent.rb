#!/usr/bin/env ruby
# Independent exact-set replay for the carrier10 + three signature ideals.

require 'digest'
require 'fileutils'

EMPTY = 0xcbf29ce484222325
PRIME = 0x100000001b3
MASK64 = (1 << 64) - 1

def read_pairs(path)
  rows = File.readlines(path, chomp: true).reject(&:empty?).map do |line|
    fields = line.split(',')
    raise "malformed row: #{line}" unless fields.length == 2
    fields.map { |field| Integer(field, 10) }.freeze
  end
  raise "unordered or duplicate rows: #{path}" unless rows == rows.uniq.sort
  rows
end

def bytes(rows)
  rows.map { |q, r| "#{q},#{r}\n" }.join
end

def fnv(rows)
  state = EMPTY
  rows.each do |pair|
    pair.each do |word|
      8.times do |byte|
        state ^= (word >> (8 * byte)) & 0xff
        state = (state * PRIME) & MASK64
      end
    end
  end
  state
end

def identity(label, rows)
  puts format('%s COUNT %d FNV %016x SHA256 %s', label, rows.length,
              fnv(rows), Digest::SHA256.hexdigest(bytes(rows)))
end

raise 'usage: independent LIVE CARRIER H19 H294 H372 OUT' unless ARGV.length == 6
live_path, carrier_path, h19_path, h294_path, h372_path, out_dir = ARGV
ledgers = {
  'LIVE' => read_pairs(live_path),
  'CARRIER10' => read_pairs(carrier_path),
  'H19' => read_pairs(h19_path),
  'H294' => read_pairs(h294_path),
  'H372' => read_pairs(h372_path)
}
expected = {
  'LIVE' => [22_647, 0xdf5374d4aca67677],
  'CARRIER10' => [10, 0x9926701692e6f8d4],
  'H19' => [36, 0x5c8af37cf2f002e7],
  'H294' => [21, 0xeadefa2fae582ca7],
  'H372' => [54, 0x47ab2af18f07ff59]
}
ledgers.each do |label, rows|
  raise "#{label} identity changed" unless [rows.length, fnv(rows)] == expected[label]
end
live = ledgers['LIVE'].to_h { |pair| [pair, true] }
nodes = ledgers.reject { |label, _| label == 'LIVE' }
node_sets = nodes.transform_values { |rows| rows.to_h { |pair| [pair, true] } }
node_sets.each do |label, rows|
  raise "#{label} escapes live" unless rows.keys.all? { |pair| live.key?(pair) }
end

puts 'LRC14_CARRIER10_THREE_IDEAL_PROOF_UNION_V1'
ledgers.each { |label, rows| identity(label, rows) }
labels = nodes.keys
labels.each_with_index do |left, index|
  labels[(index + 1)..].each do |right|
    overlap = (node_sets[left].keys & node_sets[right].keys).sort
    identity("OVERLAP_#{left}_#{right}", overlap)
    unless overlap.empty?
      puts "OVERLAP_ROWS_#{left}_#{right} " +
           overlap.map { |q, r| "#{q},#{r}" }.join(' ')
    end
  end
end
union = node_sets.values.flat_map(&:keys).uniq.sort
union_set = union.to_h { |pair| [pair, true] }
final = ledgers['LIVE'].reject { |pair| union_set.key?(pair) }
raise 'typed union count changed' unless union.length == 118 && final.length == 22_529
identity('UNION', union)
identity('FINAL', final)
maximum = final.map(&:last).max
top = final.select { |_, r| r == maximum }
puts "FINAL_MAX #{maximum} TOP_COUNT #{top.length} TOP " +
     top.map { |q, r| "#{q},#{r}" }.join(' ')
puts 'SCOPE FINITE_EXACT_TYPED_SET_UNION_SEPARATE_CERTIFICATE_NODES_NO_LRC14'
FileUtils.mkdir_p(out_dir)
File.binwrite(File.join(out_dir, 'carrier10_three_ideal_union.csv'), bytes(union))
File.binwrite(File.join(out_dir, 'post_carrier10_three_ideal_residual.csv'), bytes(final))
