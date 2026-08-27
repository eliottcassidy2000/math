#!/usr/bin/env ruby
# Independent Ruby replay of the Python overlap audit.

require 'digest'
require 'set'

def need(condition, message)
  raise message unless condition
end

def read_edges(path)
  raw = File.binread(path)
  need(raw.empty? || raw.end_with?("\n"), "missing terminal LF: #{path}")
  need(!raw.include?("\r"), "non-LF newline: #{path}")
  rows = raw.lines.map do |line|
    need(line.match?(/\A[0-9]+,[0-9]+\n\z/), "malformed row: #{line.inspect}")
    edge = line.chomp.split(',').map(&:to_i)
    need(edge[0] > 0 && edge[0] < edge[1], "invalid edge: #{edge.inspect}")
    edge.freeze
  end
  need(rows == rows.uniq.sort, "ledger not lexicographic/distinct: #{path}")
  rows
end

def serialize(rows)
  rows.map { |q, r| "#{q},#{r}\n" }.join
end

def fnv(rows)
  state = 0xcbf29ce484222325
  rows.each do |q, r|
    [q, r].pack('Q<Q<').bytes.each do |byte|
      state ^= byte
      state = (state * 0x100000001b3) & 0xffffffffffffffff
    end
  end
  format('%016x', state)
end

def identity(rows)
  [rows.length, fnv(rows), Digest::SHA256.hexdigest(serialize(rows))]
end

need(ARGV.length == 4,
     'usage: audit_overlap_independent.rb SURGERY GAIN PRIOR_RESIDUAL PYTHON_OUT')
surgery_rows = read_edges(ARGV[0])
gain_rows = read_edges(ARGV[1])
residual_rows = read_edges(ARGV[2])
python_out = ARGV[3]

surgery = surgery_rows.to_set
gain = gain_rows.to_set
residual = residual_rows.to_set
universe = gain | residual
carrier = residual.select { |_, r| r >= 645 }.to_set

sets = {
  'carrier90.csv' => carrier,
  'overlap_surgery_gain9.csv' => surgery & gain,
  'overlap_surgery_carrier5.csv' => surgery & carrier,
  'overlap_gain_carrier0.csv' => gain & carrier,
  'overlap_triple0.csv' => surgery & gain & carrier,
  'exclusive_surgery174.csv' => surgery - gain - carrier,
  'exclusive_gain577.csv' => gain - surgery - carrier,
  'exclusive_carrier85.csv' => carrier - surgery - gain,
  'combined_union850.csv' => surgery | gain | carrier,
  'final_residual23373.csv' => universe - surgery - gain - carrier
}

need((gain & residual).empty?, 'gain/residual overlap')
need(identity(universe.to_a.sort) ==
     [24223, '80ec0687d8c7dba7',
      '75a3c7616c982538363c7801ed2dab3fe9aa775ab601f7a7119dd9fb5d301552'],
     'reconstructed universe changed')

sets.each do |name, edge_set|
  rows = edge_set.to_a.sort
  raw = serialize(rows)
  need(File.binread(File.join(python_out, name)) == raw,
       "Ruby/Python ledger mismatch: #{name}")
  puts "#{name} COUNT_FNV_SHA #{identity(rows).join(' ')}"
end

final_rows = sets['final_residual23373.csv'].to_a.sort
maximum = final_rows.map(&:last).max
top = final_rows.select { |_, r| r == maximum }
need(maximum == 644, 'final maximum changed')
need(top == [[220, 644], [256, 644], [258, 644], [294, 644],
             [366, 644], [416, 644], [512, 644]], 'final top changed')
puts "FINAL_MAX #{maximum} TOP #{top.map { |q, r| "#{q},#{r}" }.join(';')}"
puts 'VERDICT PASS INDEPENDENT_RUBY_BYTE_IDENTITY'
