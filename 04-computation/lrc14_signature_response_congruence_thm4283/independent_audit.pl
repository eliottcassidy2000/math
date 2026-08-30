#!/usr/bin/env perl
# Independent standard-library replay of the THM-4283 quotient-purity census.

use strict;
use warnings;
use Digest::SHA qw(sha256_hex);
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Spec;

my $repo = abs_path(File::Spec->catdir($Bin, "..", ".."));
my $packet = File::Spec->catdir(
    $repo,
    "05-knowledge",
    "results",
    "lrc14_inactive_signature_deck_surgery_thm4282"
);

my %inputs = (
    signature421 => [
        File::Spec->catfile($packet, "components", "index367", "results", "post_thm4281_full421_inactive_signatures.csv"),
        "4f0e8da3fdab8bd5a0e14f3b4fa30050602025f63486aa35e0cf03374e6e3832",
        24_223,
    ],
    surgery520_367 => [
        File::Spec->catfile($packet, "components", "surgery520", "results", "combined_gain_520_367.csv"),
        "04aed4c107b3244a5e488266c46d1cae3bfffb20cddd831fc2d474c7b8c16a0e",
        586,
    ],
    surgery256 => [
        File::Spec->catfile($packet, "components", "surgery256", "results", "surgery_common.csv"),
        "5946ff653c51a74eba09a14430c494074e53b5aba87c3159bd17bafbe9e605d5",
        188,
    ],
    carrier90 => [
        File::Spec->catfile($packet, "results", "carrier90.csv"),
        "222afb7618d887f32847b4531ffedf5616f20c2196e92f52ca2c11b09e1eab1f",
        90,
    ],
    union850 => [
        File::Spec->catfile($packet, "results", "deletion_union850.csv"),
        "7ad581bccd253e1778b972e8a303207da44534e6b995fa3ba15bd34b2801505b",
        850,
    ],
);

sub slurp_raw {
    my ($path) = @_;
    open my $handle, "<:raw", $path or die "cannot open $path: $!\n";
    local $/;
    my $bytes = <$handle>;
    close $handle or die "cannot close $path: $!\n";
    return $bytes;
}

for my $name (sort keys %inputs) {
    my ($path, $expected_hash) = @{$inputs{$name}};
    die "missing input $path\n" unless -f $path;
    my $actual_hash = sha256_hex(slurp_raw($path));
    die "$name SHA-256 changed: $actual_hash\n"
        unless $actual_hash eq $expected_hash;
}

my %nibble_popcount = (
    0 => 0, 1 => 1, 2 => 1, 3 => 2,
    4 => 1, 5 => 2, 6 => 2, 7 => 3,
    8 => 1, 9 => 2, a => 2, b => 3,
    c => 2, d => 3, e => 3, f => 4,
);

sub hex_popcount {
    my ($hex) = @_;
    my $total = 0;
    $total += $nibble_popcount{$_} for split //, lc $hex;
    return $total;
}

my %fibres;
my %signature_of;
my ($signature_path, undef, $signature_count) = @{$inputs{signature421}};
open my $signature_handle, "<", $signature_path
    or die "cannot open $signature_path: $!\n";
my $header = <$signature_handle>;
chomp $header;
die "signature header changed\n"
    unless $header eq "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6";
my $line_number = 1;
while (my $line = <$signature_handle>) {
    ++$line_number;
    chomp $line;
    my @fields = split /,/, $line, -1;
    die "bad signature row $line_number\n" unless @fields == 10;
    my ($q, $r, $inactive_count, @words) = @fields;
    die "bad pair on signature row $line_number\n"
        unless $q =~ /^\d+$/ && $r =~ /^\d+$/ && 0 < $q && $q < $r;
    my $pair = "$q,$r";
    die "duplicate signature pair $pair\n" if exists $signature_of{$pair};
    for my $word (@words) {
        die "bad signature word at $pair\n"
            unless $word =~ /^[0-9a-f]{16}$/;
    }
    die "signature uses an index above 420 at $pair\n"
        unless $words[6] =~ /^000000[01][0-9a-f]{9}$/;
    my $observed_count = 0;
    $observed_count += hex_popcount($_) for @words;
    die "inactive-count mismatch at $pair\n"
        unless $inactive_count =~ /^\d+$/ && $observed_count == $inactive_count;
    my $signature = join ":", @words;
    $signature_of{$pair} = $signature;
    push @{$fibres{$signature}}, $pair;
}
close $signature_handle or die "cannot close $signature_path: $!\n";
die "signature count changed\n"
    unless scalar(keys %signature_of) == $signature_count;
for my $pairs (values %fibres) {
    @{$pairs} = sort {
        my ($aq, $ar) = split /,/, $a;
        my ($bq, $br) = split /,/, $b;
        $aq <=> $bq || $ar <=> $br;
    } @{$pairs};
}

my @target_names = qw(surgery520_367 surgery256 carrier90 union850);
my %targets;
for my $name (@target_names) {
    my ($path, undef, $expected_count) = @{$inputs{$name}};
    open my $handle, "<", $path or die "cannot open $path: $!\n";
    my %pairs;
    my $row_number = 0;
    while (my $line = <$handle>) {
        ++$row_number;
        chomp $line;
        die "bad $name row $row_number\n" unless $line =~ /^(\d+),(\d+)$/;
        my $pair = "$1,$2";
        die "duplicate $name pair $pair\n" if exists $pairs{$pair};
        die "$name pair outside signature universe $pair\n"
            unless exists $signature_of{$pair};
        $pairs{$pair} = 1;
    }
    close $handle or die "cannot close $path: $!\n";
    die "$name count changed\n" unless scalar(keys %pairs) == $expected_count;
    $targets{$name} = \%pairs;
}

my %typed_union = (
    %{$targets{surgery520_367}},
    %{$targets{surgery256}},
    %{$targets{carrier90}},
);
die "typed union changed\n"
    unless join("\n", sort keys %typed_union) eq
           join("\n", sort keys %{$targets{union850}});

my %expected = (
    surgery520_367 => "17593,6,5,770,555,215",
    surgery256 => "17587,15,2,367,186,181",
    carrier90 => "17522,34,48,1460,42,1418",
    union850 => "17502,48,54,2351,770,1581",
);

print "THM-4283 INDEPENDENT PERL REPLAY\n";
my @nonsingletons = grep { scalar(@{$_}) > 1 } values %fibres;
my $nonsingleton_rows = 0;
$nonsingleton_rows += scalar(@{$_}) for @nonsingletons;
die "base census changed\n"
    unless scalar(keys %fibres) == 17_604 &&
           scalar(@nonsingletons) == 803 &&
           $nonsingleton_rows == 7_422;
print "BASE ROWS 24223 FIBRES 17604 NONSINGLETON_FIBRES 803 NONSINGLETON_ROWS 7422\n";

for my $name (@target_names) {
    my ($negative, $mixed, $positive) = (0, 0, 0);
    my ($mixed_rows, $mixed_positive, $mixed_negative) = (0, 0, 0);
    for my $pairs (values %fibres) {
        my $hits = 0;
        $hits += exists $targets{$name}{$_} ? 1 : 0 for @{$pairs};
        if ($hits == 0) {
            ++$negative;
        } elsif ($hits == scalar(@{$pairs})) {
            ++$positive;
        } else {
            ++$mixed;
            $mixed_rows += scalar(@{$pairs});
            $mixed_positive += $hits;
            $mixed_negative += scalar(@{$pairs}) - $hits;
        }
    }
    my $summary = join ",", $negative, $mixed, $positive,
        $mixed_rows, $mixed_positive, $mixed_negative;
    die "$name summary changed: $summary\n" unless $summary eq $expected{$name};
    print "TARGET $name ALL_NEGATIVE $negative MIXED $mixed ALL_POSITIVE $positive ",
          "MIXED_ROWS $mixed_rows POSITIVE_ROWS $mixed_positive NEGATIVE_ROWS $mixed_negative\n";
}

my %word_histogram;
my ($four_signature, $four_words);
for my $signature (keys %fibres) {
    my %words;
    for my $pair (@{$fibres{$signature}}) {
        my $word = join "", map { exists $targets{$_}{$pair} ? 1 : 0 } @target_names;
        ++$words{$word};
    }
    ++$word_histogram{scalar(keys %words)};
    ($four_signature, $four_words) = ($signature, \%words)
        if scalar(keys %words) == 4;
}
die "joint response-word histogram changed\n"
    unless join(",", map { "$word_histogram{$_}" } 1..4) eq "17556,42,5,1";
die "four-word signature changed\n"
    unless defined($four_signature) &&
           join(":", split /:/, $four_signature) eq
           "0000000000000000:0000080000000000:0000000000000000:" .
           "0000000000000000:0000000000000000:0000000000000000:" .
           "0000000000000000";
my $four_summary = join ",", map { "$_:$four_words->{$_}" } sort keys %{$four_words};
die "four-word census changed: $four_summary\n"
    unless $four_summary eq "0000:23,0101:9,1001:21,1101:9";
print "JOINT_WORD_HISTOGRAM 1:17556 2:42 3:5 4:1\n";
print "FOUR_WORD_SIGNATURE 107 0000:23 0101:9 1001:21 1101:9\n";
print "VERDICT PASS INDEPENDENT_IMPLEMENTATION_AGREES\n";
