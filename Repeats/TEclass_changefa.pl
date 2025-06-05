#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;


# Variables to store command-line arguments
my $fasta;
my $tecres;

# Configure GetOptions
GetOptions(
    'file1=s' => \$fasta,  # --file1  expects a string
    'file2=s' => \$tecres,  # --file2  expects a string
) or die "Error in command-line arguments. Usage: $0 --file1 <FILE> --file2 <FILE>\n";

# Check if files were actually provided
unless ($fasta && $tecres) {
    die "Usage: $0 --file1 <fasta file> --file2 <tab results from TEclass>\n";
}

# Open the first file
open my $fh1, '<', $fasta
    or die "Could not open '$fasta' for reading: $!\n";

# Open the second file
open my $fh2, '<', $tecres
    or die "Could not open '$tecres' for reading: $!\n";

my %records;

while (<$fh2>) {
    next if /^\s*$/;       # Skip empty lines
    next if /^#/;          # (Optional) Skip commented lines if needed

    # Split on whitespace; adjust if columns are separated differently
    my @fields = split /\s+/, $_;

    # We only care about the first 5 columns
    my @first_five = @fields[0..4];

    # The first column is our "accession"
    my $accession = $fields[0];

    $records{$accession} = \@first_five;

}

close $fh2;

my %TEclass2RM = (
    'Copia'    => 'LTR/Copia',
    'Crypton'  => 'DNA/Crypton',
    'ERV'      => 'LTR/ERV',
    'Gypsy'    => 'LTR/Gypsy',
    'hAT'      => 'DNA/hAT',
    'Helitron' => 'RC/Helitron',
    'Jockey'   => 'LINE/I-Jockey',
    'L1_L2'    => 'LINE',
    'Maverick' => 'DNA/Maverick',
    'Merlin'   => 'DNA/Merlin',
    'P'        => 'DNA/P',
    'Pao'      => 'LTR/Pao',
    'RTE'      => 'LINE/RTE',
    'SINE'     => 'SINE',
    'TcMar'    => 'DNA/Mariner',
    'Transib'  => 'DNA/CMC-Transib',
);



# Read lines from STDIN or from file(s) passed in @ARGV
while (<$fh1>) {
	chomp;
# Only process lines that start with '>'
	if(/^>(.*)/){
  my $header = $_;
# Extract the part after '>'

# Split on the first '#' into two parts
	my ($family, $post_hash) = split /#/, $header, 2;
  $family =~ s/>//g;  # replace '>' with underscore
# Trim leading/trailing whitespace on $post_hash
	$post_hash =~ s/^\s+|\s+$//g;

# Split $post_hash into two parts:
#   1) type_subtype (the first token, no spaces)
#   2) rest         (everything else)
	my ($type_subtype, $rest) = split /\s+/, $post_hash, 2;
	$rest = "" unless defined $rest;

# Now split $type_subtype on '/' into $type and $subtype
	my ($type, $subtype) = split m{/}, $type_subtype, 2;
	$subtype ||= "";  # If there's no slash, $subtype is undefined; set it to ""

    if ( exists $records{$family} ) {
        my $aref = $records{$family};  # This is an array reference        
        my $prob = $aref->[4];
        my $teclass = $aref->[3];
        my $newRepeatMasker = $TEclass2RM{$teclass};
        print ">$family#$newRepeatMasker $rest \[Classified by TEclass with prob=$prob\]\n";

    } else {
        print ">$family#$post_hash\n";
    }
  }
   else {
  print $_;
  print"\n";
  }
  }
  close($fh1);
  close($fh2);
