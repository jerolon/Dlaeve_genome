#!/usr/bin/perl
############################################
#
# Transform the output from REvolutionH-tl
# into the table needed for Cafe
#
###########################################
use strict;
use warnings;
use Getopt::Long;

# Boilerplate
my $filename;
my $delimiter = "\t"; # Default delimiter is tab
my $output_file = 'output.txt'; # Default output file

# Parse command-line options
GetOptions(
    "file=s"     => \$filename,   # File name
    "output=s"     => \$output_file, # Output file name
) or die "Error in command-line arguments.\n";

# Ensure the file name is provided
unless ($filename) {
    die "Usage: $0 --file <filename> [--output <output_file>]\n";
}

# Open the file
open my $fh, '<', $filename or die "Cannot open file '$filename': $!\n";


# Read the first line to extract species names
my $header_line = <$fh>;
chomp $header_line;
my @header_fields = split /$delimiter/, $header_line;

# Extract species names (from the 4th column onward)
my @species_names = @header_fields[3..$#header_fields];
# Prepare the output file
open my $out_fh, '>', $output_file or die "Cannot open output file '$output_file': $!\n";

# Write the header line for the output file
print $out_fh join("\t", "Desc", "Family_ID", @species_names), "\n";

# Process each subsequent line
while (my $line = <$fh>) {
    chomp $line;  # Remove the newline character

    # Skip empty lines
    next if $line =~ /^\s*$/;

    # Split the line into fields
    my @fields = split /$delimiter/, $line;

    # Extract the fixed columns
    my $orthogroup = $fields[0];
    #my $n_genes = $fields[1];
    #my $n_species = $fields[2];

    # Extract species data
    #my %species_data;
    my %gene_counts;
    for my $i (3..$#fields) {
        my $species_name = $species_names[$i - 3];
        my $gene_list = $fields[$i] // ''; # Use an empty string if undefined
        my @genes = $gene_list eq '' ? () : split /,/, $gene_list; # Handle empty lists
       #$species_data{$species_name} = \@genes;
       $gene_counts{$species_name} = scalar @genes; # Count genes
    }

    my @counts = map { $gene_counts{$_} // 0 } @species_names; # 0 for missing species
    print $out_fh join("\t", "(null)", $orthogroup, @counts), "\n";
}

close $fh;
close $out_fh;
print "Output written to $output_file\n";
