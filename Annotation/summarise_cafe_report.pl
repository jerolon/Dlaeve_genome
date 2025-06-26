#!/usr/bin/perl

 use strict;

 # tell this scrit which report file to use
 #summarise_cafe_report.pl cafe.cafereport
 my $report=$ARGV[0];

 open(IN,$report) || die $!;
 my @lines=<IN>;
 close(IN);
 	

 	my @node_ids;
 	my %id_to_pvalue;
 	my $nodetree='';
 	my @pvalues;
 	my ($id, $fam_newick, $family_pvalue, $elements_pvalues);
# 	my $expansions='';
# 	my $contractions=''; 	
my $newick='';
my $isheader = 1;
my $nodeorder;
my @elements;
foreach my $line (@lines) {
	#Parse info in header	
  if($isheader == 1){
  			if($line=~/\(node ID, node ID\): (.+)/) { 
			  $nodeorder=$1;
			$nodeorder=~s/\s+//g;
      @node_ids = $nodeorder =~ /(\d+)/g;
      #print "node order: $line\n";
      #print join(",", @node_ids);
			}
			if($line=~/^# IDs of nodes:(\(\S+)$/) { $newick=$1;
			  $newick=~s/\s+//g;
        $nodetree = parse_newick($newick);
#traverse_tree($nodetree);
			}

#		  if($line=~/^Expansion :\s+(.+)/) { $expansions=$1; }

#			if($line=~/^nDecrease :\s+(.+)/) { $contractions=$1;}
			#Start of actual data
			if($line=~/Newick/){$isheader=0;}
		}
		#data rows
			else{
			  ($id, $fam_newick, $family_pvalue, $elements_pvalues) = split /\t/, $line;

			  # Parse ElementsPvalues into a structured format
			  $elements_pvalues=~s/\s+//g;
        @pvalues = $elements_pvalues =~ /([\d.]+(?:e[+-]?\d+)?)/gi;
        
        
        for (my $i = 0; $i < @node_ids; $i++) {
          $id_to_pvalue{$node_ids[$i]} = scientific_to_decimal($pvalues[$i]);
          }
       
        if($id_to_pvalue{1} != "undefined"){
#      print "Parsed pvalues: $elements_pvalues\n";
#      print "Node order: $nodeorder \n";
      #print Dumper($nodetree);
      #print Dumper($fam_tree);
      $fam_newick =~ s/\s+//g; # Clean whitespace
      my $fam_tree = parse_newick_family($fam_newick);
      
      
      combine_trees($nodetree, $fam_tree, \%id_to_pvalue);
      calculate_changes($nodetree, undef);
      
      #print "Combined Tree for ortholog family $id:\n";
      #traverse_tree($nodetree);
      print_tab_tree($nodetree, 0, $id);
        } 
        #else{print "pvalues are undefined\n";}
  }
}	



# Parse the Newick string into a tree data structure


# Print the tree structure (optional)
use Data::Dumper;
#print Dumper($tree);

# Recursive function to parse the newick string with ids
sub parse_newick {
    my ($string) = @_;
    
    # Base case: If the string is a leaf (e.g., "species<id>"), return it as a hash
    if ($string =~ /^([a-zA-Z0-9]+)<(\d+)>$/) {
        return { name => $1, id => $2, children => [] };
    }

    # Recursive case: Parse branches enclosed in parentheses
if ($string =~ /^\((.+)\)<(\d+)>$/) {
         my ($content, $id) = ($1, $2);

my @children = split_children($content);

        # Recursively parse each child and return the node
        return {
            id       => $id,
            children => [ map { parse_newick($_) } @children ],
        };
    }

    # Handle invalid cases
    die "Invalid Newick format: $string\n";
}

# Example traversal function
sub traverse_tree {
    my ($node, $depth) = @_;
    $depth //= 0;

    # Print the node information
    print " " x (2 * $depth) . "Node ID: $node->{id}";
    print " (Leaf: $node->{name})" if exists $node->{name};
    print " (Ngenes: $node->{genes}) " if exists $node->{genes};
    print " (Change: $node->{change}) " if exists $node->{change};
    print "pvalue : $node->{pvalue}" if exists $node->{pvalue};
    print "\n";

    # Recursively traverse children
    foreach my $child (@{$node->{children}}) {
        traverse_tree($child, $depth + 1);
    }
}

# Traverse and print the tree


# Recursive function to parse a Newick string
sub parse_newick_family {
    my ($string) = @_;

    # Base case: If the string is a leaf, return it
if ($string =~ /^([a-zA-Z0-9_]+):([\d.]+)$/) {
 # print "Matched leaf: $1 with length $2\n";  # Debug statement
  my ($name, $genes) = $1 =~ /^(.*)_(\d+)$/; # Extract species and genes
  return { name => $name, genes => $genes, length => $2, children => []};
  #print $name . "has $genes genes\n";
    }

    # Recursive case: Parse branches enclosed in parentheses
    if ($string =~ /^\((.+)\)_([\d.]+):([\d.]+)$/) {
#print "Matched internal node: $string\n";  # Debug statement
my ($content, $branch_genes, $length) = ($1, $2, $3);
# print "Content: $content, Branch Genes: $branch_genes, Length: $length\n";  # Debug statement
 
        # Split the content of parentheses into child nodes
        my @children = split_children($content);

        # Recursively parse each child and return the node
        return {
            length   => $length,
            genes    => $branch_genes,
            children => [ map { parse_newick_family($_) } @children ],
        };
    }

    # Handle invalid cases
   die "Invalid Newick format: $string\n";
}


sub split_children {
    my ($content) = @_;
    my @children;
    my $depth = 0;
    my $buffer = '';
    foreach my $char (split //, $content) {
        if ($char eq ',' && $depth == 0) {
            push @children, $buffer;
            $buffer = '';
        } else {
            $buffer .= $char;
            $depth++ if $char eq '(';
            $depth-- if $char eq ')';
        }
    }
    push @children, $buffer if $buffer;
    return @children;
}


sub combine_trees {
    my ($id_tree, $gene_tree, $pvalue_map) = @_;
    
    
        # Check if children are array references
    if (ref($id_tree->{children}) ne 'ARRAY' || ref($gene_tree->{children}) ne 'ARRAY') {
        print "Tree structure mismatch at node:\n";
        print "ID Tree: " . Dumper($id_tree) . "\n";
        print "Gene Tree: " . Dumper($gene_tree) . "\n";
        die "Tree structure mismatch!";
    }
    # Ensure both trees are structurally consistent
    die "Tree structure mismatch!" unless scalar(@{$id_tree->{children}}) == scalar(@{$gene_tree->{children}});

    # Combine information for the current node
    $id_tree->{genes}        = $gene_tree->{genes}        if exists $gene_tree->{genes};
        # Assign the pvalue based on the node's ID
    if (exists $pvalue_map->{$id_tree->{id}}) {
        $id_tree->{pvalue} = $pvalue_map->{$id_tree->{id}};
    } else {
        $id_tree->{pvalue} = "undefined";
    }
     # Recursively combine children
    for (my $i = 0; $i < scalar(@{$id_tree->{children}}); $i++) {
        combine_trees($id_tree->{children}[$i], $gene_tree->{children}[$i], $pvalue_map);
    }
}

sub scientific_to_decimal {
    my ($value) = @_;
   return "undefined" if $value eq "undefined" || $value eq "-";
    return sprintf("%.10f", $value);
}

sub calculate_changes {
    my ($node, $parent_genes) = @_;

    # Calculate "change" for the current node
    if (defined $parent_genes) {
        $node->{change} = $node->{genes} - $parent_genes;
    } else {
        $node->{change} = undef; # Root node has no parent, so no change
    }

    # Recursively calculate for children
    foreach my $child (@{$node->{children}}) {
        calculate_changes($child, $node->{genes});
    }
}


sub print_tab_tree {
    my ($node, $depth, $family) = @_;
    $depth //= 0;

    # Print the node information
    print "$family\t";
    if(exists $node->{name}){
      print "$node->{name}\t";
    } else {
    print "$node->{id}\t";
    }
    print "$node->{genes}\t";
    if($node->{change} != "undefined"){
    print "$node->{change}\t";} else {print "0\t";}
    print "$node->{pvalue}";
    print "\n";

    # Recursively traverse children
    foreach my $child (@{$node->{children}}) {
        print_tab_tree($child, $depth + 1, $family);
    }
}
