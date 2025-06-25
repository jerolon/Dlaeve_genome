#!/usr/bin/bash

for alignment in *.fasta; do
    gene=$(basename "$alignment" .fasta)

#get 10 character names of sequences    
seq_names=$(grep ">" "$alignment" | cut -c2-11 | sed 's/ *$//')

# Create a temporary control file
ctl_file="${gene}_codeml.ctl"
cat > "$ctl_file" << EOF
seqfile = $alignment
treefile = tree.nwk
outfile = ${gene}_results.txt

noisy = 3
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
model = 0
NSsites = 0
icode = 4 
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.1
EOF

tree="("
for name in $seq_names; do
    tree+="$name,"
done

tree="${tree%,});"
echo "$tree" > tree.nwk

    # Run codeml
    codeml "$ctl_file"

    echo "Processed $gene"
done
