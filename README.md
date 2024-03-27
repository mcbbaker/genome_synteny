# genome-synteny
A workflow for obtaining syntenic blocks from annotated genome assemblies.
## Data Needed
- Fasta file containing peptide sequences for annotated genes 
- Annotated genes in gff3 format
- (_Optional_) Text file containing blacklisted descriptions for genes (e.g. "retrotran"), one description on each line

## 1) All-By-All Blast
Create and execute shell script that combines all fasta files into one file:
> ./cmd.sh

Create a blast database from these sequences:
> makeblastdb -in all.ahrd.pass.gl.pep.fasta -dbtype prot -title all.prot.db -parse_seqids -out all.prot.db -logfile all.prot.db.log 2> all.prot.db.err

Blast the sequences against the database:
> blastp -query all.ahrd.pass.gl.pep.fasta -db all.prot.db -out all-by-all-hits.blast -outfmt 6 -num_threads 60 -evalue 1e-20 2> blastp.err

## 2) Remove Blacklisted Genes
Use if there are certain types of genes that you want to be excluded from the syntenic blocks. For example, we keep a text file with common functional annotation descriptions for repeat-related genes that we want to be excluded. If choosing to find blacklisted genes based on descriptions, must have the functional annotation for the genes. 

Get the genes by performing a grep against the functional annotation for all genomes:
> grep -if freq-blacklist.txt all-ahrd-genes.tsv > genes-remove-description

> less genes-remove-description | cut -f 1 > genes-remove-id

Remove blast hits if they contain any blacklisted gene as either the subject or query:
> python3 remove-blast-hits.py all-by-all-hits.blast genes-remove-id > removed-all-by-all-hits.blast

## 3) Prepare .blast and .gff for MCScanX
Remove blast hits where the subject or query contains a non-primary transcript of the gene (not .1):
> python3 strip-primary.py removed-all-by-all-hits.blast > strip-removed-all-by-all-hits.blast

Get the genes for each genome:
> less Lod.1TUR.ahrd.pass.gl.gff3 | grep gene > Lod-genes

Format the genes to fit the MCScanX criteria: 
> python3 prep-gff.py Lod-genes lo > all.gff

> python3 prep-gff.py Ler-genes le >> all.gff

Make sure the prefix of .blast and .gff are the same:
> mv strip-removed-all-by-all-hits.blast all.blast

## 4) Run MCScanX
Run MCScanX with default parameters:
> MCScanX ./all > MCScanX.log 2> MCScanX.err



