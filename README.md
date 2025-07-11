# genome_synteny

A workflow for generating syntenic blocks from annotated genome assemblies using MCScanX.
## Data Needed
For each genome:
- FASTA file of protein sequences for annotated genes 
- Annotated genes in GFF3 format

_Optional_, if there is a blacklist: 
- Text file containing blacklisted terms (e.g. "retrotran"), one per line
- Functional annotation of the genes

## 1) All-By-All Blast
Create and execute a shell script that concatenates all protein FASTA files:
```
./cmd.sh
```

Create a BLAST database from these sequences:
```
makeblastdb -in all.ahrd.pass.gl.pep.fasta -dbtype prot -title all.prot.db -parse_seqids -out all.prot.db -logfile all.prot.db.log 2> all.prot.db.err
```

Blast the sequences against the database:
```
blastp -query all.ahrd.pass.gl.pep.fasta -db all.prot.db -out all-by-all-hits.blast -outfmt 6 -num_threads 60 -evalue 1e-20 2> blastp.err
```

## 2) Remove Blacklisted Genes
Use grep with the blacklist terms as patterns to extract blacklisted genes from the functional annotation:
```
grep -if freq-blacklist.txt all-ahrd-genes.tsv > genes-remove-description
```
Get the ID of blacklisted genes:
```
less genes-remove-description | cut -f 1 > genes-remove-id
```
Remove BLAST hits if they contain any blacklisted gene as either the subject or query:
```
python3 remove_blacklist_hits.py all-by-all-hits.blast genes-remove-id > removed-all-by-all-hits.blast
```

## 3) Prepare .blast and .gff for MCScanX
Remove BLAST hits where the subject or query contains a non-primary transcript of the gene (not .1). This script also strips the transcript from the subject and query to conform with MCScanX criteria:
```
python3 remove_nonprimary_hits.py removed-all-by-all-hits.blast > strip-removed-all-by-all-hits.blast
```

Extract the genes from the GFF3 for each genome:
```
less Lod.1TUR.ahrd.pass.gl.gff3 | grep gene > Lod-genes
```
For each genome, format the genes to fit the MCScanX criteria and append to an 'all.gff' file: 
```
python3 prep_gff.py Lod-genes lo > all.gff
python3 prep_gff.py Ler-genes le >> all.gff
```

Make sure the prefix of .blast and .gff are the same:
```
mv strip-removed-all-by-all-hits.blast all.blast
```

## 4) Run MCScanX
Run MCScanX with default parameters:
```
MCScanX ./all > MCScanX.log 2> MCScanX.err
```

## 5) Filter
Filter the blocks based on score thresholds for intra-chromosomal and inter-chromosomal blocks:
```
python3 intra_inter_filter.py all.collinearity 1000 5000 > 1k5k-all.collinearity
```

