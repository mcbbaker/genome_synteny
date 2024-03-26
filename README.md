# genome-synteny
A workflow for obtaining syntenic blocks from annotated genome assemblies
## Data Needed
- Fasta file containing gene peptide sequences for each genome

## 1) All-By-All Blast
Create and execute shell script that cats all fasta files into one file:
> ./cmd.sh

Create a blast database from these sequences:
> makeblastdb -in all.ahrd.pass.gl.pep.fasta -dbtype prot -title all.prot.db -parse_seqids -out all.prot.db -logfile all.prot.db.log 2> all.prot.db.err

Blast the sequences against the database:
> blastp -query all.ahrd.pass.gl.pep.fasta -db all.prot.db -out all-by-all-hits.blast -outfmt 6 -num_threads 60 -evalue 1e-20 2> blastp.err

## 2) 
