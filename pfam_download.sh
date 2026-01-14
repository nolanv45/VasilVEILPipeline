cd phidra
mkdir pfam_database
cd pfam_database

wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

gunzip -c Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat
gunzip -c Pfam-A.hmm.gz > Pfam-A.hmm
rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz

hmmpress Pfam-A.hmm

cd VasilVEILPipeline