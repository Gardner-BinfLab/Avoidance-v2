#!/bin/bash

if [ ! -n "$1" ] || [ ! -n "$2" ]; then \
  echo "" >&2
  echo "Usage: ./Feature_scoring.sh <Threads> <CDS.fa> [option]" >&2
  echo "  <Threads>       INT    Set number of threads for RNAup_avoidance_calculator.py" >&2
  echo "  <CDS.fa>        STR    Coding sequences in fasta format" >&2
  echo "" >&2
  echo "Option:" >&2
  echo "  <5UTR_sequence> STR    Use your own 5UTR sequence (71 nt) if your plasmid backbone is not pET vector" >&2
  echo "                         Default is GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT" >&2
  echo "" >&2
  echo "Dependencies: biopython, codonw, ViennaRNA" >&2
  echo "  Install through Miniconda" >&2
  echo "    conda install -c bioconda biopython codonw viennarna" >&2
  echo "  or source/precompiled binary package" >&2
  echo "    https://sourceforge.net/projects/codonw/" >&2
  echo "    https://www.tbi.univie.ac.at/RNA/" >&2
  echo "" >&2
  echo "Examples:" >&2
  echo "  ./Feature_scoring.sh 40 sequence.fasta" >&2
  echo "  ./Feature_scoring.sh 20 sequence.fasta GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT" >&2
  exit 1
fi

THREADS=${1}
CDS=${2}
UTR5=${3:-GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT}
DIR=$PWD
SRC="$(cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)"

echo ""
echo "....................................................................................."
echo "Calculating Open Energies for Accessibility(1:30)..."
echo ""
mkdir rnaplfold
cd rnaplfold
time awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${CDS} \
| awk -v U=${UTR5} '{print $1 "\n" U $2}' \
| RNAplfold -W 210 -u 210 -O
for i in *_openen; do \
  paste \
  <(echo ${i}) \
  <(awk 'BEGIN{OFS="\t"} NR>73 && NR<=102 {print $42,$43,$44,$45,$46}' ${i} \
  | awk '{for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print}')
done \
| sed 's/_openen//;s/ /\t/g' \
| sed '1i Accession\topenen41\topenen42\topenen43\topenen44\topenen45' > ${DIR}/openen.out
cd ${DIR}

echo ""
echo "....................................................................................."
echo "Calculating Minimum Free Energies for STR(-30:+30)..."
echo ""
time awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${CDS} \
| awk '{print $1 "\n" substr($NF,1,30)}' \
| awk -v U=${UTR5} '{print $1 "\n" substr(U,length(U)-29) $2}'  \
| RNAfold --noPS \
| awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print $0}' \
| awk '{print $1 "\t" $NF}' \
| sed 's/[()]//g' \
| sed '1i Accession\tSTR'> ${DIR}/rnafold.out

echo ""
echo "....................................................................................."
echo "Calculating basepair classes..."
echo ""
python ${SRC}/basepair.py -i ${CDS} -w 60 -s 6 | gzip > basepair.raw.gz
zcat basepair.raw.gz \
| sed 's/\t/,/2;s/0,0/3/;s/0,1/2/;s/0,2/1/;s/1,0/2/;s/1,1/1/;s/1,2/3/;s/2,0/1/;s/2,1/3/;s/2,2/2/' \
| sort | uniq -c \
| awk 'BEGIN{OFS="\t"}{print $2,$3,$1}' > basepair.txt
paste \
<(awk '$2==1 {print $1 "\t" $3}' basepair.txt) \
<(awk '$2==2 {print $3}' basepair.txt) \
<(awk '$2==3 {print $3}' basepair.txt) \
| join -t$'\t' \
<(sort -k1,1 -) \
<(awk '{a[$1] += $3} END{for (i in a) print i "\t" a[i]}' basepair.txt | sort -k1,1) \
| sed '1i Accession\tc1\tc2\tc3\tTotal_c' > basepair.out

echo ""
echo "....................................................................................."
echo "Calculating Codon Usage Indices..."
time codonw ${CDS} codonw.out codonw.blk \
-cai -fop -cbi -enc -gc -gc3s -sil_base \
-nomenu -silent
awk 'BEGIN{OFS="\t"}{sub("title","Accession"); print $1,$2,$3,$4,$5,$6,$9,$10,$11}' codonw.out > codonw.out2
mv codonw.out2 codonw.out

echo ""
echo "....................................................................................."
echo "Calculating Biosynthetic Cost and tAI..."
time python ${SRC}/CodonMuSe/CodonMuSe.py \
-f ${CDS} \
-tscan ${SRC}/CodonMuSe/eschColi_BL21DE3_1-tRNAs.out \
-ind -par -m GC_Sc_St \
-tc 11 -GC 0.2075 -Sc -0.0741 -St 0.2222
mv ${CDS%.*}_OptimisationResults.txt codonmuse.out
sed -i 's/%//g' codonmuse.out

echo ""
echo "....................................................................................."
echo "Calculating Minimum Free Energies for Avoidance of mRNA:ncRNA Intereactions(1:30)..."
echo ""
mkdir ${DIR}/avoidance_scores
cd ${DIR}/avoidance_scores
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${CDS} \
| awk '{print $1 "\n" substr($NF,1,30)}' > w1.30.fa
time python ${SRC}/RNAup_avoidance_calculator.py \
-mrna w1.30.fa \
-ncrna ${SRC}/data/ncrna.ril.fa \
-parallel ${THREADS}
for i in *.result; do \
  sequence=$(echo ${i%.*})
  avoidance_mfe=$(cat ${i%.*}.result \
  | awk 'match($0,/.*\s\(([-]*[0-9]+.[0-9]+)\s=\s.*/,m){avoid+=m[1]}END{print avoid}')
  echo -e $sequence"\t"$avoidance_mfe
done \
| sed '1i Accession\tAvoidance' > ${DIR}/avoidance.out
#avoidance for each ncRNA
for i in *.result; do \
  paste \
  <(echo ${i%.*}) \
  <(awk 'match($0,/.*\s\(([-]*[0-9]+.[0-9]+)\s=\s.*/,m){print m[1]}' ${i} | transpose -t)
done \
| cat \
<(grep ">" ${SRC}/data/ncrna.ril.fa | sed 's/>//;s/-/_/g' | transpose -t | sed 's/^/Accession\t/') - \
> ${DIR}/avoidance.ncrna.out
cd ${DIR}

#rm -r avoidance_scores rnaplfold ${CDS%.*}_*.txt *.blk basepair.txt