#!/bin/bash

if [ ! -n "$1" ]; then \
  echo "" >&2
  echo "Usage: ./featureCalc.sh <CDS.fa> [options...]" >&2
  echo "Options:" >&2
  echo "  <Threads>       INT    Set number of threads for Avoidance calculator" >&2
  echo "                         Default is 40 threads" >&2
  echo "  <5UTR_sequence> STR    Use your own 5UTR sequence if your plasmid backbone is not pET vector" >&2
  echo "                         Default is GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT" >&2
  echo "" >&2
  echo "Dependencies: biopython, codonw, ViennaRNA" >&2
  echo "  Install from Miniconda: conda install -c bioconda biopython codonw viennarna" >&2
  echo "  Or from source: https://sourceforge.net/projects/codonw/" >&2
  echo "" >&2
  exit 1
fi

THREADS=${2:-40}
UTR5=${3:-GGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT}
DIR=$PWD
SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

echo ""
echo "....................................................................................."
echo "Calculating Codon Usage Indices..."
/usr/bin/time codonw ${1} codonw.out codonw.blk \
-cai -fop -cbi -enc -gc -gc3s -sil_base \
-nomenu -silent
awk 'BEGIN{OFS="\t"}{sub("title","Accession"); print $1,$2,$3,$4,$5,$6,$9,$10,$11}' codonw.out > codonw.out2
mv codonw.out2 codonw.out

echo ""
echo "....................................................................................."
echo "Calculating Biosynthetic Cost and tAI..."
/usr/bin/time python ${SRC}/CodonMuSe/CodonMuSe.py \
-f ${1} \
-tscan ${SRC}/CodonMuSe/eschColi_BL21DE3_1-tRNAs.out \
-ind -par -m GC_Sc_St \
-tc 11 -GC 0.2075 -Sc -0.0741 -St 0.2222
mv ${1%.*}_OptimisationResults.txt codonmuse.out
sed -i 's/%//g' codonmuse.out

echo ""
echo "....................................................................................."
echo "Calculating Minimum Free Energies for Avoidance of mRNA:ncRNA Intereactions(1:30)..."
echo ""
mkdir ${DIR}/avoidance_scores
cd ${DIR}/avoidance_scores
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${DIR}/${1} \
| awk '{print $1 "\n" substr($NF,1,30)}' > w1.30.fa
/usr/bin/time python ${SRC}/RNAup_avoidance_calculator.py \
-mrna w1.30.fa \
-ncrna ${SRC}/data/ncrna.ril.fa \
-parallel ${THREADS}
for i in *.result; do \
  sequence=$(echo ${i%.*})
  avoidance_mfe=$(cat ${i%.*}.result \
  | gawk 'match($0,/.*\s\(([-]*[0-9]+.[0-9]+)\s=\s.*/,m){avoid+=m[1]}END{print avoid}')
  echo -e $sequence"\t"$avoidance_mfe
done \
| sed '1i Accession\tAvoidance' > ${DIR}/avoidance.out
cd ${DIR}

echo ""
echo "....................................................................................."
echo "Calculating Open Energies for Accessibility(1:30)..."
echo ""
mkdir rnaplfold
cd rnaplfold
/usr/bin/time awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${DIR}/${1} \
| awk -v U=${UTR5} '{print $1 "\n" U $2}' \
| RNAplfold -W 210 -u 210 -O
for i in *_openen; do \
  paste \
  <(echo ${i}) \
  <(awk 'BEGIN{OFS="\t"} NR>73 && NR<=101 {print $42,$43,$44,$45,$46}' ${i} \
  | awk '{for (i=1;i<=NF;i++){sums[i]+=$i;maxi=i}} END{for(i=1;i<=maxi;i++){printf("%s ",sums[i])} print}') \
  | sed 's/_openen//'
done \
| sed '1i Accession\topenen42\topenen43\topenen44\topenen45\topenen46' > ${DIR}/openen.out
cd ${DIR}

echo ""
echo "....................................................................................."
echo "Calculating Minimum Free Energies for STR(-30:+30)..."
echo ""
/usr/bin/time awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' ${1} \
| awk '{print $1 "\n" substr($NF,1,30)}' \
| awk -v U=${UTR5} '{print $1 "\n" substr(U,length(U)-29) $2}'  \
| RNAfold --noPS \
| awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print $0}' \
| awk '{print $1 "\t" $NF}' \
| sed 's/[()]//g' \
| sed '1i Accession\tSTR'> rnafold.out

rm -r avoidance_scores rnaplfold *.txt *.blk