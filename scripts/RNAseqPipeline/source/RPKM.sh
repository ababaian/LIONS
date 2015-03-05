#!/bin/sh
set -e
set -o pipefail

# RPKM.sh

usage()
{
echo "Script to calculate RPKM (total, average, max, min)"
echo ""
echo "USAGE:"
echo "/home/jlever/bin/RPKM.sh <name> <read length> <type>[G/T, where G is for genes' and T is for 'transcripts'] <metric>[A/M/X, where A=average exon coverage, M=median exon coverage, X=95% value of exon coverage] <species>[mm9/hg18/hg19] <strand specific S/R>"
echo ""
echo "OUTPUT FORMAT:"
echo "exn/tot (exonic,total) files: id rpkm1 rpkm2 z12 a12 [tab separated]  (z- z-score, a-asymmetry)"
echo "ave (average) files: id min_rpkm1 ave_rpkm1 max_rpkm1 min_rpkm2 ave_rpkm2 max_rpkm2 min_z12 ave_z12 maz_z12 min_a12 ave_a12 max_a12 [tab separated]  (z- z-score, a-asymmetry)"
}

if [ "$#" -ne 6 ] ; then
        usage
        exit 1
fi

date
echo "Running - $0"
svn info https://svn01.bcgsc.ca/svn/Solexa_Shell/RPKM.sh | awk '/Revision/ || /Last Changed Date/'

# read lengths
name=$1
length=$2

type=$3
metric=$4

if [ "$metric" == "A" ]; then
	cut=1-6
elif [ "$metric" == "M" ]; then
	echo "Not ready"; exit
	cut=1-5,10
elif [ "$metric" == "X" ]; then
	echo "Not ready"; exit
	cut=1-5,11
else
	echo "Wrong value for 'metric'-parameter"
	usage
	exit 1
fi

species=$5
echo "Species $species"

sr=$6

SHELL_BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
reschrs=$(sh $SHELL_BASE/RNAgetRes.sh $species)
res=$(echo $reschrs | cut -f 1 -d ',')
chrs=$(echo $reschrs | cut -f 2 -d ',')

pcid=$res$species"_genes.pc.EnsID"
ncid=$res$species"_genes.nc.EnsID"

echo "Library: $name"
rm -rf $name".rpkm.log"

# get total reads in gene
# We need another file: genes for coverage
#less *_genes.$name*coverage* | awk '{l=$3-$2+1; r=r+$6*l} END{print "Total number of reads in genes (entire gene boundary): "int(r/"'$length'"+0.5)}'

# Note in case of strand specific there are two files; some of the coverage exons will overlap, but when the leakage is low it should not matter
less *_exons_for_coverage*$name*coverage* | awk '{print $1":"$2":"$3"\t"$4 }' | sed 's/chr//' | sort -k1,1 > f

# get total reads in all exons and exon/intron ratio
reads_in_genes=$(less *_genes_for_coverage*$name*coverage* | awk '{t=t+($3-$2+1)*$4} END {print int(t/("'$length'"+0)+0.5)}');
reads_in_exons=$(less f | sed 's/:/\t/g' | awk '{l=$3-$2+1; r=r+$4*l} END {print int(r/("'$length'"+0)+0.5)}');

echo "Total number of reads within gene boundaries:	$reads_in_genes" >> $name".rpkm.log"
echo "Total number of reads in all (collapsed) exons:	$reads_in_exons" >> $name".rpkm.log"
reads_in_introns=`expr $reads_in_genes - $reads_in_exons`
echo "Total number of reads in all introns:	$reads_in_introns" >> $name".rpkm.log"
ieratio=$(echo "$reads_in_introns" | awk '{print $0/("'$reads_in_exons'"+0)}') 
echo "Intron/Exon ratio:	$ieratio" >> $name".rpkm.log"

EXPRESSION_THRESHOLD=0.005

if [ "$sr" == "S" ]; then
	echo "Strand specific RNA-seq"
	join f $res$species"_exons_for_coverage.ss.pc" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (protein coding, no MT, no ribo proteins, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.ss.nc" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (non coding, no MT, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.ss.rb" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (ribosomal proteins, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.ss.mt" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (MT, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.ss.pc" > fpc
else
	echo "Regular RNA-seq"
	join f $res$species"_exons_for_coverage.pc" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (protein coding, no MT, no ribo proteins, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.nc" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (non coding, no MT, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.rb" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (ribosomal proteins, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.mt" | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; r=r+$4*l} END{print "Total number of exonic reads (MT, "c" exons):\t"int(r/("'$length'"+0)+0.5)}' >> $name".rpkm.log"
	join f $res$species"_exons_for_coverage.pc" > fpc
fi
echo $res$species"_exons_for_coverage.ss.pc"

count=$(less fpc | wc -l)
echo "Number of protein coding exons for coverage: $count"

th=$(echo "$count*$EXPRESSION_THRESHOLD" | bc | awk '{print int($0)}');
echo "Number of top $EXPRESSION_THRESHOLD expressed exons to exclude: $th"

reads=$(less fpc | sort -k2,2nr | sed 's/:/\t/g' | awk '{c=c+1; l=$3-$2+1; if(c>("'$th'"+0)){r=r+$4*l}} END{print int(r/("'$length'"+0)+0.5)}')
echo "Total number exonic reads (protein coding only; excluded: 1) ribosomal proteins genes excluded, 2) MT reads, 3) top expressed $EXPRESSION_THRESHOLD fraction of exons): $reads"
echo "NOTE: this number is going to be used in normalization for the RPKM calculations"
echo "Total number of exonic reads for RPKM (protein coding; no MT, no ribo proteins, top expressed $EXPRESSION_THRESHOLD exons excluded):	$reads" >> $name".rpkm.log"

if [ "$type" == "G" ]; then
	prefix=genes;
elif [ "$type" == "T" ]; then
	prefix=transcripts;
else
	echo "Wrong value for the type: G/T"
	usage
	exit 1
fi

cat *_exons_for_$prefix*$name*coverage* | cut -f$cut | awk '{print $1":"$2"-"$3"<"$4"|"$5"\t"$3-$2+1"\t"$6}' | sort -k1,1 > f

echo "reads=$reads"
echo "length=$length"

echo "less f | 
awk 'BEGIN{n=1000000000/(\"'$reads'\"+0)/(\"'$length'\"+0)} { print $1\"\t\"$2*$3\"\t\"$3*n}' | sed 's/|/\t/' | sort -k2,2 -k1,1"

########
less f | 
awk 'BEGIN{n=1000000000/("'$reads'"+0)/("'$length'"+0)} { print $1"\t"$2*$3"\t"$3*n}' | sed 's/|/\t/' | sort -k2,2 -k1,1 > $name.$type.exn.$metric.rpkm

echo "pcid=$pcid"
echo "ncid=$ncid"

########
echo "Calculating $prefix RPKM (total)..."
less f | cut -f2- -d"|" | sort -k1,1 | 
awk 'BEGIN{n=1000000000/("'$reads'"+0)/("'$length'"+0);nr=1/("'$length'"+0);rmin=1000000000000000; rmax=0;} {if($1==id) {t=t+$3*$2; len=len+$2; c=c+1; rave=rave+$3; if($3>rmax) {rmax=$3}; if($3<rmin) {rmin=$3};} else {if(id!=null) {ir=int(t*nr+0.5); r=t/len*n; print id"\t"ir"\t"r"\t"rave/c*n"\t"rmin*n"\t"rmax*n}; t=$3*$2; len=$2; rave=$3; rmax=$3; rmin=$3; c=1; id=$1}} END{ir=int(t*nr+0.5); r=t/len*n; print id"\t"ir"\t"r"\t"rave/c*n"\t"rmin*n"\t"rmax*n}' | sort -k1,1 > tot
join $pcid tot | sed 's/ /\t/g' > $name.$type.$metric.rpkm.pc
join $ncid tot | sed 's/ /\t/g' > $name.$type.$metric.rpkm.nc

# total RPKM distribution summary
echo ""
echo "less $name.$type.$metric.rpkm.pc"

less $name.$type.$metric.rpkm.pc | awk 'BEGIN{l="Total protein coding:";	l1="RPKM>1:"; l5="RPKM>5:"; l10="RPKM>10:"; l50="RPKM>50:";	l100="RPKM>100:"}{c0=c0+1; if($3>1) {c1=c1+1}; if($3>5){c5=c5+1}; if($3>10){c10=c10+1}; if($3>50){c50=c50+1}; if($3>100){c100=c100+1}} END{f1=int(c1/c0*10000+0.5)/100; f5=int(c5/c0*10000+0.5)/100; f10=int(c10/c0*10000+0.5)/100; f50=int(c50/c0*10000+0.5)/100; f100=int(c100/c0*10000+0.5)/100; print l"\t"c0"\n"l1"\t"c1"("f1"%)\n"l5"\t"c5"("f5"%)\n"l10"\t"c10"("f10"%)\n"l50"\t"c50"("f50"%)\n"l100"\t"c100"("f100"%)"}' >> $name".rpkm.log"

echo "less $name.$type.$metric.rpkm.nc"
less $name.$type.$metric.rpkm.nc | awk 'BEGIN{l="Total non coding:";	l1="RPKM>1:"; l5="RPKM>5:"; l10="RPKM>10:"; l50="RPKM>50:";	l100="RPKM>100:"}{c0=c0+1; if($3>1) {c1=c1+1}; if($3>5){c5=c5+1}; if($3>10){c10=c10+1}; if($3>50){c50=c50+1}; if($3>100){c100=c100+1}} END{f1=int(c1/c0*10000+0.5)/100; f5=int(c5/c0*10000+0.5)/100; f10=int(c10/c0*10000+0.5)/100; f50=int(c50/c0*10000+0.5)/100; f100=int(c100/c0*10000+0.5)/100; if(c0==0){c0=1}; print l"\t"c0"\n"l1"\t"c1"("f1"%)\n"l5"\t"c5"("f5"%)\n"l10"\t"c10"("f10"%)\n"l50"\t"c50"("f50"%)\n"l100"\t"c100"("f100"%)"}' >> $name".rpkm.log"

#rm -rf tot ave f fpc 
