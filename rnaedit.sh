#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64GB
#SBATCH --account=TG-MED210017
#SBATCH --partition=shared
#SBATCH -t 47:00:00
#SBATCH -o /home/a1mark/logs/%j.%N.out
#SBATCH -e /home/a1mark/logs/%j.%N.err

export PATH=/expanse/projects/qstore/data/bcbio/anaconda/bin:/expanse/projects/qstore/data/bcbio/tools/bin:/home/a1mark/.conda/envs/adam/bin:$PATH
echo "Sample: $sample"

# Environment
genome=/expanse/projects/qstore/data/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
SNPiR=/expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/software/SNPiR
DB=/expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/db
RepeatMasker=$DB/RepeatMasker.hg38.chr.bed
AluRegions=$DB/Alu.hg38.chr.bed
REDIportal=$DB/TABLE1_hg38.txt
REDIportal_bed=$DB/TABLE1_hg38.chr.bed
gene_annotation=$DB/SNPiR_annotation.chr.hg38
vep_cache=/expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/annotation/vep/GRCh38
workspace=/expanse/lustre/projects/csd691/kfisch/2022-01_Scripps_T-ALL_OC_seq/final/$sample
mkdir -p $workspace/logs

if [ ! $workspace/$workspace/"$sample"-gatk-haplotype.vcf.gz.tbi ]; then
	tabix -p vcf $workspace/"$sample"-gatk-haplotype.vcf.gz
fi

if [ ! $workspace/"$sample"-ready.bam.bai ]; then
	samtools index "$sample"-ready.bam
fi

SelectSampleEdits () {
	if [ ! $(find $workspace/"$sample".VarFilt.pass.dp.vcf  -type f -size +0c 2>/dev/null) ]; then
	echo "converting GVCF to VCF"
	gatk GenotypeGVCFs \
	-R $genome \
	--variant $workspace/"$sample"-gatk-haplotype.vcf.gz \
	-O $workspace/"$sample".vcf
	bgzip -f $workspace/"$sample".vcf
	tabix -p vcf $workspace/"$sample".vcf.gz

	echo "\nextracting edits from REDIportal database"
	bcftools view \
	-R $REDIportal_bed \
	$workspace/"$sample".vcf.gz | \
	bcftools sort -Oz -T $workspace \
	-o $workspace/"$sample"-rediportal.vcf.gz
	tabix -p vcf $workspace/"$sample"-rediportal.vcf.gz
	echo "\ndone extracting edits from REDIportal database"


    # GATK filter
        echo "Filtering variants from "$sample"-rediportal.vcf.gz"
        gatk VariantFiltration \
          -R $genome \
          -V $workspace/"$sample"-rediportal.vcf.gz \
          -O $workspace/"$sample".VarFilt.vcf \
          --window 35 \
          --cluster 3 \
          --filter-name "FS30" \
          --filter "FS > 30.0" \
          --filter-name "QD2" \
          --filter "QD < 2.0" \
          --genotype-filter-name "DP5" \
          --genotype-filter-expression "DP < 5"

        gatk SelectVariants \
             -V $workspace/"$sample".VarFilt.vcf \
             -O $workspace/"$sample".VarFilt.pass.vcf \
             --set-filtered-gt-to-nocall \
             --exclude-filtered

        grep "^#" $workspace/"$sample".VarFilt.pass.vcf > $workspace/"$sample".VarFilt.pass.dp.vcf
        grep -v "^#" $workspace/"$sample".VarFilt.pass.vcf | grep -v DP5 >> $workspace/"$sample".VarFilt.pass.dp.vcf

        echo "Done."

    fi
 }


AluEdits () {
    if [ ! $(find $workspace/"$sample"_alu_rnaedit_sites.vcf -type f -size +0c 2>/dev/null) ]; then
        # # alu RNA edit sites 
        echo "Performing SNPir filtering steps, except for RNA editing site removal"
        # convert vcf format into custom SNPiR format and filter variants with quality <20
        $SNPiR/convertVCF.sh $workspace/"$sample".VarFilt.pass.dp.vcf  $workspace/"$sample"_raw_rnaedits.txt 20
        # filt for total depth of at least 5 reads and 2 reads with alternate allele
        awk -F '\t|,' '{if (($3 >= 5) && ($4 >= 2)) {print}}' $workspace/"$sample"_raw_rnaedits.txt > $workspace/"$sample"_raw_rnaedits.mincov.txt

        echo "Filtering mismatches at read ends"
        # note: add the -illumina option if your reads are in Illumina 1.3+ quality format
        $SNPiR/filter_mismatch_first6bp.pl \
            -infile $workspace/"$sample"_raw_rnaedits.mincov.txt \
            -outfile $workspace/"$sample"_raw_rnaedits.rmhex.txt \
            -bamfile $workspace/"$sample"-ready.bam

        echo "Filtering for Alu edits" 
        awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/"$sample"_raw_rnaedits.rmhex.txt | \
            intersectBed \
            -wa \
            -header \
            -a stdin \
            -b $AluRegions \
            > $workspace/"$sample"_alu_rnaedit_sites.txt

        # # retrieve alu RNA edit sites from vars
        awk '{OFS="\t";print $1,$2,$2}' $workspace/"$sample"_alu_rnaedit_sites.txt | \
            intersectBed \
            -wa \
            -header \
            -a $workspace/"$sample".VarFilt.pass.dp.vcf  \
            -b stdin \
            | uniq \
            > $workspace/"$sample"_alu_rnaedit_sites.vcf

        echo "Done."
    fi
}


##############################################################################################################
NonAluEdits () {

    echo "Filtering for Nonalu edits" 
    if [ ! $(find $workspace/"$sample"_nonalu_rnaedit_sites.vcf -type f -size +0c 2>/dev/null) ]; then
        # Continue on to nonalu sites
        # filter out alu regions
        awk '{OFS="\t";$3=$2+1;print $0}' $workspace/"$sample"_raw_rnaedits.rmhex.txt | \
            intersectBed \
            -v \
            -a stdin \
            -b $AluRegions \
            > $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.txt

        # Remove edits in repeat regions
        awk '{OFS="\t";$2=$2-1"\t"$2;print $0}' $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.txt | \
            intersectBed \
            -v \
            -a stdin \
            -b $RepeatMasker | \
            cut -f1,3-7 > $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.txt

        # filter intronic sites that are within 4bp of splicing junctions
        # make sure your gene annotation file is in UCSC text format and sorted by chromosome and 
        # transcript start position
        $SNPiR/filter_intron_near_splicejuncts.pl \
            -infile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.txt \
            -outfile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.txt \
            -genefile $gene_annotation

        # filter variants in homopolymers
        $SNPiR/filter_homopolymer_nucleotides.pl \
            -infile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.txt \
            -outfile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
            -refgenome $genome

        # filter variants that were caused by mismapped reads
        # this may take a while depending on the number of variants to screen and the size of the reference genome
        # note: add the -illumina option if your reads are in Illumina 1.3+ quality format
        $SNPiR/BLAT_candidates.pl \
            -infile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.rmhom.txt \
            -outfile $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt \
            -bamfile $workspace/"$sample"-ready.bam \
            -refgenome $genome

        awk '{OFS="\t";print $1,$2-1,$2}' $workspace/"$sample"_raw_rnaedits.rmhex.nonalu.rmsk.rmintron.rmhom.rmblat.txt | \
        intersectBed \
            -wa \
            -header \
            -a $workspace/"$sample".VarFilt.pass.dp.vcf \
            -b stdin \
            | uniq \
            > $workspace/"$sample"_nonalu_rnaedit_sites.vcf

        echo "Done."
    fi
}

AnnotateFilter () {
    if [ ! $(find $workspace/"$sample"_alu_rnaedit_sites.vep.vcf -type f -size +0c 2>/dev/null) ]; then
        echo "Annotations alu edits"
        `which vep` \
            -i $workspace/"$sample"_alu_rnaedit_sites.vcf \
            -o $workspace/"$sample"_alu_rnaedit_sites.vep.vcf \
            --vcf \
            --verbose \
            --cache_version 105 \
            --cache \
	    --dir_cache $vep_cache \
            --everything \
            --fork 4 \
            --fasta $genome
    fi

    if [ ! $(find $workspace/"$sample"_alu_rnaedit_sites.vep.maf -type f -size +0c 2>/dev/null) ]; then
        perl `which vcf2maf.pl` \
            --tumor-id $sample \
            --normal-id $sample \
            --ref-fasta $genome \
            --ncbi-build GRCh38 \
            --inhibit-vep \
            --input-vcf $workspace/"$sample"_alu_rnaedit_sites.vep.vcf \
            --output-maf $workspace/"$sample"_alu_rnaedit_sites.vep.maf
    fi

    if [ ! $(find $workspace/"$sample"_nonalu_rnaedit_sites.vep.vcf -type f -size +0c 2>/dev/null) ]; then
        echo "Annotating nonalu edits"
        `which vep` \
            -i $workspace/"$sample"_nonalu_rnaedit_sites.vcf \
            -o $workspace/"$sample"_nonalu_rnaedit_sites.vep.vcf \
            --vcf \
            --verbose \
            --cache \
            --cache_version 105 \
            --dir_cache $vep_cache \
            --everything \
            --fork 4 \
            --fasta $genome
    fi

    if [ ! $(find $workspace/"$sample"_nonalu_rnaedit_sites.vep.maf -type f -size +0c 2>/dev/null) ]; then
        perl `which vcf2maf.pl` \
            --tumor-id $sample \
            --normal-id $sample \
            --ref-fasta $genome \
            --ncbi-build GRCh38 \
            --inhibit-vep \
            --input-vcf $workspace/"$sample"_nonalu_rnaedit_sites.vep.vcf \
            --output-maf $workspace/"$sample"_nonalu_rnaedit_sites.vep.maf
    fi

    if [ ! $(find $workspace/"$sample"_rnaedit_sites_final.tsv -type f -size +0c 2>/dev/null) ]; then
        echo "Filtering for final RNA edits"
        Rscript /expanse/lustre/projects/csd691/kfisch/RNA_editing_pipeline/scripts/filterRnaEdits.R \
            $workspace/"$sample"_alu_rnaedit_sites.vep.maf \
            $workspace/"$sample"_nonalu_rnaedit_sites.vep.maf \
            $workspace/"$sample"_rnaedit_sites_final.tsv
    fi
    echo "Done."
}

SelectSampleEdits > $workspace/logs/"$sample".log.1.filter.log 2>&1
echo "$sample variant filter exit code: $?" >> $workspace/logs/"$sample".exit.log

AluEdits  > $workspace/logs/"$sample".log.2.alu.log 2>&1
echo "$sample alu filtering exit code: $?" >> $workspace/logs/"$sample".exit.log

NonAluEdits  > $workspace/logs/"$sample".log.3.nonalu.log 2>&1
echo "$sample nonalu filtering exit code: $?" >> $workspace/logs/"$sample".exit.log

AnnotateFilter  > $workspace/logs/"$sample".log.4.annotate.log 2>&1
echo "$sample annotation and final filtering exit code: $?" >> $workspace/logs/"$sample".exit.log

grep -v "exit code: 0" $workspace/logs/"$sample".exit.log > $workspace/logs/"$sample".error.log

