# Correct Nanopore reads using assembled contigs.

# Number of threads.
t=16

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
time=command time -v -o $@.time

.DELETE_ON_ERROR:
.SECONDARY:

all:

ifndef ref
psitchensiscpmt_8.%.paf:
	$(MAKE) ref=psitchensiscpmt_8 $@
endif

# Rename scaffolds for minimap2, which requires the length is less than 255 characters.
Q903-ARCS_c4_l4_a0.5-8.rename.fa: %.rename.fa: %.fa
	sed 's/,/ /' $< >$@

# minimap2

# Index a FASTA file.
%.fa.mmi: %.fa
	$(time) minimap2 -d $@ $<

# Index a FASTQ file.
%.fastq.mmi: %.fastq
	$(time) minimap2 -d $@ $<

# Align a FASTA file to a FASTQ file and produce a PAF file.
$(ref).%.paf: $(ref).fastq.mmi %.fa
	$(time) minimap2 $^ >$@

# Align a FASTA file to a FASTQ file and produce a SAM file.
$(ref).%.sam: $(ref).fastq.mmi %.fa
	$(time) minimap2 -a $^ >$@

# Align a FASTQ file to a FASTA file and produce a PAF file.
$(ref).%.paf: $(ref).fa.mmi %.fastq
	$(time) minimap2 $^ >$@

# Align a FASTQ file to a FASTA file and produce a SAM file.
$(ref).%.sam: $(ref).fa.mmi %.fastq
	$(time) minimap2 -a $^ >$@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam
	samtools sort -@$t -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index $<
