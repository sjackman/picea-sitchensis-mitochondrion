# Correct Nanopore reads using assembled contigs.

# Number of threads.
t=16

# Parallel compression with pigz.
gzip=pigz -p$t

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
time=command time -v -o $@.time

.DELETE_ON_ERROR:
.SECONDARY:

all: miniasm

miniasm: FAH26843.minimap2.miniasm.gfa.png

ifndef ref
psitchensiscpmt_8.%.paf.gz:
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
%.fq.mmi: %.fq.gz
	$(time) minimap2 -d $@ $<

# Align a FASTA file to a FASTQ file and produce a PAF file.
$(ref).%.paf.gz: $(ref).fq.mmi %.fa
	$(time) minimap2 $^ | $(gzip) >$@

# Align a FASTA file to a FASTQ file and produce a SAM file.
$(ref).%.sam.gz: $(ref).fq.mmi %.fa
	$(time) minimap2 -a $^ | $(gzip) >$@

# Align a FASTQ file to a FASTA file and produce a PAF file.
$(ref).%.paf.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 $^ | $(gzip) >$@

# Align a FASTQ file to a FASTA file and produce a SAM file.
$(ref).%.sam.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -a $^ | $(gzip) >$@

# Overlap reads with Minimap2 and produce a PAF file.
%.minimap2.paf.gz: %.fq.gz
	$(time) minimap2 -xava-ont $< $< | $(gzip) >$@

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

# Miniasm

# Assemble reads with Miniasm
%.minimap2.miniasm.gfa: %.fq.gz %.minimap2.paf.gz
	$(time) miniasm -f $^ >$@

# Convert GFA to FASTA.
%.fa: %.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Bandage

# Render a GFA file to PNG using Bandage.
%.gfa.png: %.gfa
	Bandage image $< $@

# Render a GFA file to SVG using Bandage.
%.gfa.svg: %.gfa
	Bandage image $< $@
