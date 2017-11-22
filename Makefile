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

miniasm: \
	Q903_4.minimap2.miniasm.gfa.png \
	Q903_4.minimap2.miniasm.minimap2.psitchensiscpmt_8.paf.gz

ifndef ref
%.psitchensiscpmt_8.paf.gz:
	$(MAKE) ref=psitchensiscpmt_8 $@

%.psitchensiscpmt_8.sam.gz:
	$(MAKE) ref=psitchensiscpmt_8 $@

%.Q903_ARCS.paf.gz:
	$(MAKE) ref=Q903_ARCS $@

%.Q903_ARCS.sam.gz:
	$(MAKE) ref=Q903_ARCS $@
endif

# Compress the data.
%.fq.gz: data/%-cleaned.fastq
	$(gzip) -c $< >$@

# Concatenate and compress the data.
Q903_4.fq.gz: data/FAH26226-cleaned.fastq data/FAH26689-cleaned.fastq data/FAH26719-cleaned.fastq data/FAH26843-cleaned.fastq
	$(gzip) -c $^ >$@

# Rename scaffolds for minimap2, which requires the length is less than 255 characters.
Q903-ARCS_c4_l4_a0.5-8.rename.fa: %.rename.fa: %.fa
	sed 's/,/ /' $< >$@

# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align a FASTA file to the reference genome using BWA-MEM.
%.bwa.$(ref).sam.gz: %.fa $(ref).fa.bwt
	bwa mem $(ref).fa $< | $(gzip) >$@

# minimap2

# Index a FASTA file.
%.fa.mmi: %.fa
	$(time) minimap2 -xmap-ont -d $@ $<

# Index a FASTQ file.
%.fq.mmi: %.fq.gz
	$(time) minimap2 -xava-ont -d $@ $<

# Align a FASTA file to the indexed reference genome and produce a PAF file.
%.minimap2.$(ref).paf.gz: $(ref).fa.mmi %.fa
	$(time) minimap2 -xmap-ont $^ | $(gzip) >$@

# Align a FASTA file to the indexed reference genome and produce a SAM file.
%.minimap2.$(ref).sam.gz: $(ref).fa.mmi %.fa
	$(time) minimap2 -xmap-ont -a $^ | $(gzip) >$@

# Align a FASTQ file to the indexed reference genome and produce a PAF file.
%.minimap2.$(ref).paf.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -xmap-ont $^ | $(gzip) >$@

# Align a FASTQ file to the indexed reference genome and produce a SAM file.
%.minimap2.$(ref).sam.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -xmap-ont -a $^ | $(gzip) >$@

# Overlap reads with Minimap2 and produce a PAF file.
%.minimap2.paf.gz: %.fq.gz
	$(time) minimap2 -xava-ont $< $< | $(gzip) >$@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam.gz
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

# Racon

# The draft assembly to correct.
draft=FAH26843.minimap2.miniasm

# Align a FASTA file to the indexed draft genome and produce a SAM file.
$(draft).minimap2.%.sam.gz: $(draft).fa.mmi %.fa
	$(time) minimap2 -a -xmap-ont $^ | samtools view -h -F4 | $(gzip) >$@

# Align a FASTQ file to the indexed draft genome and produce a SAM file.
$(draft).minimap2.%.sam.gz: $(draft).fa.mmi %.fq.gz
	$(time) minimap2 -a -xmap-ont $^ | samtools view -h -F4 | $(gzip) >$@

# Add fake quality values to a SAM file.
%.q.sam.gz: %.sam.gz
	gunzip -c $< \
		| awk -vOFS='\t' '/^@/ { print; next } { $$11 = $$10; gsub(".", "I", $$11); print }' \
		| $(gzip) >$@

# Call the consensus sequence using Racon.
# Add fake quality values for Racon.
$(draft).%.racon.fa: $(draft).%.q.sam.gz $(draft).fa
	racon --sam NA $< $(draft).fa $@
