# Correct Nanopore reads using assembled contigs.

# The long reads
reads=Q903_7

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
	Q903_7.minimap2.miniasm.gfa.png \
	Q903_7.minimap2.miniasm.minimap2.psitchensiscpmt_8.paf.gz

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
Q903_7.fq.gz: data/FAH26226-cleaned.fastq data/FAH26318-cleaned.fastq data/FAH26380-cleaned.fastq data/FAH26689-cleaned.fastq data/FAH26719-cleaned.fastq data/FAH26768-cleaned.fastq data/FAH26843-cleaned.fastq
	$(gzip) -c $^ >$@

# Rename scaffolds for minimap2, which requires the length is less than 255 characters.
Q903-ARCS_c4_l4_a0.5-8.rename.fa: %.rename.fa: %.fa
	sed 's/,/ /' $< >$@

# Convert FASTA to FASTQ.
%.fa.fq: %.fa
	bioawk -cfastx '{ $$qual = $$seq; gsub(".", "I", $$qual); print "@" $$name "\n" $$seq "\n+\n" $$qual }' $< >$@

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
	$(time) minimap2 -xmap-ont -r50000 $^ | $(gzip) >$@

# Align a FASTQ file to the indexed reference genome and produce a SAM file.
%.minimap2.$(ref).sam.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -xmap-ont -r50000 -a $^ | $(gzip) >$@

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
	$(time) miniasm -c1 -f $^ >$@

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

# Align the reads to the draft genome and produce a PAF file.
%.minimap2.$(reads).paf.gz: %.fa $(reads).fq.gz
	$(time) minimap2 -xmap-ont -w5 $^ | $(gzip) >$@

# Polish the assembly using Racon.
%.racon.fa: %.minimap2.$(reads).paf.gz $(reads).fq.gz %.fa
	gunzip -c $< | $(time) racon $(reads).fq.gz - $*.fa $@

# Bedtools

# Compute statistics on the depth of coverage of a BAM file.
%.bam.genomecov.tsv: %.bam $(ref).fa.fai
	(printf "Rname\tDepth\tCount\tRsize\tFraction\n"; \
		bedtools genomecov -g $(ref).fa.fai -ibam $<) >$@

# Compute statistics on the depth of coverage of a PAF file.
%.paf.genomecov.tsv: %.paf.gz $(ref).fa.fai
	(printf "Rname\tDepth\tCount\tRsize\tFraction\n"; \
		gunzip -c $< | cut -f 6,8,9 \
		| bedtools sort \
		| bedtools genomecov -g $(ref).fa.fai -i -) >$@

# Calculate depth of coverage statistics from bedtools genomecov.
%.genomecov.stats.tsv: %.genomecov.tsv
	mlr --tsvlite \
		then filter '$$Rname == "genome" && $$Depth > 0' \
		then step -a rsum -f Fraction \
		then put -q '@Depth_count += $$Count; if (is_null(@p25) && $$Fraction_rsum >= 0.25) { @p25 = $$Depth }; if (is_null(@p50) && $$Fraction_rsum >= 0.50) { @p50 = $$Depth }; if (is_null(@p75) && $$Fraction_rsum >= 0.75) { @p75 = $$Depth } end { emitf @Depth_count, @p25, @p50, @p75 }' \
		then rename p25,Depth_p25,p50,Depth_p50,p75,Depth_p75 \
		then put '$$Depth_IQR = $$Depth_p75 - $$Depth_p25' \
		$< >$@

# Compute the depth of coverage in bedgraph format.
%.paf.bedgraph: %.paf.gz $(ref).fa.fai
	gunzip -c $< | cut -f 6,8,9 \
	| bedtools sort \
	| bedtools genomecov -bga -g $(ref).fa.fai -i - >$@

# Select putative mitochondrial contigs

# Convert PAF to TSV.
%.paf.tsv: %.paf.gz
	gunzip -c $< \
	| mlr --nidx --fs tab --otsvlite \
		then rename 1,Qname,6,Tname,10,Matches,11,Length \
		then filter '$$Matches > 5000' \
		then stats1 -a sum -g Qname,Tname -f Matches,Length \
		then put '$$Identity = $$Matches_sum / $$Length_sum' \
		then sort -nr Identity >$@

# Select putative mitochondrial contigs.
%.paf.mt.id: %.paf.tsv
	mlr --itsvlite --onidx filter '$$Tname != "KU215903"' then cut -f Qname then uniq -f Qname then sort -f Qname $< >$@

# Generate a FASTA file of putative mitochondrial contigs.
%.minimap2.miniasm.minimap2.psitchensiscpmt_8.paf.mt.fa: %.minimap2.miniasm.minimap2.psitchensiscpmt_8.paf.mt.id %.minimap2.miniasm.fa
	samtools faidx $*.minimap2.miniasm.fa `<$<` | seqtk seq >$@

# GraphViz

# Render a graph to PDF using dot.
%.gv.dot.pdf: %.gv
	dot -Tpdf -o $@ $<
