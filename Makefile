# Correct Nanopore reads using assembled contigs.

# Long reads
reads=Q903_18.porechop

# Linked reads
lr=HYN5VCCXX_4.trimadap

# Number of threads.
t=16

# Parallel compression with pigz.
gzip=pigz -p$t

# Parameters of Miniasm
miniasm_c=3

# Parameters of Unicycler
k=59

# Parameters of ARCS
c=2
e=50000
r=0.05

# Parameters of LINKS
a=0.7
l=10

# Paths to programs.
pilon_jar=/gsc/btl/linuxbrew/Cellar/pilon/1.22/pilon-1.22.jar

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
time=command time -v -o $@.time

.DELETE_ON_ERROR:
.SECONDARY:

all: Q903_18.porechop.minimap2.c3.miniasm.racon.pglaucamt.minimap2.Q903_18.porechop.paf.m5000.bold.unicycler.minimap2.Q903_18.porechop.paf.m5000.flye.racon.al5ki90.pglaucamt.gfa
all: Q903_18.porechop.minimap2.c3.miniasm.racon.pglaucamt.minimap2.Q903_18.porechop.paf.m5000.bold.unicycler.minimap2.Q903_18.porechop.paf.m5000.flye.racon.al5ki90.pglaucamt.compact.renumber.unicycler-polish.gfa
all: Q903_18.porechop.minimap2.c3.miniasm.racon.pglaucamt.minimap2.Q903_18.porechop.paf.m5000.bold.unicycler.minimap2.Q903_18.porechop.paf.m5000.flye.racon.al5ki90.pglaucamt.compact.renumber.unicycler-polish.HYN5VCCXX_4.trimadap.bx.sort.bam.bai

assemblies:
	miniasm \
	miniasm_racon \
	canu \
	unicycler \
	unicycler_canu \
	unicycler_flye

arcs: \
	miniasm_arcs \
	canu_contigs_arcs \
	canu_unitigs_arcs \
	flye_arcs \
	unicycler_arcs \
	unicycler_canu_arcs

miniasm: Q903_18.porechop.minimap2.c$(miniasm_c).miniasm.gfa

miniasm_racon: Q903_18.porechop.minimap2.c$(miniasm_c).miniasm.racon.racon.fa

miniasm_arcs: Q903_18.porechop.minimap2.c$(miniasm_c).miniasm.racon.racon.arcs.fa

canu: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.canu.stamp

canu_contigs_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.canu.contigs.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

canu_unitigs_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.canu.unitigs.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.k$k.unicycler.fa

unicycler_canu: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k$k.unicycler.fa

unicycler_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.unicycler.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler_canu_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k$k.unicycler.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler_canu_tigmint: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k$k.unicycler.tigmint.fa

unicycler_canu_tigmint_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k$k.unicycler.tigmint.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler_flye: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.flye.k$k.unicycler.fa

flye_arcs: Q903_18.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

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

# NCBI

# Download the Picea sitchensis plastid FASTA.
psitchensiscp.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=KU215903.2&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Download the Picea glauca mitochondrion FASTA.
pglaucamt.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/AM/LKAM01/LKAM01.1.fsa_nt.gz | gunzip -c | seqtk seq >$@

# Download the Picea glauca mitochondrion GBF (GenBank flat file).
pglaucamt.gbf:
	curl ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/AM/LKAM01/LKAM01.1.gbff.gz | gunzip -c >$@

# Download the Pinus taeda mitochondrion FASTA.
ptaedamt.fa:
	curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=MF991879.1&db=nucleotide&rettype=fasta' | seqtk seq >$@

# Data

# Separate the largest mitochondrial scaffolds from the organellar assembly.
psitchensismt_8.fa: psitchensiscpmt_8.fa
	seqtk seq -L150000 $< >$@

# Compress the data.
%.fq.gz: data/%-cleaned.fastq
	$(gzip) -c $< >$@

# Concatenate and compress the data.
Q903_18.fq.gz: \
		data/FAH26226-cleaned.fastq \
		data/FAH26318-cleaned.fastq \
		data/FAH26380-cleaned.fastq \
		data/FAH26689-cleaned.fastq \
		data/FAH26719-cleaned.fastq \
		data/FAH26768-cleaned.fastq \
		data/FAH26843-cleaned.fastq \
		data/FAH44324-cleaned.fastq \
		data/FAH44332-cleaned.fastq \
		data/FAH44409-cleaned.fastq \
		data/FAH44455-cleaned.fastq \
		data/FAH44462-cleaned.fastq \
		data/FAH44750-cleaned.fastq \
		data/FAJ02862-cleaned.fastq \
		data/FAJ04671-cleaned.fastq \
		data/FAJ04750-cleaned.fastq \
		data/FAK02037-cleaned.fastq \
		data/FAK02424-cleaned.fastq
	$(gzip) -c $^ >$@

# Rename scaffolds for minimap2, which requires the length is less than 255 characters.
Q903-ARCS_c4_l4_a0.5-8.rename.fa: Q903-ARCS_c4_l4_a0.5-8.fa
	sed 's/,/ /' $< >$@

# Convert FASTA to FASTQ.
%.fa.fq: %.fa
	bioawk -cfastx '{ $$qual = $$seq; gsub(".", "I", $$qual); print "@" $$name "\n" $$seq "\n+\n" $$qual }' $< >$@

# Trimadap

# Trim adapter sequences using trimadap.
%.trimadap.fq.gz: %.fq.gz
	trimadap-mt -p$t -t1 $< | sed 's/^X$$/N/' | $(gzip) >$@

# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align a FASTA file to the reference genome using BWA-MEM.
%.bwa.$(ref).sam.gz: %.fa $(ref).fa.bwt
	bwa mem $(ref).fa $< | $(gzip) >$@

# Align linked reads to the draft genome and do not sort.
%.$(lr).bx.sortn.bam: %.fa.bwt $(lr).bx.fq.gz
	bwa mem -t$t -pC $*.fa $(lr).bx.fq.gz | samtools view -@$t -h -F4 -o $@

# EMA

# Map linked reads to the draft genome using EMA.
# Filter out reads without barcodes.
%.$(lr).bx.ema.sortn.bam: $(lr).bx.fq.gz %.fa.bwt
	gunzip -c $< | paste - - - - - - - - | grep "BX:Z:" | tr '\t' '\n' \
	| $(time) ema align -t$t -r $*.fa -1 /dev/stdin | samtools view -@$t -h -F4 -o $@

# minimap2

# Index a FASTA file.
%.fa.mmi: %.fa
	$(time) minimap2 -t$t -xmap-ont -d $@ $<

# Index a FASTQ file.
%.fq.mmi: %.fq.gz
	$(time) minimap2 -t$t -xava-ont -d $@ $<

# Align a FASTA file to the indexed reference genome and produce a PAF file.
%.minimap2.$(ref).paf.gz: $(ref).fa.mmi %.fa
	$(time) minimap2 -t$t -xmap-ont $^ | $(gzip) >$@

# Align a FASTA file to the indexed reference genome and produce a SAM file.
%.minimap2.$(ref).sam.gz: $(ref).fa.mmi %.fa
	$(time) minimap2 -t$t -xmap-ont -a $^ | $(gzip) >$@

# Align a FASTQ file to the indexed reference genome and produce a PAF file.
%.minimap2.$(ref).paf.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -r50000 $^ | $(gzip) >$@

# Align a FASTQ file to the indexed reference genome and produce a SAM file.
%.minimap2.$(ref).sam.gz: $(ref).fa.mmi %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -r50000 -a $^ | $(gzip) >$@

# Overlap reads with Minimap2 and produce a PAF file.
%.minimap2.paf.gz: %.fq.gz
	$(time) minimap2 -t$t -xava-ont -I100G $< $< | $(gzip) >$@

# Compute the depth of coverage of each target sequence.
%.paf.depth.tsv: %.paf.gz
	gunzip -c $< | sed 's/Consensus_//' | \
		mlr --nidx --fs tab --otsvlite \
		then rename 1,Qname,6,Tname,7,Tlength,10,Matches,11,Length \
		then filter '$$Matches >= 5000' \
		then stats1 -a count,median,sum -g Tname -f Tlength,Length \
		then cut -f Tname,Tlength_median,Length_count,Length_sum \
		then rename Tlength_median,Tlength,Length_count,RC \
		then put '$$Depth = $$Length_sum / $$Tlength' \
		then sort -f Tname >$@

# Sambamba

# Filter a BAM file by alignment score and number of mismatches.
%.as100.nm5.bam: %.bam
	sambamba view -t $t -F '[AS] >= 100 and [NM] < 5' -fbam -o $@ $<

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a query-name-sorted BAM file by target.
%.sort.bam: %.sortn.bam
	samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@ $<

# Sort a query-name-sorted BAM file by BX tag.
%.sortbx.bam: %.sortn.bam
	samtools sort -@$t -tBX -T$$(mktemp -u -t $@.XXXXXX) -o $@ $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam.gz
	samtools sort -@$t -T$$(mktemp -u -t $@.XXXXXX) -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# Convert a BAM file to FASTQ.
%.sortn.bam.fq.gz: %.sortn.bam
	samtools fastq -@$t -TBX $< | $(gzip) -p$t >$@

# Miniasm

# Assemble long reads with Miniasm.
%.minimap2.c$(miniasm_c).miniasm.gfa: %.fq.gz %.minimap2.paf.gz
	$(time) miniasm -c$(miniasm_c) -f $^ >$@

# Assemble long reads with Miniasm, filtering out contained reads.
%.minimap2.Rc$(miniasm_c).miniasm.gfa: %.fq.gz %.minimap2.paf.gz
	$(time) miniasm -R -c$(miniasm_c) -f $^ >$@

# Convert GFA to FASTA.
%.miniasm.fa: %.miniasm.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Canu

# Assemble reads with Canu.
%.canu.stamp: %.fq.gz
	canu -d $*.canu -p canu genomeSize=6m -nanopore-raw $<
	touch $@

# Symlink the Canu unitigs FASTA.
%.canu.unitigs.fa: %.canu.stamp
	ln -sf $*.canu/canu.unitigs.fasta $@

# Symlink the Canu unitigs GFA.
%.canu.unitigs.gfa: %.canu.stamp
	ln -sf $*.canu/canu.unitigs.gfa $@

# Symlink the Canu contigs FASTA.
%.canu.contigs.fa: %.canu.stamp
	ln -sf $*.canu/canu.contigs.fasta $@

# Symlink the Canu contigs GFA.
%.canu.contigs.gfa: %.canu.stamp
	ln -sf $*.canu/canu.contigs.gfa $@

# Add the sequence to a Canu GFA file.
%.seq.gfa: %.fa %.gfa
	seqtk seq $< | gawk -vOFS='\t' 'ARGIND == 1 { id = substr($$1, 2); getline; x[id] = $$1; next } $$1 == "S" { $$3 = x[$$2] } 1' - $*.gfa >$@

# Flye

# Assemble reads with Flye.
%.flye.stamp: %.fq.gz
	flye --version
	flye -t$t -g6m --nano-raw=$< -o $*.flye
	touch $*.flye.stamp

# Copy the FASTA file of edges.
%.flye.fa: %.flye.stamp
	seqtk seq $*.flye/2-repeat/graph_final.fasta >$@

# Symlink the GFA file of edges.
%.flye.gfa: %.flye.stamp
	ln -sf $*.flye/2-repeat/graph_final.gfa $@

# Align the long reads to the Flye assembly using minimap2.
%.flye.minimap2.long.paf.gz: %.flye.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont $^ | $(gzip) >$@

# Align the long reads to the Flye assembly using minimap2.
%.flye.minimap2.long.sam.gz: %.flye.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -a $^ | $(gzip) >$@

# Align the long reads to the Flye assembly using unicycler_align.
%.flye.unicycler-align.long.sam.gz: %.flye.fa %.fq.gz
	unicycler_align --threads $t --ref $*.flye.fa --reads $*.fq.gz --sam $*.flye.unicycler-align.long.sam
	$(gzip) $*.flye.unicycler-align.long.sam

# Align the long reads to the Flye assembly graph using Bandage querypaths.
%.flye.bandage-querypaths.long.tsv: %.flye.gfa %.fa
	Bandage querypaths $^ $*.flye.bandage-querypaths.long

# Convert a Bandage querypaths file to GraphViz format.
%.flye.bandage-querypaths.long.gv: %.flye.bandage-querypaths.long.tsv
	cut -f2 $< | awk 'NF == 4 { \
			sub(",", ""); \
			++n["\"" $$2 "\" -> \"" $$3 "\""]; \
			++f["\"" $$2 "\" -> \"" $$3 "\""]; \
			if (!sub("+$$", "-", $$2)) sub("-$$", "+", $$2); \
			if (!sub("+$$", "-", $$3)) sub("-$$", "+", $$3); \
			++n["\"" $$3 "\" -> \"" $$2 "\""] \
			++r["\"" $$3 "\" -> \"" $$2 "\""] \
		} \
		END { print "digraph {"; \
			for (e in n) { print e, "[n=" n[e], "label=\"n=" f[e] "+" r[e] "=" n[e] "\"]" }; \
			print "}" \
		}' >$@

# Polish the Flye assembly with long reads using Racon.
%.flye.racon.fa: %.fq.gz %.flye.minimap2.long.sam.gz %.flye.fa
	$(time) racon -t $t $^ >$@

# Polish the Racon-polished Flye assembly with short reads using unicycler_polish.
# Do not correct local misassemblies.
# The --no_fix_local option is added by PR https://github.com/rrwick/Unicycler/pull/136
%.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.unicycler-polish.stamp: \
		%.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.fa \
		%.HYN5VCCXX_4.trimadap.bx.sort.mt.1.fq.gz \
		%.HYN5VCCXX_4.trimadap.bx.sort.mt.2.fq.gz
	rm -rf $*.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.unicycler-polish
	mkdir $*.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.unicycler-polish
	cd $*.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.unicycler-polish \
	&& $(time) unicycler_polish --threads $t --no_fix_local \
		-a ../$*.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_18.porechop.paf.mt.flye.racon.fa \
		-1 ../$*.HYN5VCCXX_4.trimadap.bx.sort.mt.1.fq.gz \
		-2 ../$*.HYN5VCCXX_4.trimadap.bx.sort.mt.2.fq.gz \
		--pilon $(pilon_jar)
	touch $@

# Create a GFA file from the FASTA file polished by Racon.
%.racon.gfa: %.racon.fa %.gfa
	seqtk seq $< | gawk -vOFS='\t' 'ARGIND == 1 { id = substr($$1, 2); rc[id] = $$3; getline; s[id] = $$1; next } $$1 == "S" { $$3 = s[$$2]; $$4 = "LN:i:" length($$3); $$5 = rc[$$2] } /^[SL]/' - $*.gfa >$@

# Create a GFA file from the FASTA file polished by Unicycler-polish.
%.unicycler-polish.gfa: %.unicycler-polish.fa %.gfa
	seqtk seq $< | gawk -vOFS='\t' 'ARGIND == 1 { id = substr($$1, 2); getline; x[id] = $$1; next } $$1 == "S" { $$3 = x[$$2] } 1' - $*.gfa >$@

# Extract reads with split alignments.
%.flye.minimap2.long.split.fa: %.flye.minimap2.long.paf.gz %.fq.gz
	seqtk subseq $*.fq.gz <(gunzip -c $< | awk '$$10 >= 1' | cut -f1 | uniq -d) | seqtk seq -A >$@

# Align reads with split alignments to the assembly using Bandage and Blast.
%.flye.minimap2.long.split.bandage-querypaths.tsv: %.flye.gfa %.flye.minimap2.long.split.fa
	Bandage querypaths $^ $*.flye.minimap2.long.split.bandage-querypaths

# Align reads with split alignments to the assembly using unicycler_align.
%.flye.minimap2.long.split.unicycler-align.sam.gz: %.flye.fa %.flye.minimap2.long.split.fa
	unicycler_align --threads $t --ref $*.flye.fa --reads $*.flye.minimap2.long.split.fa --sam $*.flye.minimap2.long.split.unicycler_align.sam
	$(gzip) $*.flye.minimap2.long.split.unicycler_align.sam

# Unicycler-polish

# Select the first read of the read pair.
%.1.fq.gz: %.fq.gz
	seqtk dropse $< | seqtk seq -1 | $(gzip) >$@

# Select the second read of the read pair.
%.2.fq.gz: %.fq.gz
	seqtk dropse $< | seqtk seq -2 | $(gzip) >$@

%.unicycler-polish.stamp: %.fa %.$(lr).bx.sortn.bam.1.fq.gz %.$(lr).bx.sortn.bam.2.fq.gz
	rm -rf $*.unicycler-polish
	mkdir $*.unicycler-polish
	cd $*.unicycler-polish \
	&& $(time) unicycler_polish --threads $t --no_fix_local --pilon $(pilon_jar) \
		-a ../$*.fa \
		-1 ../$*.$(lr).bx.sortn.bam.1.fq.gz \
		-2 ../$*.$(lr).bx.sortn.bam.2.fq.gz
	touch $@

# Copy the final FASTA file from unicycler_polish.
%.unicycler-polish.fa: %.unicycler-polish.stamp
	seqtk seq $*.unicycler-polish/???_final_polish.fasta | sed 's/ LN:i:[0-9]*//' >$@

# Porechop

# Trim adapter sequence using Porechop.
%.porechop.fq.gz: %.fq.gz
	porechop -t $t -i $< | $(gzip) >$@

# Unicycler

# Assemble long reads using Unicycler.
%.bold.unicycler.stamp: %.fq.gz
	unicycler -t$t --mode bold -o $*.bold.unicycler -l $*.fq.gz
	touch -r $*.bold.unicycler/assembly.gfa $@

# Copy the FASTA file.
%.unicycler.fa: %.unicycler.stamp
	seqtk seq $*.unicycler/assembly.fasta >$@

# Symlink the GFA file.
%.unicycler.gfa: %.unicycler.stamp
	ln -sf $*.unicycler/assembly.gfa $@

# The long reads for Unicycler
unicycler_long=$(reads).minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.$(reads).paf.mt

# Convert a BAM file to FASTQ files.
%.1.fq.gz %.2.fq.gz %.s.fq.gz: %.bam
	samtools sort -@$t -n -tBX $< | samtools fastq -@$t -TBX -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -

# Symlink the long reads for Unicycler.
$(reads).minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.long.fq.gz: $(unicycler_long).fq.gz
	ln -s $< $@

# Assemble short and long reads using Unicycler.
%.unicycler.fa: %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --keep 3 -o $*.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz
	seqtk seq $*.unicycler/assembly.fasta >$@

# Assemble short and long reads using Unicycler for a fixed value of k.
%.k$k.unicycler.fa: %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --kmers=$k -o $*.k$k.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz
	seqtk seq $*.unicycler/assembly.fasta >$@

# Assemble short and long reads with Canu unitigs of the long reads.
%.canu.unitigs.bold.unicycler.fa: $(unicycler_long).canu.unitigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold -o $*.canu.unitigs.bold.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.unitigs.bold.unicycler/assembly.fasta >$@

# Assemble short and long reads with Canu contigs of the long reads.
%.canu.contigs.bold.unicycler.fa: $(unicycler_long).canu.contigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold -o $*.canu.contigs.bold.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.contigs.bold.unicycler/assembly.fasta >$@

# Assemble short and long reads with Canu contigs of the long reads for a specified value of k.
%.canu.contigs.k$k.normal.unicycler.fa: $(unicycler_long).canu.contigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode normal --kmers=$k -o $*.canu.contigs.k$k.normal.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.contigs.k$k.normal.unicycler/assembly.fasta >$@

# Assemble short and long reads with Canu unitigs of the long reads for a specified value of k.
%.canu.unitigs.k$k.unicycler.fa: $(unicycler_long).canu.unitigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --kmers=$k -o $*.canu.unitigs.k$k.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.unitigs.k$k.unicycler/assembly.fasta >$@

# Assemble short and long reads with Canu contigs of the long reads for a specified value of k.
%.canu.contigs.k$k.unicycler.fa: $(unicycler_long).canu.contigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --kmers=$k -o $*.canu.contigs.k$k.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.contigs.k$k.unicycler/assembly.fasta >$@

# Assemble short and long reads with Flye contigs of the long reads for a specified value of k.
%.flye.k$k.unicycler.fa: $(unicycler_long).flye.gfa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --kmers=$k -o $*.flye.k$k.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.flye.k$k.unicycler/assembly.fasta >$@

# Align the long reads to the Unicycler assembly.
%.canu.contigs.k$k.unicycler.long.paf.gz: %.canu.contigs.k$k.unicycler.fa %.long.fq.gz
	$(time) minimap2 -t$t -xmap-ont $^ | $(gzip) >$@

# Extract reads with split alignments.
%.canu.contigs.k$k.unicycler.long.split.fq.gz: %.canu.contigs.k$k.unicycler.long.paf.gz %.long.fq.gz
	seqtk subseq $*.long.fq.gz <(gunzip -c $< | awk '$$10 >= 5000' | cut -f1 | uniq -d) | $(gzip) >$@

# Bandage

# Separate the largest component of a GFA file.
%.comp1.gfa: %.gfa
	Bandage reduce $< $@ --scope aroundnodes --nodes 1 --distance 100

# Filter the assembly graph by segment length.
%.l150k.gfa: %.gfa
	Bandage reduce $< $@ --scope aroundnodes --nodes $$(awk '$$1 == "S" && length($$3) >= 150000 { printf "," $$2 }' $<)

# Convert GFA to FASTA.
%.miniasm.l150k.fa: %.miniasm.l150k.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Filter the assembly graph by segment length, and retain adjacent vertices.
%.l150kd1.gfa: %.gfa
	Bandage reduce $< $@ --scope aroundnodes --distance 1 --nodes $$(awk '$$1 == "S" && length($$3) >= 150000 { printf "," $$2 }' $<)

# Convert GFA to FASTA.
%.miniasm.l150kd1.fa: %.miniasm.l150kd1.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Filter the assembly graph by depth of coverage.
%.c0.6-1.2.gfa: %.gfa
	Bandage reduce $< $@ --scope depthrange --mindepth 0.5 --maxdepth 1.2

# Convert GFA to FASTA.
%.c0.6-1.2.fa: %.c0.6-1.2.gfa
	awk '/^S/ { print ">" $$2 " " $$4 " " $$5 "\n" $$3 }' $< >$@

# Separate a GFA file of putative mitochondrial contigs.
%.mt.gfa: %.gfa psitchensismt_8.fa
	Bandage reduce $< $@ --scope aroundblast --query psitchensismt_8.fa

# Separate a subgraph of segments that align well to Picea glauca mitochondrion.
%.pglaucamt.gfa: %.gfa pglaucamt.fa
	Bandage reduce $< $@ --scope aroundblast --query pglaucamt.fa --alfilter 5000

# Separate a subgraph of segments that align well to Picea glauca mitochondrion.
%.al5ki90.pglaucamt.gfa: %.gfa pglaucamt.fa
	Bandage reduce $< $@ --scope aroundblast --query pglaucamt.fa --alfilter 5000 --ifilter 90

# Convert GFA to FASTA.
%.pglaucamt.fa: %.pglaucamt.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Convert GFA to FASTA.
%.compact.fa: %.compact.gfa
	awk '/^S/ { print ">" $$2 "\n" $$3 }' $< >$@

# Extract reads with split alignments.
%.paf.m5000.split.fa: %.paf.m5000.fq.gz %.paf.gz
	seqtk subseq $< <(gunzip -c $*.paf.gz | awk '$$10 >= 100' | cut -f1 | uniq -d) | seqtk seq -A >$@

# Align split long reads to the assembly graph using Bandage querypaths.
%.minimap2.$(reads).paf.m5000.split.qp.tsv: %.gfa %.minimap2.$(reads).paf.m5000.split.fa
	Bandage querypaths $^ $*.minimap2.$(reads).paf.m5000.split.qp

# Convert a Bandage querypaths file to GraphViz format.
%.qp.gv: %.qp.tsv
	cut -f2 $< | awk 'NF == 4 { \
			sub(",", ""); \
			++n["\"" $$2 "\" -> \"" $$3 "\""]; \
			++f["\"" $$2 "\" -> \"" $$3 "\""]; \
			if (!sub("+$$", "-", $$2)) sub("-$$", "+", $$2); \
			if (!sub("+$$", "-", $$3)) sub("-$$", "+", $$3); \
			++n["\"" $$3 "\" -> \"" $$2 "\""] \
			++r["\"" $$3 "\" -> \"" $$2 "\""] \
		} \
		END { print "digraph {"; \
			for (e in n) { print e, "[n=" n[e], "label=\"n=" f[e] "+" r[e] "=" n[e] "\"]" }; \
			print "}" \
		}' >$@

# Render a GFA file to PNG using Bandage.
%.gfa.png: %.gfa
	Bandage image $< $@

# Render a GFA file to SVG using Bandage.
%.gfa.svg: %.gfa
	Bandage image $< $@

# Compute graph metrics using Bandage.
%.gfa.info: %.gfa
	Bandage info $< >$@

# Convert plain text to TSV.
%.gfa.info.tsv: %.gfa.info
	sed 's/: */:/g;s/ /_/g' $< \
	| mlr --ixtab --ips : --otsvlite put '$$Inverse_score = $$Node_count * (1 + $$Dead_ends)**2; $$Score = 1/$$Inverse_score; $$Filename = "$(*F).gfa"' >$@

# Aggregate SPAdes graph metrics.
%/assembly_graph.comp1.gfa.info.tsv: \
		%/k025_assembly_graph.comp1.gfa.info.tsv \
		%/k045_assembly_graph.comp1.gfa.info.tsv \
		%/k059_assembly_graph.comp1.gfa.info.tsv \
		%/k073_assembly_graph.comp1.gfa.info.tsv \
		%/k085_assembly_graph.comp1.gfa.info.tsv \
		%/k093_assembly_graph.comp1.gfa.info.tsv \
		%/k101_assembly_graph.comp1.gfa.info.tsv \
		%/k109_assembly_graph.comp1.gfa.info.tsv \
		%/k115_assembly_graph.comp1.gfa.info.tsv \
		%/k121_assembly_graph.comp1.gfa.info.tsv
	mlr --tsvlite cat $^ >$@

# Racon

# Correct the reads using Racon.
%.racon.fa: %.fq.gz %.minimap2.paf.gz
	$(time) racon -t $t -f $^ $< | tr '_' ' ' >$@

# Align the reads to the Canu assembly and produce a SAM file.
%.canu.unitigs.minimap2.reads.sam.gz: %.canu.unitigs.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 -a $^ | $(gzip) >$@

# Align the reads to the Canu assembly and produce a PAF file.
%.canu.unitigs.minimap2.reads.paf.gz: %.canu.unitigs.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 $^ | $(gzip) >$@

# Polish the Canu assembly using Racon.
%.canu.unitigs.racon.fa: %.fq.gz %.canu.unitigs.minimap2.reads.sam.gz %.canu.unitigs.fa
	$(time) racon -t $t $^ | tr '_' ' ' >$@

# Align the reads to the draft genome and produce a PAF file.
%.minimap2.c$(miniasm_c).miniasm.minimap2.reads.paf.gz: %.minimap2.c$(miniasm_c).miniasm.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 $^ | $(gzip) >$@

# Polish the assembly using Racon.
%.minimap2.c$(miniasm_c).miniasm.racon.fa: %.fq.gz %.minimap2.c$(miniasm_c).miniasm.minimap2.reads.paf.gz %.minimap2.c$(miniasm_c).miniasm.fa
	$(time) racon -t $t $^ | tr '_' ' ' >$@

# Align the reads to the filtered miniasm assembly and produce a PAF file.
%.minimap2.c$(miniasm_c).miniasm.l150kd1.minimap2.reads.paf.gz: %.minimap2.c$(miniasm_c).miniasm.l150kd1.fa %.fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 $^ | $(gzip) >$@

# Polish the filtered miniasm assembly using Racon.
%.minimap2.c$(miniasm_c).miniasm.l150kd1.racon.fa: %.fq.gz %.minimap2.c$(miniasm_c).miniasm.l150kd1.minimap2.reads.paf.gz %.minimap2.c$(miniasm_c).miniasm.l150kd1.fa
	$(time) racon -t $t $^ | tr '_' ' ' >$@

# Align the reads to the draft genome and produce a PAF file.
%.minimap2.$(reads).paf.gz: %.fa $(reads).fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 $^ | $(gzip) >$@

# Align the reads to the draft genome and produce a SAM file.
%.minimap2.$(reads).sam.gz: %.fa $(reads).fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 -a $^ | $(gzip) >$@

# Polish the assembly using Racon.
%.racon.fa: $(reads).fq.gz %.minimap2.$(reads).sam.gz %.fa
	$(time) racon -t $t $^ >$@

# Tigmint
tigmint_n=10

# Correct misassemblies using Tigmint.
%.tigmint.fa: %.fa
	tigmint-make tigmint t=$t span=$(tigmint_n) window=1000 draft=$* reads=$(lr).bx

# Convert BED to BAM.
%.$(lr).bx.as0.65.nm5.molecule.size2000.bed.bam: %.$(lr).bx.as0.65.nm5.molecule.size2000.bed %.fa.fai
	bedtools bedtobam -i $< -g $*.fa.fai | samtools sort -@$t -Obam -o $@

# ARCS

# Create a graph of linked contigs using ARCS.
%.$(lr).c$c_e$e_r$r.arcs_original.gv %.$(lr).c$c_e$e_r$r.arcs.dist.gv %.$(lr).c$c_e$e_r$r.arcs.dist.tsv: %.$(lr).bx.sortn.bam %.fa
	bin/arcs -s98 -c$c -l0 -z500 -m4-20000 -d0 -e$e -r$r -v \
		-f $*.fa \
		-b $*.$(lr).c$c_e$e_r$r.arcs \
		-g $*.$(lr).c$c_e$e_r$r.arcs.dist.gv \
		--tsv=$*.$(lr).c$c_e$e_r$r.arcs.dist.tsv \
		--barcode-counts=$<.barcode-counts.tsv \
		$<

# Convert the ARCS graph to LINKS TSV format.
%.$(lr).c$c_e$e_r$r.arcs.links.tsv: %.$(lr).c$c_e$e_r$r.arcs_original.gv %.fa
	bin/arcs-makeTSVfile $< $@ $*.fa

# Scaffold the assembly using the ARCS graph and LINKS.
%.$(lr).c$c_e$e_r$r.arcs.a$a_l$l.links.scaffolds.fa %.$(lr).c$c_e$e_r$r.arcs.a$a_l$l.links.assembly_correspondence.tsv: %.$(lr).c$c_e$e_r$r.arcs.links.tsv %.fa
	cp $< $*.$(lr).c$c_e$e_r$r.arcs.a$a_l$l.links.tigpair_checkpoint.tsv
	LINKS -k20 -l$l -t2 -a$a -x1 -s /dev/null -f $*.fa -b $*.$(lr).c$c_e$e_r$r.arcs.a$a_l$l.links

# Rename the scaffolds.
%.links.fa: %.links.scaffolds.fa
	gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' $< >$@

# Convert an ARCS distance estimate graph to GFA.
%.$(lr).c$c_e$e_r$r.arcs.dist.n$n.gv.gfa: %.fa.fai %.$(lr).c$c_e$e_r$r.arcs.dist.n$n.gv
	abyss-todot -e --gfa $^ >$@

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

# Filter alignments and select reads

# Identify long reads with a good alignment score.
%.paf.m2000.id: %.paf.gz
	gunzip -c $< | awk '$$10 >= 2000 {print $$1}' | sort -u >$@

# Select long reads with a good alignment score.
%.paf.m2000.fq.gz: %.paf.m2000.id $(reads).fq.gz
	seqtk subseq $(reads).fq.gz $< | $(gzip) >$@

# Identify long reads with a good alignment score.
%.paf.m5000.id: %.paf.gz
	gunzip -c $< | awk '$$10 >= 5000 {print $$1}' | sort -u >$@

# Select long reads with a good alignment score.
%.paf.m5000.fq.gz: %.paf.m5000.id $(reads).fq.gz
	seqtk subseq $(reads).fq.gz $< | $(gzip) >$@

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

# Select putative mitochondrial long reads.
%.porechop.paf.mt.id: %.porechop.paf.tsv
	mlr --itsvlite --onidx filter '$$Tname != "KU215903"' then cut -f Qname then uniq -f Qname then sort -f Qname $< >$@

# Select putative mitochondrial contigs.
%.psitchensiscpmt_8.paf.mt.id: %.psitchensiscpmt_8.paf.tsv
	mlr --itsvlite --onidx filter '$$Matches_sum >= 50000 && $$Tname != "KU215903"' then cut -f Qname then uniq -f Qname then sort -f Qname $< >$@

# Select putative mitochondrial contigs.
%.minimap2.psitchensiscpmt_8.mt.fa: %.minimap2.psitchensiscpmt_8.paf.mt.id %.fa
	samtools faidx $*.fa `<$<` | seqtk seq >$@

# Select putative mitochondrial long reads.
%.paf.mt.fq.gz: %.paf.mt.id $(reads).fq.gz
	seqtk subseq $(reads).fq.gz $< | $(gzip) >$@

# Select putative mitochondrial short reads.
%.HYN5VCCXX_4.trimadap.bx.sort.mt.bam: %.HYN5VCCXX_4.trimadap.bx.sort.bam %.HYN5VCCXX_4.trimadap.bx.sort.bam.bai %.minimap2.psitchensiscpmt_8.paf.mt.id
	samtools view -@$t -o $@ $< `<$*.minimap2.psitchensiscpmt_8.paf.mt.id`

# Renumber the segments of a GFA file, sorted by length.
%.renumber.gfa: %.gfa
	bin/gfapy-renumber -o $@ $<

# Convert GFA to FASTA.
%.renumber.fa: %.renumber.gfa
	awk '/^S/ { print ">" $$2 " " $$4 " " $$5 "\n" $$3 }' $< >$@

# bcftools

# Call variants using samtools mpileup.
%.HYN5VCCXX_4.trimadap.bx.vcf.gz: %.fa %.HYN5VCCXX_4.trimadap.bx.sort.bam
	samtools mpileup -vf $^ >$@

# Call variants using samtools mpileup and the quality filtered BAM file.
%.HYN5VCCXX_4.trimadap.bx.as100.nm5.vcf.gz: %.fa %.HYN5VCCXX_4.trimadap.bx.sort.as100.nm5.bam
	samtools mpileup -vf $^ >$@

# Call genotypes.
%.gt.vcf.gz: %.vcf.gz
	bcftools call -mvf GQ,GP -Oz $< >$@

# Filter high-quality biallelic SNV.
%.qual.vcf.gz: %.vcf.gz
	bcftools filter -S . -e 'QUAL < 50 || FMT/GQ < 50' $< | bcftools view -a | bcftools view -m2 -M2 -v snps -Oz >$@

# Filter high-quality homozygous ALT SNV.
%.homo.vcf.gz: %.vcf.gz
	bcftools filter -i 'FMT/GT="1/1"' $< -Oz >$@

# Compute stats of a VCF file.
%.vcf.gz.stats: %.vcf.gz
	bcftools stats $< >$@

# GraphViz

n=3

# Label vertices and edges
%.dist.label.gv: %.dist.gv
	sed 's/Consensus_//g' $< | gvpr -c 'N{label = sprintf("%s\\n%u bp", name, l)} E{label = sprintf("n=%u", n)}' >$@

# Filter scaffolds by length using gvpr.
%.l20k.gv: %.gv
	gvpr -i 'N[l >= 20000]' -o $@ $<

# Filter scaffolds by length using gvpr.
%.l50k.gv: %.gv
	gvpr -i 'N[l >= 50000]' -o $@ $<

# Filter edges by their attribute n.
%.n$n.gv: %.gv
	gvpr 'E[n >= $n]' -o $@ $<

# Render a graph to PDF using dot.
%.gv.dot.pdf: %.gv
	dot -Tpdf -o $@ $<

# Miller

# Convert TSV to CSV.
%.csv: %.tsv
	mlr --itsvlite --ocsv cat $< >$@

# BLAST

# Align sequences to the nt database using blastn.
%.nt.blastn.txt: %.fa
	blastn -num_threads $t -db nt -query $< -out $@

# Align sequences to the nt database using blastn and report TSV.
%.nt.blastn.tsv: %.fa
	(printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tsskingdom\tssciname\tstitle\n"; \
		blastn -num_threads $t -db nt -outfmt '7 std staxid sskingdom ssciname stitle' -query $<) >$@

# Align sequences to the nt database using blastn and report SAM.
%.nt.blastn.sam: %.fa
	blastn -num_threads $t -db nt -outfmt 17 -query $< >$@

# Align sequences to the nt database using blastn and output an organism report.
%.nt.blastn.organism: %.fa
	blastn -num_threads $t -db nt -outfmt 18 -query $< -out $@

# ABySS

# Calculate assembly contiguity metrics using abyss-fac.
%.fac.tsv: %.fa
	/gsc/btl/linuxbrew/Cellar/abyss/2.1.0_1-k128/bin/abyss-fac -t100000 -G5500000 $< >$@

# Aggregate assembly contiguity metrics.
%.fac.tsv: \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.bold.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k41.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k43.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k45.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k47.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k49.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k51.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k53.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k55.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k57.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k59.normal.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k59.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k61.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k63.unicycler.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k51.unicycler.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k53.unicycler.tigmint.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k55.unicycler.tigmint.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k57.unicycler.tigmint.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k59.unicycler.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k61.unicycler.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k61.unicycler.tigmint.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k63.unicycler.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv \
		%.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k63.unicycler.tigmint.HYN5VCCXX_4.trimadap.c2_e50000_r0.05.arcs.a0.7_l10.links.fac.tsv
	mlr --tsvlite cat $^ >$@
