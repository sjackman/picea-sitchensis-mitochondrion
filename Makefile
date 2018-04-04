# Correct Nanopore reads using assembled contigs.

# Long reads
reads=Q903_12.porechop

# Linked reads
lr=HYN5VCCXX_4.trimadap

# Number of threads.
t=16

# Parallel compression with pigz.
gzip=pigz -p$t

# Parameters of Miniasm
miniasm_c=2

# Parameters of ARCS
c=2
e=50000
r=0.05

# Parameters of LINKS
a=0.7
l=10

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
time=command time -v -o $@.time

.DELETE_ON_ERROR:
.SECONDARY:

all: \
	miniasm \
	miniasm_racon \
	canu \
	unicycler \
	unicycler_canu

arcs: \
	miniasm_arcs \
	canu_contigs_arcs \
	canu_unitigs_arcs \
	unicycler_arcs \
	unicycler_canu_arcs

miniasm: Q903_12.porechop.minimap2.c$(miniasm_c).miniasm.gfa

miniasm_racon: Q903_12.porechop.minimap2.c$(miniasm_c).miniasm.racon.racon.fa

miniasm_arcs: Q903_12.porechop.minimap2.c$(miniasm_c).miniasm.racon.racon.arcs.fa

canu: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_12.porechop.paf.mt.canu.contigs.fa

canu_contigs_arcs: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_12.porechop.paf.mt.canu.contigs.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

canu_unitigs_arcs: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.Q903_12.porechop.paf.mt.canu.unitigs.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.k51.unicycler.fa

unicycler_canu: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k51.unicycler.fa

unicycler_arcs: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.unicycler.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

unicycler_canu_arcs: Q903_12.porechop.minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.canu.contigs.k51.unicycler.HYN5VCCXX_4.trimadap.c$c_e$e_r$r.arcs.a$a_l$l.links.fa

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

# Separate the largest mitochondrial scaffolds from the organellar assembly.
psitchensismt_8.fa: psitchensiscpmt_8.fa
	seqtk seq -L150000 $< >$@

# Compress the data.
%.fq.gz: data/%-cleaned.fastq
	$(gzip) -c $< >$@

# Concatenate and compress the data.
Q903_12.fq.gz: \
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
		data/FAH44462-cleaned.fastq
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
	samtools sort -@$t -o $@ $<

# Sort a SAM file and produce a sorted BAM file.
%.sort.bam: %.sam.gz
	samtools sort -@$t -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index -@$t $<

# Convert a BAM file to FASTQ.
%.bam.fq.gz: %.bam
	samtools fastq -@16 -TBX $< | $(gzip) -p16 >$@

# Miniasm

# Assemble reads with Miniasm
%.minimap2.c$(miniasm_c).miniasm.gfa: %.fq.gz %.minimap2.paf.gz
	$(time) miniasm -c$(miniasm_c) -f $^ >$@

# Convert GFA to FASTA.
%.miniasm.fa: %.miniasm.gfa
	awk '/^S/ { print ">" $$2 " " $$4 "\n" $$3 }' $< >$@

# Canu

# Assemble reads with Canu.
%.canu.stamp: %.fq.gz
	canu -d $*.canu -p canu genomeSize=6m -nanopore-raw $<
	touch $@

# Symlink the Canu unitigs.
%.canu.unitigs.fa: %.canu.stamp
	ln -sf $*.canu/canu.unitigs.fasta $*.canu.unitigs.fa
	ln -sf $*.canu/canu.unitigs.gfa $*.canu.unitigs.gfa

# Symlink the Canu contigs.
%.canu.contigs.fa: %.canu.stamp
	ln -sf $*.canu/canu.contigs.fasta $*.canu.contigs.fa
	ln -sf $*.canu/canu.contigs.gfa $*.canu.contigs.gfa

# Porechop

# Trim adapter sequence using Porechop.
%.porechop.fq.gz: %.fq.gz
	porechop -t $t -i $< | $(gzip) >$@

# Unicycler

# The long reads for Unicycler
unicycler_long=$(reads).minimap2.c2.miniasm.racon.racon.minimap2.psitchensiscpmt_8.mt.minimap2.$(reads).paf.mt

# Convert a BAM file to FASTQ files.
%.1.fq.gz %.2.fq.gz %.s.fq.gz: %.bam
	samtools sort -@16 -n -tBX $< | samtools fastq -@16 -TBX -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -

# Symlink the long reads for Unicycler.
$(reads).minimap2.c2.miniasm.racon.racon.HYN5VCCXX_4.trimadap.bx.sort.mt.long.fq.gz: $(unicycler_long).fq.gz
	ln -s $< $@

# Assemble short and long reads using Unicycler.
%.unicycler.fa: %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --keep 3 -o $*.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz
	seqtk seq $*.unicycler/assembly.fasta >$@

# Assemble short and long reads using Unicycler for a fixed value of k.
%.k51.unicycler.fa: %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --keep 3 --kmer_count=2 --min_kmer_frac=0.4 --max_kmer_frac=0.4 -o $*.k51.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz
	seqtk seq $*.unicycler/assembly.fasta >$@

# Assemble short and long reads with a Canu assembly of the long reads.
%.canu.contigs.k51.unicycler.fa: $(unicycler_long).canu.contigs.fa %.1.fq.gz %.2.fq.gz %.s.fq.gz %.long.fq.gz
	unicycler -t$t --mode bold --keep 3 --kmer_count=2 --min_kmer_frac=0.4 --max_kmer_frac=0.4 -o $*.canu.contigs.k51.unicycler -1 $*.1.fq.gz -2 $*.2.fq.gz -s $*.s.fq.gz -l $*.long.fq.gz --existing_long_read_assembly $<
	seqtk seq $*.canu.contigs.k51.unicycler/assembly.fasta >$@

# Bandage

# Separate the largest component of a GFA file.
%.comp1.gfa: %.gfa
	Bandage reduce $< $@ --scope aroundnodes --nodes 1 --distance 100

# Separate a GFA file of putative mitochondrial contigs.
%.mt.gfa: %.gfa psitchensismt_8.fa
	Bandage reduce $< $@ --scope aroundblast --query psitchensismt_8.fa

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

# Align the reads to the draft genome and produce a PAF file.
%.minimap2.$(reads).paf.gz: %.fa $(reads).fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 $^ | $(gzip) >$@

# Align the reads to the draft genome and produce a SAM file.
%.minimap2.$(reads).sam.gz: %.fa $(reads).fq.gz
	$(time) minimap2 -t$t -xmap-ont -w5 -a $^ | $(gzip) >$@

# Polish the assembly using Racon.
%.racon.fa: $(reads).fq.gz %.minimap2.$(reads).sam.gz %.fa
	$(time) racon -t $t $^ >$@

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
%.HYN5VCCXX_4.trimadap.bx.sort.mt.bam: %.HYN5VCCXX_4.trimadap.bx.sort.bam %.minimap2.psitchensiscpmt_8.paf.mt.id
	samtools view -@$t -o $@ $< `<$*.minimap2.psitchensiscpmt_8.paf.mt.id`

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
