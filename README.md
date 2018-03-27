# Experiments with Nanopore sequencing of Sitka spruce

- Map long reads to the organelles with Minimap2
- Map organelles to the long reads with Minimap2
- Map long reads to the nuclear genome with Minimap2
- Map nuclear genome to the long reads with Minimap2
- Classify the reads as plastid, mitochondrion, or nuclear
- Assemble the reads with Miniasm

# Data

- Nanopore sequencing `/projects/spruceup_scratch/psitchensis/Q903/data/reads/nanopore/*-cleaned.fastq`
- Nuclear genome `/projects/btl_scratch/hackathons/october2017/scaffolding/abyss-longscaff/abyss-longScaffold/Q903-ARCS_c4_l4_a0.5-8.fa`
- Organellar genomes v8 `/projects/btl/sjackman/picea-sitchensis-mitochondrion/psitchensiscpmt_8.fa`

| Run | Flowcell |
|----:|----------|
|   1 | FAH26843 |
|   2 | FAH26226 |
|   3 | FAH26689 |
|   4 | FAH26719 |
|   5 | FAH26768 |
|   6 | FAH26380 |
|   7 | FAH26318 |
|   8 | FAH44324 |
|   9 | FAH44332 |
|  10 | FAH44409 |
|  11 | FAH44462 |
|  12 | FAH44455 |
