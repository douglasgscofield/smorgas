smorgas
=======

`smorgas` eats the rich tasty treat that is `samtools mpileup` output and does a number of things with it.  The current target of `smorgas` is the assessment of variation in genome assembly projects, not for general SNP/indel/CNV calling.  Currently `smorgas` requires raw pileup with position-specific mapping quality (`-s` but neither `-g` nor `-u`).

**Nothing here is ready for production yet.  It may not even compile :-)**

Currently underway or planned are:

* Empirical mappability and mapping-quality tracts: mappability can be assessed from a k-mer standpoint, but empirical results from read mapping can tell us much more: variation in mapping quality, multiple mappings (using output from a planned [`yoruba`][yoruba] command), coverage variation, etc.; as well as where empirical and theoretical assessments do and do not match.

* Mapping-quality aware variant calls: if you are building a genome for a diploid organism from haplotypes drawn from the sequenced individual, say from [BACs][], [fosmid pools][] or a plant [gametophyte][], and you map reads from diploid tissue from that same individual, then half the pool of reads, those drawn from the haplotype used for assembly, should map close to perfectly, while the other half may map less well due to variants not present in the assembly.  Unsurprisingly this is an idealized explanatory model, in practice it is messier.

* Read-depth aware variant calls: in repeat- and/or CNV-heavy genomes, you expect read depth to be higher in assembled regions that collapse these variants.  Where are these regions and what does the combined assessment of read depth, read continuity, mapping quality and read contiguity (this last bit is likely to be the hardest) tell us about the number of variants and the variation within among them?

* Library assessment: If you have mapped reads from several libraries and have tagged reads from different libraries with different read group tags (hey, why not use [`yoruba readgroup`][yorubareadgroup] for that?) or simply have them in different BAM files, `smorgas` can split up the assessments above by library so you can see how well they are performing.

`smorgas` is written in C++ and uses a new `PileupParser` class to ingest and serve pileup.  Development of both is moving forward pretty quickly.  A major usability goals is low memory usage regardless of reference genome size, fragmentation or read mapping depth, [which should be goals common to every bioinformatics project][rikerdictionary].

[BACs]:            http://en.wikipedia.org/wiki/Bacterial_Artificial_Chromosome
[fosmid pools]:    http://en.wikipedia.org/wiki/Fosmid
[gametophyte]:     http://en.wikipedia.org/wiki/Gametophyte
[yoruba]:          https://github.com/douglasgscofield/yoruba
[yorubareadgroup]: https://github.com/douglasgscofield/yoruba#readgroup
[rikerdictionary]: https://github.com/douglasgscofield/riker/blob/master/README.md#createsequencedictionary

