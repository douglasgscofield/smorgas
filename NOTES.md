PileupParser
============

To Do
-----

1. Implement a read stack that tracks read-specific characteristics

This will allow handle read quality encoded in read starts, and will allow us to manage pileup with multiple samples by tracking sample origin for each read.

Read-specific characteristics would be read mapping quality, read group, etc.  Internally mapping quality values attached to read starts will be tracked here rather than in each strata unless `-s` is specifically specified, in which case the `-s` values will be trackec included in the strata and the read-start values will continue to be tracked in the read stack.

Smorgas
=======

To Do
-----

1. Implement a megagametophyte-assuming SNP caller

2. Implement a megagametophyte-assuming coverage tracker

This will (hopefully) enable the detection of recombination points between megagametophytes