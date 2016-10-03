#!/bin/bash

gawk '{ l=length($2); L+=l; print l} END{ print "TOTAL: "L;}' all.cds.tbl
#800,322  (markov5)

gawk '{ l=length($2); L+=l; print l} END{ print "TOTAL: "L;}' all.intron.tbl

