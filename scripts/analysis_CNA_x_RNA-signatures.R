#!/usr/bin/env R

# see if the CNAs can explain (some part of) the expression of the DGE signatures

if(!exists('cnv') | !exists('cnv.metadata')) {
  source('scripts/load_genomic_alterations.R')
}



