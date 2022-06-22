#!/usr/bin/env R

# see if the CNAs can explain (some part of) the expression of the DGE signatures


# load data ----


if(!exists('cnv') | !exists('cnv.metadata')) {
  source('scripts/load_genomic_alterations.R')
}




