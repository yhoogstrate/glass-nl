#!/usr/bin/env R


if(!exists("metadata.glass.per.resection")) {
  warning('metadata was not loaded')
  
  source('scripts/load_metadata.R')
}
