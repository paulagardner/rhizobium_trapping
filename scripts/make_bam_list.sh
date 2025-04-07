#!/bin/bash

cd /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/bamfiles || exit  #change shells to where the sorted bam files are
find "$PWD" -type f -name "*.sorted" > /mnt/c/Users/Paula/Desktop/TrappingCh2/rhizobium_trapping/data/2025sorted_bamfiles.txt