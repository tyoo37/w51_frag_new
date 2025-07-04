#!/bin/bash
for reg in W51-E W51-IRS2
do
    for band in B3 B6
    do
        python run_initial_dendro.py --region $reg --band $band
    done
done