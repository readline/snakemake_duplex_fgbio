#!/usr/bin/env python
import os,sys
import gzip
import numpy as np

def allocated(resource, rule, lookup, default="__default__"):
    """Pulls resource information for a given rule. If a rule does not have any information 
    for a given resource type, then it will pull from the default. Information is pulled from
    definitions in the cluster.json (which is used a job submission). This ensures that any 
    resources used at runtime mirror the resources that were allocated.
    :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
    :param rule <str>: rule to lookup its information
    :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
    :param default <str>: default information to use if rule information cannot be found
    :return allocation <str>: 
        allocation information for a given resource type for a given rule
    """
    try: 
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]
    
    return allocation

def ignore(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is met (i.e. True) then it will not try to run that rule. This function is the 
    inverse to provided(). 
    """
    if condition:
        # If condition is True, returns an empty list to prevent rule from running
        samplelist = []
    return samplelist

def get_read_length(file_path, num_reads=100):
    sequence_lengths = []
    try:
        # Check if the file is in gzip format
        with gzip.open(file_path, 'rt') as file:
            for _ in range(num_reads):
                # Skip the Sequence identifier line
                next(file)
                # Read the Sequence line and get its length
                sequence_line = next(file).strip()
                sequence_lengths.append(len(sequence_line))
                # Skip the Quality score identifier and Quality scores lines
                next(file)
                next(file)
    except OSError:
        # If not gzip, assume regular fastq format
        with open(file_path, 'r') as file:
            for _ in range(num_reads):
                # Skip the Sequence identifier line
                next(file)
                # Read the Sequence line and get its length
                sequence_line = next(file).strip()
                sequence_lengths.append(len(sequence_line))
                # Skip the Quality score identifier and Quality scores lines
                next(file)
                next(file)

    return int(np.mean(sequence_lengths))