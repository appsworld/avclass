#!/usr/bin/env python3
'''
AVClass2 labeler
'''

import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, os.path.join(script_dir, 'lib/'))
sys.path.insert(1, os.path.join(script_dir, '../shared/'))
import argparse
from avclass2_common import AvLabels
from operator import itemgetter
import evaluate_clustering as ec
import json
import traceback
import gzip

# Default tagging file
default_tag_file = os.path.join(script_dir, "data/default.tagging")
# Default expansion file
default_exp_file = os.path.join(script_dir, "data/default.expansion")
# Default taxonomy file
default_tax_file = os.path.join(script_dir, "data/default.taxonomy")

def guess_hash(h):
    ''' Given a hash string, guess the hash type based on the string length '''
    hlen = len(h)
    if hlen == 32:
        return 'md5'
    elif hlen == 40:
        return 'sha1'
    elif hlen == 64:
        return 'sha256'
    else:
        return None

def format_tag_pairs(l, taxonomy=None):
    ''' Return ranked tags as string '''
    if not l:
        return ""
    if taxonomy is not None:
        p = taxonomy.get_path(l[0][0])
    else:
        p = l[0][0]
    out = "%s|%d" % (p, l[0][1])
    for (t,s) in l[1:]:
        if taxonomy is not None:
            p = taxonomy.get_path(t) 
        else:
            p = t
        out += ",%s|%d" % (p, s)
    return out

def list_str(l, sep=", ", prefix=""):
    ''' Return list as a string '''
    if not l:
        return ""
    out = prefix + l[0]
    for s in l[1:]:
        out = out + sep + s
    return out

def find_family(hash_ident="md5", vt=""):
    # Select hash used to identify sample, by default MD5
    hash_type = hash_ident

    # If ground truth provided, read it from file
    gt_dict = {}
    ifile_l = []
    av_labels = AvLabels(default_tag_file, default_exp_file, default_tax_file,None, False)
    if vt:
        ifile_l.append(vt)
        ifile_are_vt = True
    get_sample_info = av_labels.get_sample_info_vt_v3

    # Initialize state
    first_token_dict = {}
    token_count_map = {}
    pair_count_map = {}
    vt_all = 0
    avtags_dict = {}
    stats = {'samples': 0, 'noscans': 0, 'tagged': 0, 'maltagged': 0,
            'FAM': 0, 'CLASS': 0, 'BEH': 0, 'FILE': 0, 'UNK': 0}
    families_tags = []
    # Process each input file
    for ifile in ifile_l:
        fd = open(ifile, 'r')
        for line in fd:
            # If blank line, skip
            if line == '\n':
                continue

            # Read JSON line
            vt_rep = json.loads(line)

            # Extract sample info
            sample_info = get_sample_info(vt_rep)
            # If no sample info, log error and continue
            if sample_info is None:
                try:
                    name = vt_rep['md5']
                    sys.stderr.write('\nNo scans for %s\n' % name)
                except KeyError:
                    sys.stderr.write('\nCould not process: %s\n' % line)
                sys.stderr.flush()
                stats['noscans'] += 1
                continue

            # Sample's name is selected hash type (md5 by default)
            name = getattr(sample_info, hash_type)

            # If the VT report has no AV labels, output and continue
            if not sample_info.labels:
                sys.stdout.write('%s\t-\t[]\n' % (name))
                continue

            # Compute VT_Count (using list of AV engines if provided)
            vt_count = av_labels.get_sample_vt_count(sample_info)

            # Get the distinct tokens from all the av labels in the report
            # And print them. 
            try:
                av_tmp = av_labels.get_sample_tags(sample_info)
                tags = av_labels.rank_tags(av_tmp)

                # Collect stats
                # FIX: should iterate once over tags, 
                # for both stats and aliasdetect
                if tags:
                    stats["tagged"] += 1

                families_tags = format_tag_pairs(tags,av_labels.taxonomy)


            except:
                traceback.print_exc(file=sys.stderr)
                continue
    return families_tags


# if __name__=='__main__':
#     families_tags = find_family(vt="metadata.json")
