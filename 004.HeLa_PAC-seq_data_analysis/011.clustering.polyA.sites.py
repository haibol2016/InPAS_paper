# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 12:02:14 2022

@author: liuh

"""
import collections
import copy, os, re
def makehash():
    return collections.defaultdict(makehash)

path = r"C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa\PAC-seq.polyA.sites"
files = os.listdir(path)
os.chdir(path)
out_files = [re.sub(r'simplified', "simplified.clustered", f) for f in files]

for f in range(len(files)):
    file = files[f]
    out_file = out_files[f]
    polya_sites = makehash()
    clustered_polya = makehash()
    with open(file, "rt") as f:
        for line in f:
            entry = line.rstrip().split(sep = "\t")
            # nested dictionarys: chr: strand: location: count
            polya_sites[entry[0]][entry[3]][entry[1]] = int(entry[2])
    
    ## cluster polyA sites by dominant ones and within 12 flanking sequences
    for chr in polya_sites:
        for strand in ["+", "-"]:
            if polya_sites[chr][strand]:
                site_rank_count = {k: v for k, v in sorted(polya_sites[chr][strand].items(), key=lambda item: item[1],  reverse = True)}
                # print(site_rank_count['153534600'])
                rank_site_count = {k: v for k, v in sorted(polya_sites[chr][strand].items(), key=lambda item: int(item[0]))}
                for site_dominant in site_rank_count:
                    if site_dominant in rank_site_count:
                        clustered_polya[chr][strand][site_dominant] = site_rank_count[site_dominant]
                        del rank_site_count[site_dominant]
                        rank_site_count_copy = copy.deepcopy(rank_site_count)
                        for site_sorted in rank_site_count:
                            dist = abs(int(site_dominant) - int(site_sorted))
                            if dist <= 12:
                                clustered_polya[chr][strand][site_dominant] += rank_site_count[site_sorted]
                                del rank_site_count_copy[site_sorted]
                        rank_site_count = copy.deepcopy(rank_site_count_copy)
    
    ## output
    with open(out_file, "w+") as of:
        for chr in clustered_polya:
            print(chr)
            for strand in clustered_polya[chr]:
                for loc in clustered_polya[chr][strand]:
                    print("\t".join([chr, loc, str(clustered_polya[chr][strand][loc]), strand]), 
                          file = of)
                    

                        
            

        
        