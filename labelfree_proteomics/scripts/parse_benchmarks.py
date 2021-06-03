#!/usr/bin/env python

import glob
keep_first_header = True
with open('benchmarks.txt','w') as outfile:

    for filename in glob.iglob('./**/*.benchmark.txt', recursive=True):
 
        print(filename)
        
        with open(filename, 'r') as infile:
            
            for line in infile:
                if "h:m:s" in line and keep_first_header == False:
                    continue
                
                keep_first_header = False
                if line.strip(): 
                    outfile.write(filename+"\t"+line)
