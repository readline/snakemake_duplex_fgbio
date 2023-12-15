#!/usr/bin/env python
import os,sys

def main():
    try:
        inputs = sys.argv[1:]
    except:
        sys.exit(sys.argv[0]+' [input 1] [input 2] ... [input N]')
    
    if len(inputs) == 0:
        sys.exit(sys.argv[0]+' [input 1] [input 2] ... [input N]')
    
    with open(inputs[0]) as infile:
        print(infile.read().replace(',','\t').rstrip())
    for i in inputs[1:]:
        with open(i) as infile:
            infile.readline()
            print(infile.read().replace(',','\t').rstrip())
if __name__ == '__main__':
    main()
    