#!/usr/bin/env python

"""
 Typical usage (example):

    dlpdbfile.py < list_of_pdb_codes.txt

Where the input file contains a list of pdb-codes, one per line:
100d
101d
102d
103d
104d
107d
108d
109d
 Every line in this file contains a 4-letter pdb identifier.
Each 4-letter pdb code is compared with a list of pdb codes downloaded so far
(which is initially read from the file "pdbs_most_recent.txt", if present).
If it is new 4-letter pdb code, then the corresponding pdb file is downloaded,
and appended to the end of the "pdbs_most_recent.txt" file.
The corresponding DSSP file is also downloaded, if available.

"""

# author: Andrew Jewett
g_program_name = __file__.split('/')[-1]
g_date_str = '2021-1-24'
g_version_str = '0.6.0'


import sys, urllib, urllib.request, random, time, gzip


def FileExists(fname):
    try:
        file = open(fname, 'r')
    except IOError:
        exists = False
    else:
        exists = True
    return exists



def ExtractFileName(url):
    # The name of the file we create will be the same as the last few 
    # characters of the URL. (Whatever comes after the last '/' character.)
    out_fname = url
    slash_loc = out_fname.rfind('/')
    if (slash_loc >= 0):  # if '/' was found, then remove text preceeding '/'
        out_fname = out_fname[slash_loc+1:]
        if (out_fname == ''):
            sys.stderr.write('Error: filename difficult to determine from URL:\n' +\
                                 '       \"'+url+'\"\n')
            sys.exit(-1)
    return out_fname



# This function downloads a file from a URL, and (if requested by the user)
# prints out a warning message if there was a problem.

import os

def DownloadFileTo(url, file_name, verbose_mode = True):
    # Open the input file (at url) for reading
    # Simple way:
    #   in_file = urllib.urlopen(url)
    # But the web site may return garbage if the file doesn't exist.
    # Instead, we check for errors first:
    try: in_file = urllib.request.urlopen(url)
    except urllib.error.URLError as e:
        if (verbose_mode):
            sys.stderr.write(str(e)+'\n')
            sys.stderr.write('   omitting file \"'+file_name+'\"\n')
    else:
        if (verbose_mode):
            sys.stderr.write('downloading file \"'+file_name+'\"\n')
        # Open the output file for writing
        out_file = open(file_name, 'wb')
        out_file.write(in_file.read())
        out_file.close()



# This version infers the file name from the URL.
def DownloadFile(url, verbose_mode = True):
    DownloadFileTo(url, ExtractFileName(url), verbose_mode)




# Read in the list of PDB files already downloaded:

pdbs_old = set([])
try:
    pdbs_old_file = open('pdbs_old.txt', 'r')
    for pdb_code in pdbs_old_file:
        pdb_code = pdb_code.strip()#eliminate trailing and preceeding whitespace
        pdb_code = pdb_code.lower()#pdb codes are case-insensitive
        if (len(pdb_code) != 4):
            sys.stderr.write('Error in in \"pdbs_old.txt\":\n')
            sys.stderr.write('      Invalid PDB-code: \"'+pdb_code+'\"\n')
            exit(-1)
        pdbs_old.add(pdb_code)
    pdbs_old_file.close()
except IOError:
    pass




def main():
    # Here we keep track of the list of pdb files which
    #  i) are in the current list of pdb files requested in 
    #     sys.stdin
    # ii) are in the current list, but are not downloaded yet 
    #     (new pdb files)
    pdbs_current_file = open('pdbs_most_recent.txt', 'w')
    pdbs_current = set([]) #entire list of pdb codes requested
    pdbs_new  = set([]) #a list of new pdb codes requested that were not in the old list
    pdbs_old_file = open('pdbs_old.txt', 'a')


    for line in sys.stdin:
        line = line.strip()
        if len(line) == 0:
            continue
        pdb_code = line.lower()  #(pdb codes are case-insensitive)

        if (pdb_code in pdbs_new):
            sys.stderr.write(pdb_code+' appears redundantly. skipping\n')
        elif (pdb_code in pdbs_old):
            if (not (pdb_code in pdbs_current)):
                sys.stderr.write(pdb_code+' downloaded already. skipping.\n')
                pdbs_current.add(pdb_code)
                pdbs_current_file.write(pdb_code+'\n')
                pdbs_current_file.flush()  #<- necessary in case we get interrupted
            else:
                sys.stderr.write(pdb_code+' appears redundantly. skipping\n')

        else:
            #Download the corresponding pdb file
            file_name = pdb_code+'.pdb.gz' # <- these are compressed files
            url = 'http://www.rcsb.org/pdb/files/'+file_name
            # Note, if the URL above files, try this one instead:
            # url = 'https://files.rcsb.org/download/'+file_name
            DownloadFileTo(url, file_name)

            #Unzip the pdb file:
            with gzip.open(file_name, 'rb') as f:
                file_content = f.read()
                f.close()
            f = open(pdb_code+'.pdb', 'wb')
            f.write(file_content)
            if not FileExists(pdb_code+'.pdb'):
                sys.stderr.write('Error: A problem occured when trying to download PDB code \"'+pdb_code+'\"\n'
                                 '       Delete this entry from the file \"pdbs_old.txt\", and rerun '+g_program_name+'\n')
                sys.exit(-1)

            #Optional: Download the corresponding DSSP file
            url = 'ftp://ftp.cmbi.ru.nl/pub/molbio/data/dssp/'+pdb_code+'.dssp'
            DownloadFile(url)
            if not FileExists(pdb_code+'.dssp'):
                sys.stderr.write("    (The old DSSP PDB is server flaking out again.\n"
                                 "     Don't worry.  DSSP files are not needed.)\n")

            #Keep track of the the pdbs we have downloaded so far:

            pdbs_current.add(pdb_code)
            pdbs_current_file.write(pdb_code+'\n')
            pdbs_current_file.flush()  #<- necessary in case we get interrupted
            pdbs_old_file.write(pdb_code+'\n')
            pdbs_old_file.flush()      #<- necessary in case we get interrupted
            pdbs_new.add(pdb_code)


    pdbs_current_file.close()

    def ChainFromPDBfile(chainID, pdb_file):
        for line in pdb_file:
            line_type = line[0:6]
            if line_type in set(['ATOM  ', 'HETATM', 'ANISOU', 'SIGATM', 'SIGUIJ']):
                if line[21:22] == chainID:
                    yield line
            elif line_type == "HET   ":
                if line[12:13] == chainID:
                    yield line
            elif line_type == "HELIX ":
                initChainID = line[19:20]
                endChainID  = line[31:32]
                if (initChainID == chainID) or (endChainID == chainID):
                    yield line
            elif line_type == "SHEET ":
                initChainID = line[21:22]
                endChainID  = line[32:33]
                if (initChainID == chainID) or (endChainID == chainID):
                    yield line
            elif line_type == "TURN  ":
                initChainID = line[19:20]
                endChainID  = line[30:31]
                if (initChainID == chainID) or (endChainID == chainID):
                    yield line
            elif line_type == "SEQRES":
                if line[11:12] == chainID:
                    yield line
            elif line_type == "TER   ":
                ter_chainID = line[21:22]
                if (ter_chainID == chainID) or (ter_chainID == " "):
                    yield line
            else:
                yield line
        return


if __name__ == "__main__":
    main()
