#!/usr/bin/env python

"""
 Typical usage (example):

    dlpisces.py < cullpdb_pc40_res2.5_R1.0_d090206_chains9151

 This program processes lists of "culled pdb" files generated by the 
 PISCES server at 
   http://dunbrack.fccc.edu/PISCES.php
   http://dunbrack.fccc.edu/Guoli/pisces_download.php
 This server uses a method described in:
 G. Wang and R. L. Dunbrack, Jr. 
 "PISCES: a protein sequence culling server."
  Bioinformatics, 19:1589-1591, 2003

 These lists contain non-redundant subsets of all the PDB files
 in the rcsb repository (obtained by crystallography?) whose
 resolution (in Angstroms) and R-values lie below some 
 desired threshold.  Different lists correspond to
 different levels of tolerable sequence homology, as well
 as different tolerated cutoffs for resolution and R-values.

 Here is an except from one of these files:
"cullpdb_pc40_res2.5_R1.0_d090206_chains9151"
 -----------------------
IDs         length Exptl.  resolution  R-factor FreeRvalue
7ODCA       424  XRAY        1.600    0.20    0.23
2AXPA       173  XRAY        2.500    0.26    0.29
2VKEA       207  XRAY        1.620    0.25    0.30
1IHJA        98  XRAY        1.800    0.22    0.24
2NUTC       196  XRAY        2.300    0.22    0.27
1DQTA       117  XRAY        2.000    0.23    0.25

 Every line in this file (other than the comment line at the beginning)
 begins with a 4-letter pdb identifier, immediately followed by a 1-letter
 chain identifier.  Thus it's possible for the same 4-letter PDB indentifier
 to appear twice in this list (with different chains).
 Each 4-letter pdb code is compared with a list of pdb codes downloaded so far
(which is initially read from the file "pdbs_most_recent.txt", if present).
If it is new 4-letter pdb code, then the corresponding pdb file is downloaded,
and appended to the end of the "pdbs_most_recent.txt" file.
The corresponding DSSP file is also downloaded, if available.

Finally, the relevant chain is extracted from each PDB file and placed
in a new PDB file like 7odc_chainA.pdb

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





def main():
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




    # Here we keep track of the list of pdb files which
    #  i) are in the current list of pdb files requested in 
    #     sys.stdin
    # ii) are in the current list, but are not downloaded yet 
    #     (new pdb files)
    pdbs_current_file = open('pdbs_most_recent.txt', 'w')
    pdbs_current = set([]) #entire list of pdb codes requested
    pdbs_new  = set([]) #a list of new pdb codes requested that were not in the old list
    pdbs_old_file = open('pdbs_old.txt', 'a')


    pisces_list = sys.stdin.readlines()
    for line in pisces_list:

        #Check to make sure the line does not contain comments.
        #Comments include the string "R-factor".

        if (line.find('R-factor') == -1):
            pdb_code = line[0:4]
            if (len(pdb_code) != 4):
                sys.stderr.write('Error:  Invalid PDB-code: \"'+pdb_code+'\"\n')
                exit(-1)

            pdb_code = pdb_code.lower()  #(pdb codes are case-insensitive)

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

                #Unzip the pdb file
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


    # Now loop through the list of chains in the pisces_list again
    # and exctract each chain from it's corresponding PDB file.

    pdb_files_file = open('pdb_files.txt','wb')

    for line in pisces_list: #<-reading the html file generated by www.rcsb.edu

        #Check to make sure the line does not contain comments.
        #Comments include the string "R-factor".
        if (line.find('R-factor') == -1):
            pdb_code = line[0:4]
            pdb_code = pdb_code.lower()  #(pdb codes are case-insensitive)
            assert(len(pdb_code) == 4)
            chainID  = line[4:5]
            try:
                pdb_file = open(pdb_code+'.pdb', 'r')
            except IOError:
                sys.stderr.write('Error: A problem occured when trying to download PDB code \"'+pdb_code+'\"\n'
                                 '       Delete this entry from the file \"pdbs_old.txt\", and rerun dlpisces.py\n')
                sys.exit(-1)
            # If not already present, then create a new PDB 
            # file containing only the requested chain.
            new_filename = pdb_code+'_chain'+chainID+'.pdb'
            try:
                pdb_file_chain = open(new_filename, 'r')
                pdb_file_chain.close()
            except IOError:
                sys.stderr.write('creating file \"'+new_filename+'\"\n')
                pdb_file_chain = open(new_filename, 'w')
                for line in ChainFromPDBfile(chainID, pdb_file):
                    pdb_file_chain.write(line)
                pdb_file_chain.close()
            pdb_file.close()
            pdb_files_file.write(new_filename+'\n')
    pdb_files_file.close()

    # Here, we help the user keep track of the pdb files which 
    # were in the original list, but are not in the current list:

    pdbs_not_needed_file = open('pdbs_not_needed.txt', 'wb')
    for pdb_code in pdbs_old:
        if (not (pdb_code in pdbs_current)):
            pdbs_not_needed_file.write(pdb_code+'\n')
    pdbs_not_needed_file.close()



if __name__ == "__main__":
    main()
