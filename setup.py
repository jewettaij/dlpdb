from setuptools import setup

setup(
  
  name='dlpdb',

  packages=['dlpdb'],

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  description='collect statistics from the entire PDB library',

  long_description='A collection of tiny scripts to help automate the process of downloading and extracting coordinates, sequences, and secondary-structure information from (a representative subset of) structures from the entire PDB library.  These scripts were originally used to estimate distances distances and angles between specific atoms in proteins and DNA.  The distribution of these distances and angles can be used to estimate force-field parameters for coarse grained models for these biomolecules.',

  url='https://github.com/jewettaij/dlpdb',

  download_url='https://github.com/jewettaij/dlpdb/archive/v0.5.0.zip',

  version='0.5.0',

  license='MIT',

  classifiers=['Development Status :: 4 - Beta',
               'License :: OSI Approved :: MIT License',
               'Environment :: Console',
               'Topic :: Scientific/Engineering :: Chemistry',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows',
               'Intended Audience :: Science/Research'],

  scripts=['dlpdb/scripts/extract_angles.sh',
           'dlpdb/scripts/extract_dihedrals.sh',
           'dlpdb/scripts/extract_projected_dihedrals.sh',
           'dlpdb/scripts/extract_distances.sh',
           'dlpdb/scripts/extract_helix_angles.sh',
           'dlpdb/scripts/extract_helix_dihedrals.sh',
           'dlpdb/scripts/extract_helix_projected_dihedrals.sh',
           'dlpdb/scripts/extract_helix_distances.sh',
           'dlpdb/scripts/extract_helix_resAveDistances.sh',
           'dlpdb/scripts/extract_resAveDistances.sh',
           'dlpdb/scripts/extract_sheet_angles.sh',
           'dlpdb/scripts/extract_sheet_dihedrals.sh',
           'dlpdb/scripts/extract_sheet_projected_dihedrals.sh',
           'dlpdb/scripts/extract_sheet_distances.sh',
           'dlpdb/scripts/extract_sheet_resAveDistances.sh',
           'dlpdb/scripts/extract_turn_angles.sh',
           'dlpdb/scripts/extract_turn_dihedrals.sh',
           'dlpdb/scripts/extract_turn_projected_dihedrals.sh',
           'dlpdb/scripts/extract_turn_distances.sh',
           'dlpdb/scripts/extract_turn_resAveDistances.sh',
           'dlpdb/scripts/move_membrane_proteins.sh',
           'dlpdb/scripts/move_missing_dna_heavy_atoms.sh',
           'dlpdb/scripts/move_missing_protein_heavy_atoms.sh',
           'dlpdb/scripts/move_missing_secondary_str.sh',
           'dlpdb/scripts/move_nmr_structures.sh',
           'dlpdb/scripts/move_non-dna.sh',
           'dlpdb/scripts/replace_all_secondary_str.sh',
           'dlpdb/scripts/replace_missing_secondary_str.sh'],

  entry_points={
      'console_scripts': ['coords2angles.py=dlpdb.coords2angles:main',
                          'coords2dihedrals.py=dlpdb.coords2dihedrals:main',
                          'coords2projected_dihedrals.py=dlpdb.coords2projected_dihedrals:main',
                          'coords2distances.py=dlpdb.coords2distances:main',
                          'coords2helixAngleOmega.py=dlpdb.coords2helixAngleOmega:main',
                          'dlpisces.py=dlpdb.dlpisces:main',
                          'download_pdbs.py=dlpdb.download_pdbs:main',
                          'dssp2pdb.py=dlpdb.dssp2pdb:main',
                          'fix_dna_residue_order_pdb.py=dlpdb.fix_dna_residue_order_pdb:main',
                          'has_dna_heavy_atoms.py=dlpdb.has_dna_heavy_atoms:main',
                          'has_helices.py=dlpdb.has_helices:main',
                          'has_protein_heavy_atoms.py=dlpdb.has_protein_heavy_atoms:main',
                          'has_rna_heavy_atoms.py=dlpdb.has_rna_heavy_atoms:main',
                          'has_secondary_str.py=dlpdb.has_secondary_str:main',
                          'has_sheets.py=dlpdb.has_sheets:main',
                          'has_turns.py=dlpdb.has_turns:main',
                          'helixAngleOmega.py=dlpdb.helixAngleOmega:main',
                          'merge_lines_periodic.py=dlpdb.merge_lines_periodic:main',
                          'pdb2coords_ave.py=dlpdb.pdb2coords_ave:main',
                          'pdb2coords.py=dlpdb.pdb2coords:main',
                          'pdb2helix.py=dlpdb.pdb2helix:main',
                          'pdb2sequence.py=dlpdb.pdb2sequence:main',
                          'pdb2sheet.py=dlpdb.pdb2sheet:main',
                          'pdb2turn.py=dlpdb.pdb2turn:main',
                          'pdb_interleave_residues.py=dlpdb.pdb_interleave_residues:main',
                          'select_chains_with_dna.py=dlpdb.select_chains_with_dna:main',
                          'select_interval.py=dlpdb.select_interval:main',
                          'strip_secondary_str.py=dlpdb.strip_secondary_str:main',
                          'truncate_chars.py=dlpdb.truncate_chars:main',
                          'truncate_tokens.py=dlpdb.truncate_tokens:main']},
  zip_safe=True,
  include_package_data=True
)

