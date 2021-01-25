# The code in this repository was not originally meant to be invoked from within
# python.  Instead I expected users to invoke the scripts from the shell.
# Hence, only a few of the modules contain useful functions that can be
# accessed from within python.  They are below:
from .resid import ResID
from .closest_line_points import ClosestLinePoints
from .coords2angles import Coords2AnglesLengths, Coords2Angles
from .coords2dihedrals import Coords2DihedralsAnglesLengths, Coords2Dihedrals
from .coords2projected_dihedrals import Coords2ProjectedDihedralsLengths, Coords2ProjectedDihedrals
from .helixAngleOmega import CalcOmegaFromThetaPhi, CalcOmega

# I no longer remember why I import "main" from the executable scripts.
# Perhaps these next few lines are unnecessary, but they seem to do no harm:
from .coords2angles import main
from .coords2projected_dihedrals import main
from .coords2dihedrals import main
from .coords2distances import main
from .coords2helixAngleOmega import main
from .dlpdbfile import main
from .dlpisces import main
from .dna_interleave_residues import main
from .download_pdbs import main
from .dssp2pdb import main
from .has_dna_heavy_atoms import main
from .has_helices import main
from .has_protein_heavy_atoms import main
from .has_rna_heavy_atoms import main
from .has_secondary_str import main
from .has_sheets import main
from .has_turns import main
from .merge_lines_periodic import main
from .pdb2coords_ave import main
from .pdb2coords import main
from .pdb2helix import main
from .pdb2sequence import main
from .pdb2sheet import main
from .pdb2turn import main
from .select_chains_with_dna import main
from .select_interval import main
from .strip_secondary_str import main
from .truncate_chars import main
from .truncate_tokens import main


__all__ = ['closest_points',
           'coords2angles',
           'coords2projected_dihedrals',
           'coords2dihedrals',
           'coords2distances',
           'coords2helixAngleOmega',
           'dlpisces',
           'dna_interleave_residues',
           'download_pdbs',
           'dssp2pdb',
           'has_dna_heavy_atoms',
           'has_helices',
           'has_protein_heavy_atoms',
           'has_rna_heavy_atoms',
           'has_secondary_str',
           'has_sheets',
           'has_turns',
           'helixAngleOmega',
           'merge_lines_periodic',
           'pdb2coords_ave',
           'pdb2coords',
           'pdb2helix',
           'pdb2sequence',
           'pdb2sheet',
           'pdb2turn',
           'select_chains_with_dna',
           'select_interval',
           'strip_secondary_str',
           'truncate_chars.',
           'truncate_tokens']
