from simtk.openmm.app import PDBFile
import pmx
import itertools
import mdtraj
import pdbfixer.pdbfixer
import tempfile
import numpy as np


def build_pdb(sequence, filename, n_cap=None, c_cap=None, pH=7.0):
    """Build a PDB from a sequence and save to disk.
    
    Parameters
    ----------
    sequence : str
        String representation of protein sequence as 1 letter codes.
    filename : str
        name of output filename
    n_cap : str, optional, default=None
        Either None or "ACE"
    c_cap : str, optional, default=None
        Either None, "NME", or "NH2"
    pH : float, optional, default=7.0
        pH to use when building amino acids.
    """
    chain = pmx.Chain().create(sequence)
    
    if c_cap is not None:
        chain.add_cterm_cap()
    
    if n_cap is not None:
        chain.add_nterm_cap()

    temp_file = tempfile.NamedTemporaryFile(suffix=".pdb")
    temp_file.close

    chain.write(temp_file.name)
    
    # Now fix errors in element entries in CAP atoms
    # Also convert 
    traj = mdtraj.load(temp_file.name)
    top, bonds = traj.top.to_dataframe()

    if n_cap == "ACE":
        ind = np.where((top.name == "H3")&(top.resName == "ACE"))[0][0]
        top.element.ix[ind] = "H"

    if c_cap in ["NME", "NH2"]:
        ind = np.where((top.name == "H3")&(top.resName == "NME"))[0][0]
        top.element.ix[ind] = "H"

    if c_cap == "NH2":
        # Keep all atoms except the 3 NME methyl protons
        keep_ind = np.where((top.resName != "NME") | ((top.name != "H1") & (top.name != "H2") & (top.name != "H3")))[0]
        
        #Convert the NME carbon into a proton
        convert_ind = np.where((top.resName == "NME") & (top.name == "C"))[0][0]
        top.element.ix[convert_ind] = "H"
        top.name.ix[convert_ind] = "HN2"
        
        convert_ind = np.where((top.resName == "NME") & (top.name == "H"))[0][0]
        top.name.ix[convert_ind] = "HN1"
        
        
        top.resName.ix[np.where((top.resName == "NME"))[0]] = "NH2"
        
        traj._topology = mdtraj.Topology.from_dataframe(top, bonds)
        traj.restrict_atoms(keep_ind)

        top, bonds = traj.top.to_dataframe()

    if n_cap or c_cap:
        traj._topology = mdtraj.Topology.from_dataframe(top, bonds)
        
    traj.save(temp_file.name)  # Save output with fixed element names in caps.

    # Now fix missing charged termini.
    #structure = pdbfixer.pdbfixer.PdbStructure(open(temp_file.name))
    fixer = pdbfixer.pdbfixer.PDBFixer(temp_file.name)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(filename, 'w'))
