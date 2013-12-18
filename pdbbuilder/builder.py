from simtk.openmm.app import PDBFile
import pmx
import itertools
import mdtraj
import pdbfixer
import tempfile
import numpy as np

def build_pdb(sequence, filename, n_capped=False, c_capped=False, ph=7.0):
    chain = pmx.Chain().create(sequence)
    
    if c_capped:
        chain.add_cterm_cap()
    
    if n_capped:
        chain.add_nterm_cap()

    temp_file = tempfile.NamedTemporaryFile(suffix=".pdb")
    temp_file.close

    chain.write(temp_file.name)
    traj = mdtraj.load(temp_file.name)

    if n_capped or c_capped:
        top, bonds = traj.top.to_dataframe()

    if n_capped:
        ind = np.where((top.name == "H3")&(top.resName == "ACE"))[0][0]
        top.element.ix[ind] = "H"

    if c_capped:
        ind = np.where((top.name == "H3")&(top.resName == "NME"))[0][0]
        top.element.ix[ind] = "H"

    if n_capped or c_capped:
        traj._topology = mdtraj.Topology.from_dataframe(top, bonds)
        
    traj.save(temp_file.name)  # Save output with fixed element names in caps.

    # Now fix missing charged termini.
    structure = pdbfixer.PdbStructure(open(temp_file.name))
    fixer = pdbfixer.PDBFixer(structure)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(filename, 'w'))
