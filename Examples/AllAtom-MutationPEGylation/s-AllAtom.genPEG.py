import KPrepInfo as KPI
import KPrepRead as KPR
import KPrepWrite as KPW
import KPrepFunctions as KPF

# USER VARIABLES
protName = "1xpb.91pdc"
pegRes = 91
numPEG = 113  # 113 = 5kDa

# PDB-BASED VARIABLES
# These will most likely not need to be changed. If you choose to use a different PDB as a monomer for
# polymer creation you will need to change them, or if you use an amino acid other than PAZ-DBCO for 
# polymer attachment then you will need to change them.
protRemoveAtom = ' H27'
protBondAtom = ' C08'
pegRemoveAtomFront = ' H1C'
pegRemoveAtomBack = ' H2C'
pegBondAtomFront = ' C1 '
pegBondAtomBack = ' C2 '


# Files to read
protFile = protName + ".pdb"  # Needs an accurate SEQRES section; HIS will be re-defined if necessary
pegFile = "PEG.pdb"
# Name of the files to write
outFile = protName + ".peg" + str(numPEG) + ".pdb"

###################################################
##### Do not change variables past this point #####
###################################################

#Read pdb and modify information
KPF.errorCheck_ATOM(protFile)
protAtom = KPR.read_pdb(protFile)
pegAtom = KPR.read_pdb(pegFile)

atomNames = {
    "protRemoveAtom": protRemoveAtom,
    "protBondAtom": protBondAtom,
    "pegRemoveAtomFront": pegRemoveAtomFront,
    "pegRemoveAtomBack": pegRemoveAtomBack,
    "pegBondAtomFront": pegBondAtomFront,
    "pegBondAtomBack": pegBondAtomBack
}

#Generate PEG
protAtom_PEG = KPF.pegGenAllAtom(protAtom, pegAtom, pegRes, atomNames, numPEG)

#Write new PDB
KPW.write_pdb(outFile, protAtom_PEG)
