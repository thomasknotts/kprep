import KPrepInfo as KPI
import KPrepRead as KPR
import KPrepWrite as KPW
import KPrepFunctions as KPF

# USER VARIABLES
protName = "1xpb" # The protion that will be mutated
resName = "PDC"   # The mutant residue name found in the PDB
resSeq = 91       # The residue in the protein you will mutate
secondResAtomName = ' CB '

# Files to read
protFile = protName + ".pdb"  # Needs an accurate SEQRES section; HIS will be re-defined if necessary
resFile = resName + ".pdb"
# Name of the files to write
outFile = protName + "." + str(resSeq) + resName.lower() + ".pdb"

###################################################
##### Do not change variables past this point #####
###################################################

#Read pdb and modify information
KPF.errorCheck_SEQREStoATOM(protFile)
KPF.errorCheck_ATOM(protFile)
protAtom = KPR.read_pdb(protFile)
resAtom = KPR.read_pdb(resFile)

#Generate PEG
protAtom_mutated = KPF.mutateResidue(protAtom, resAtom, resSeq, resName, secondResAtomName)

#Write new PDB
KPW.write_pdb(outFile, protAtom_mutated)
