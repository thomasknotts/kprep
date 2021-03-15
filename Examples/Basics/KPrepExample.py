import KPrepInfo as KPI
import KPrepRead as KPR
import KPrepWrite as KPW
import KPrepFunctions as KPF

# -------------- Files to read ----------------
#All-Atom
allAtomPDB = "1bdd.pdb"
#Corse grain Go model
corseGrainPDB = "1bdd.go.pdb"
corseGrainPARAM = "1bdd.go.param"
corseGrainTOP = "1bdd.go.top"

# -------------- Files to write ---------------
#All-Atom
allAtomCharmmPDB = "NEW-1bdd.charmm.pdb"
#Corse Grain Go model - LAMMPS compatible
corseGrainDATA = "NEW-1bdd.go.data" # LAMMPS-readable data file
corseGrainMOD = "NEW-1bdd.nc.mod"   # LAMMPS-readable mod file, important for native contact interactions


# ------------------ Making a CHARMM-compatible pdb --------------------------------------
# Generate array of Atom objects 
arrAtom = KPR.read_pdb(allAtomPDB) # The read function creates an array of "atom" object.
                                   #   These object definitions are defined in KPrepInfo
                                   #   With these objects we can do a lot of things

# CHARMM won't read pre-protinated histidines - named "HIS" - we need to re-name to "HSD"
#Iterate over all atoms, find the histidines and re-name
for atom in arrAtom :
    if atom.resType == "HIS" : atom.resType = "HSD"


# CHARMM also wont read if the first atom doesn't correspond with the first molecule
#KPrepWrite can 'reset' the numbering of the molecule count if you want
KPW.write_pdb(allAtomCharmmPDB, arrAtom, reset=True)
# ----------------------------------------------------------------------------------------

# ---------------- Coarse-grain prep and manipulation ------------------------------------
# Read the Go-model pdb, generate an array of Go particles
arrGo = KPR.readGO_pdb(corseGrainPDB)

# Corse-grain data from the web server do not include residue amino acid names, added here
arrGo = KPF.addAAgo(arrGo, allAtomPDB, renameH="HSD")

# Web server also doen't include MW info
arrGo = KPF.addMWgo(arrGo, corseGrainTOP)

# To generate a LAMMPS data file we need to know all the bond, angle, dihedral, and nonbonded
#  interactions. Defined in KPI are objects for each of these forces.
arrB, arrA, arrD, arrN, arrNb = KPR.readGO_param(corseGrainPARAM)

# We now have all the information to create a LAMMPS data file
KPW.writeLAMMPS_GO_data(corseGrainDATA,arrAtom=arrGo, arrBond=arrB, \
                        arrAngle=arrA, arrDihedral=arrD, arrNonbonded=arrN)
KPW.writeLAMMPS_GO_nc(corseGrainMOD, arrNb)
# ----------------------------------------------------------------------------------------

