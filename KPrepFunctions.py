#
# KPrepFunctions.py
# --------------------
# Authors: Addison Smith, Josh Wilkerson, Derek Bush
# Contributors: Anthony Gillepsie
# Written to be run on python 2 or 3

#########################################################
#--------------- Index of Function ----------------------
#--------- Prep Scripts for Pre-processing --------------
# addAAgo             Adds-in proper AA info into your Go Atom structure array
# addMWgo             Adds-in MW data into your Go Atoms structure attay from mmtsb-generated top file
# pegGenGo            Generates a Go-model PEG starting on a desired AA resID
# mutateResidue       Mutates a given amino acid into given protein Atom array at specified AA resID
# pegGenAllAtom       Generates an all atom PEG starting on the specfied AA resID
# catPRM              concatenates all inputted .prm files into one .prm file
# catRTF              concatenates all inputted .rtf files into one .rtf file
# charmmPDBprep       Takes an un-modified .pdb from the website and generates one ready for CHARMM input
# lammpsRepX_inputWorld   Creates world variables for temperature and boxIDs necessary for RepX analysis
# errorCheck_SEQREStoATOM Opens a pdb file and compares the sequences in the SEQRES and ATOM sections 
# errorCheck_ATOM         Opens a pdb file and makes sure the ATOM section atoms are numbered correctly
#--------- Prep Scripts for Post-processing -------------
# mbarNC_screenDCD_to_nc
# mbarEner_screen_to_box
# mbarPrep_box_to_mbar
# mbarPrep_box_and_nc_to_mbar
# mbarPrep_box_and_nc_to_mbar_domain
# shiftDCD
# catDCD
# catLog
# excelPrep
#--------------------------------------------------------
# Note: to toggle VIM folds: put cursor at your location -> command mode -> "za"
#       to create a fold: command mode -> ex. ":12,34fold" the numbers are line numbers
#--------------------------------------------------------
#########################################################


#########################################################
import KPrepInfo as KPI
import KPrepRead as KPR
import KPrepWrite as KPW
import numpy as np
import re
import random
import copy
import glob as gb
from datetime import datetime
import os
#########################################################


############### glob sort patterns #######################
# Function used for sorting the glob on the penultimate element#{{{
def penult_elem(elem) :
    return int(elem.split(".")[-2])

# Function used for sorting the glob on the second element
def second_elem(elem) :
    return int(elem.split(".")[1])#}}}
##########################################################
    

#--------- Prep Scripts for Pre-processing --------------
##################### addAAgo ###########################
def addAAgo(arrAtom, fileNamePdb, renameH='') :#{{{
    # Description:  Add amino acid information into the Atom class for a Go-model protein

    # Inputs:
    # FileNameVar   This must be the KPI.FileName class object variable where
    #                    FleNameVar.pdb has been defined previoustly
    # AtomArray     This is the array of "Go" Atoms (ex. Atoms that come from pdb_go)
    # renameH       (optional) Meaning "rename histidine" If the indicated label is accpeted
    #                    any hisidines in the pdb will be renamed to the desired name
    #                    Current acceptable inputs: "HIS" "HSD"

    # Actions:
    # Will read the pdb file and find the sequence information in the SEQRES section
    # Will then rename the resType Atoms from arrAtom (these Atoms don't have correct 
    #       amino acid labeling) to be their appropriate name.
    # Renaming histidine to indicated label will not be done unless specified
    #       and it said label is within the library of accepted label:
    #       acceptable labels: HIS, HSD

    # Outputs:
    # An array of "KPI.Atoms" classvbariables that conatain appropriatley renamed resType
    
    # Set-up and Pre-define variables
    naturalAminoAcids = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY', \
                         'HIS','HSD','ILE','LEU','LYS','MET','PHE','PRO', \
                         'SER','THR','TRP','TYR','VAL']
    acceptableHistidines = ['HIS', 'HSD']
    AAline = ''
    i = 0   # Current Atom indentifier
    renameCheck = False
    orig_arrAtom = arrAtom.copy()   # Make a copy in case restructure was incorrect
    numAtoms = len(arrAtom)

    # Check to see if we are to rename histidine
    if renameH == '' : renameCheck = False
    elif acceptableHistidines.count(renameH) == 1 : renameCheck = True
    else :
        print("\nWARNING: The given Histidine rename variable (renameH='%s') is not accepted." %renameH)
        print("Continueing the method but without Histidine replacement.\n")
        renameCheck = False

    # Begin reading the pdb
    PDBfile = open(fileNamePdb, 'r')
    for line in PDBfile :
        comment = line[0:6]

        if comment == "SEQRES" :
            AAline = line[19:].rstrip()
            print("Reading Sequnce %s" % AAline)
            for aa in AAline.split() :
                # Redefine the resType
                if i >= numAtoms :
                    print("ERROR the number of amino acids in %s SEQRES section" % (fileNamePdb))
                    print("do not match the number of amino acids in you Go .pdb file.")
                    quit()
                if naturalAminoAcids.count(aa) == 1 :
                    arrAtom[i].resType = aa
                    i += 1
                else: 
                    print("\nWARNING residue %s was skipped and not counted;" %aa)
                    print("order might be compromised.\n")
                # Redefine again if renameing histidine
                if acceptableHistidines.count(aa) == 1 and renameCheck == True :
                    arrAtom[i-1].resType = renameH  #i-1 for correct referenceing
    # End For
    PDBfile.close()
    
    # Check to see if number of Atoms equals number of changes
    if len(arrAtom) == i :
        return arrAtom
    else :
        print("\nWARNING: Residue IDs do not match. No modifications were made!\n")
        return orig_arrAtom#}}}
#########################################################


##################### addMWgo ###########################
def addMWgo (arrAtom, fileNameTop) :#{{{
    # Description:  Add molecular weight information into the Atom class for a Go-model protein
    
    # Inputs:
    # arrAtom       This is the array of KPI.Atoms 
    # fileNameTop   This must be the top file as generated by mmtsb server

    # Actions:
    # Will read the .top file and find the MASS information
    # Will then rename each Atom's mw from arrAtom to be the appropriate mass.

    # Outputs:
    # An array of "KPI.Atoms" class variables that conatain appropriatley defined mw
    
    # Set-up and Pre-define variables
    orig_arrAtom = arrAtom.copy()
    i = 0   # Current Atom identifier

    # Begin reading the top
    TOPfile = open(fileNameTop, 'r')
    for line in TOPfile :
        if re.search('^MASS', line) :
            # Define the mw
            mwVal = float(line.split()[3])
            arrAtom[i].mw = mwVal
            i += 1
    # End For
    TOPfile.close()

    # Check to see if number of Atoms equals number of changes
    if len(arrAtom) == i :
        return arrAtom
    else :
        print("\nWARNING: Residue IDs do not match. No modifications were made!\n")
        return orig_arrAtom#}}}
#########################################################


#################### pegGenGo ###########################
def pegGenGo (arrAtom, arrBond, arrNonbonded, startAtom, kDa=-1, numPEG=-1, DBCOlinker=False, maxIt=50, fittingFactor=0) :#{{{
    # Description: Go-model PEG generation w/ Purely-Repusive parameters. Will NOT modify function inputs.
    
    # Inputs:
    # arrAtom       This is the array of "Go" Atoms (ex. Atoms that come from pdb_go)
    # arrBond       This is the array of atom Bonds read-in from param_go
    # arrNonbonded  This is the array of Nonbnonded interaction read-in from param_go
    # startAtom     This is the atomNum of the site you want the PEG attached to
    # kDa           Length of PEG chain in kDa units (kDa OR numPEG must be defined)
    # numPEG        Number of PEG monomer units to be added (kDa OR numPEG must be defined)
    # DBCOlinker    If your PEG has a bulky DBCO linker define this variabe to be True.
    # maxIt         Defines the maximum iteration, if the atomNum is not on the surface, you will not generate PEG

    # Actions:
    # There are many model parameters defined suring set-up
    # The script will randomly propose a PEG coordinate comming off the startAtom
    # After each successful PEG coordinate proposal, arrAtom, arrBond and arrNonbonded are appended
    # The script will repeate this process unitll all PEGs have been generated

    # Outputs:
    # A modified FileName variable and arrays of modified KPI.Atom, KPI.Bond, KPI.Nonbnoded that contain
    #   all the purely-repuslive data needed for the Go-model
    
    # Set-up and Pre-define variables
    new_arrAtom      = arrAtom.copy()
    new_arrBond      = arrBond.copy()
    new_arrNonbonded = arrNonbonded.copy()
    
    totAtom = len(new_arrAtom)  # Total number of Atoms passed-in
    totPEG  = 0                 # Total number of verified PEG Atoms 
    last_atomNum = new_arrAtom[-1].atomNum # The last used atomNum in your Atom array
    last_atomID = 0             # The highest atomID in your Atom array (fully defined in next two lines)
    for k in range(totAtom) :
        if last_atomID < new_arrAtom[k].atomID : last_atomID = new_arrAtom[k].atomID
    last_resNum = new_arrAtom[-1].resNum + 1 # The last used resNum in your Atom array
    last_resID = 0              # The highest resID in your Atom array (fully defined in next two lines)
    for k in range(totAtom) :
        if last_resID < new_arrAtom[k].resID : last_resID = new_arrAtom[k].resID

    totBond = len(new_arrBond)  # Total number of Bonds passed-in
    last_bondNum = new_arrBond[-1].bondNum # The last used bondNum in your Bond arrayy
    last_bondID = 0             # The highest bondID in your Bond array (fully defined in next two lines)
    for k in range(totBond) :
        if last_bondID < new_arrBond[k].bondID : last_bondID = new_arrBond[k].bondID

    totNonbonded = len(new_arrNonbonded)
    rPEG   = 3.7/2          # ang units, PEG monomer radius. PEG bond-length as defined by GaussView: 3.699
    rDBCO  = 5.2            # ang units, DBCO linker molecule radius. Linker radius defined by GaussView: 5.20
    mwPEG  = 44             # gm/mol, molecular weight of a single PEG monomer
    mwDBCO = 303            # gm/mol, molecular weight of a single DBCO linker
    bondEpsilPEG = 378.000      # Bond class LJ epsilon parameter; probably stronger, but untested. Value sufficient
    bondEpsilDBCO = 378.000     # Bond class LJ epsilon parameter; probably stronger, but untested. Value sufficient
    nonbondedEpsilPEG  = -0.000132  # Nonbonded class LJ epsilon parameter; a Go-model assumption
    nonbondedSigmaPEG  = 4.3 / KPI.nbRminToEten  # Nonbonded class LJ sigma parameter in the CHARMM form Rmin/2; 
                                    # based on an all-atom research paper
                                    # Citation:
    nonbondedEpsilDBCO = -0.000132  # Nonbonded class LJ epsilon parameter; a Go-model assumption
    nonbondedSigmaDBCO = (rDBCO + 1) / KPI.nbRminToEten  # Nonbonded class LJ sigma parameter in the CHARMM form Rmin/2; 
                                    # an approximation based on modeling trends

    PEG1 = [new_arrAtom[startAtom-1].coordinates[0], new_arrAtom[startAtom-1].coordinates[1], \
            new_arrAtom[startAtom-1].coordinates[2], new_arrNonbonded[startAtom-1].rmin2]
    PEG2 = [0.0, 0.0, 0.0, 0.0]     # Variable to track the "next PEG" potential coordinates
    Nbr  = [0.0, 0.0, 0.0, 0.0]     # Variable used when checking if "next PEG" clashes/overlaps with existing "neighbor" Atoms

    # Error checks
    if kDa == -1 and numPEG == -1 :
        print("\nERROR: Please define kDa OR numPEG. Did not generate a PEG.\n")
        return arrAtom, arrBond, arrNonbonded
    if kDa != -1 and numPEG != -1 :
        print("WARNING: Both kDa AND numPEG were defined, discarding kDa information")
        kDa = -1
    if kDa != -1 and numPEG == -1 :
        numPEG = int(kDa*(113/5))  # Conversion: 113 PEG monomers per 5 kDa PEG chain
    
    # Start PEG generation
    print("Generating a PEG chain %i monomers long starting on residue %i..." % (numPEG, startAtom))
    for i in range(numPEG) :
        countIt = 1     # Reset iteration counter
        clash   = True  # Reset bool. If True, the proposed PEG location clashes with a neighbor;
                        #  assume a randomly generated PEG clashes until proven otherwise

        # Define radius of next PEG monomer
        if i == 0 and DBCOlinker == True : PEG2[3] = rDBCO
        else :                             PEG2[3] = rPEG
        
        # Find an acceptable randomized position for the PEG
        while clash == True and countIt <= maxIt :
            # Define spherical coordinate in a theta and phi random direction
            theta = random.uniform(0,1)*np.pi   # Range of theta = [0, pi]
            phi   = random.uniform(0,1)*np.pi*2 # Range of phi [0, 2pi]
            r     = PEG1[3] + PEG2[3]
            if i == 0 and fittingFactor != 0 : r += fittingFactor
            
            # Convert to cartesian and shift into correct location
            PEG2[0] = PEG1[0] + r*np.sin(theta)*np.cos(phi)
            PEG2[1] = PEG1[1] + r*np.sin(theta)*np.sin(phi)
            PEG2[2] = PEG1[2] + r*np.cos(theta)

            # Check to see any clashing in the system
            for j in range(totAtom + totPEG-1) :
                # Define neighbor Atom to check
                Nbr = [new_arrAtom[j].coordinates[0], new_arrAtom[j].coordinates[1], \
                       new_arrAtom[j].coordinates[2], 0.0]
                if j < totAtom :
                    Nbr[3] = new_arrNonbonded[j].rmin2
                elif DBCOlinker == True and j == totAtom :
                    Nbr[3] = nonbondedSigmaDBCO
                else :
                    Nbr[3] = nonbondedSigmaPEG
                # Distance between PEG2 and Nbr
                dist1 = np.sqrt((PEG2[0]-Nbr[0])**2 + (PEG2[1]-Nbr[1])**2 + (PEG2[2]-Nbr[2])**2) 
                # Smallest distance between PEG2 and Nbr before clashing occurs
                dist2 = PEG2[3] + Nbr[3]
                # Clash check
                if dist2 > dist1 :
                    clash = True
                    break
                else :
                    clash = False

            # Increment iteration
            countIt += 1
        # End While

        # If no clashes, create the new PEG atom
        if clash == False and countIt <= maxIt :
            # Create a new KPI.Atom for the PEG
            n = totAtom + i
            new_arrAtom.append(KPI.Atom())
            new_arrAtom[n].type = 'ATOM  '
            new_arrAtom[n].coordinates = [PEG2[0], PEG2[1], PEG2[2]]
            new_arrAtom[n].atomNum = last_atomNum+1 + i
            if   DBCOlinker == True and i == 0 : new_arrAtom[n].atomID = last_atomID+1
            elif DBCOlinker == True and i >= 1 : new_arrAtom[n].atomID = last_atomID+2
            elif DBCOlinker == False           : new_arrAtom[n].atomID = last_atomID+1
            new_arrAtom[n].atomName = ' CCO'
            new_arrAtom[n].resNum = last_resNum+1 + i 
            new_arrAtom[n].resSeq = last_resNum+1 + i # True because Go-model
            if   DBCOlinker == True and i == 0 : new_arrAtom[n].resID = last_resID+1
            elif DBCOlinker == True and i >= 1 : new_arrAtom[n].resID = last_resID+2
            elif DBCOlinker == False           : new_arrAtom[n].resID = last_resID+1
            new_arrAtom[n].resType = 'PEG'
            new_arrAtom[n].chainType = new_arrAtom[startAtom-1].chainType    # This may have to change!!!
            if i == 0 and DBCOlinker == True : new_arrAtom[n].mw = mwDBCO
            else :                             new_arrAtom[n].mw = mwPEG
            new_arrAtom[n].atomType = "PEG"

            # Create a new KPI.Bond for the PEG
            n = totBond + i
            new_arrBond.append(KPI.Bond())
            new_arrBond[n].bondNum = last_bondNum+1 + i
            if   DBCOlinker == True and i == 0 : 
                new_arrBond[n].bondID = last_bondID+1
                new_arrBond[n].sigma  = rDBCO + new_arrBond[startAtom-1].sigma
                new_arrBond[n].epsil  = bondEpsilDBCO     # If paramiterized, may need to be recalculated
            elif DBCOlinker == True and i == 1 : 
                new_arrBond[n].bondID = last_bondID+2
                new_arrBond[n].sigma  = rDBCO + rPEG
                new_arrBond[n].epsil  = bondEpsilDBCO     # If paramiterized, may need to be recalculated
            elif DBCOlinker == True and i >= 2 : 
                new_arrBond[n].bondID = last_bondID+3
                new_arrBond[n].sigma  = rPEG + rPEG
                new_arrBond[n].epsil  = bondEpsilPEG     # If paramiterized, may need to be recalculated
            elif DBCOlinker == False :                            
                new_arrBond[n].bondID = last_bondID+1
                new_arrBond[n].sigma  = rPEG + rPEG
                new_arrBond[n].epsil  = bondEpsilPEG      # If paramiterized, may need to be recalculated
            if i == 0 : new_arrBond[n].bond = [startAtom, last_atomNum+1+i]
            else :      new_arrBond[n].bond = [last_atomNum+i, last_atomNum+1+i]

            # Create a new KPINonbonded for the PEG
            n = totNonbonded + i
            if   DBCOlinker == True and i == 0 :
                new_arrNonbonded.append(KPI.Nonbonded())
                new_arrNonbonded[n].atomID = last_atomID+1
                new_arrNonbonded[n].rmin2 = nonbondedSigmaDBCO
                new_arrNonbonded[n].epsil = nonbondedEpsilDBCO
            elif DBCOlinker == True and i == 1 :
                new_arrNonbonded.append(KPI.Nonbonded())
                new_arrNonbonded[n].atomID = last_atomID+2
                new_arrNonbonded[n].rmin2 = nonbondedSigmaPEG
                new_arrNonbonded[n].epsil = nonbondedEpsilPEG
            elif DBCOlinker == False and i == 0 :
                new_arrNonbonded.append(KPI.Nonbonded())
                new_arrNonbonded[n].atomID = last_atomID+1
                new_arrNonbonded[n].rmin2 = nonbondedSigmaPEG
                new_arrNonbonded[n].epsil = nonbondedEpsilPEG
            
            # Print status
            if (i+1) % 50 == 0 or i+1 == numPEG : print("Generated %i PEGs" % (i+1))

            # Shift to next PEG
            PEG1 = PEG2.copy()
            totPEG += 1
        # End create new PEG If statment
        
        if countIt >= maxIt : 
            print("\nCRITICAL ERROR: Max iterations hit! No PEG was generated!\n")
            quit()
            return arrAtom, arrBond, arrNonbonded
    # End For
    
    print("Successfully generated a PEG %i monomers long." % (i+1))
    return new_arrAtom, new_arrBond, new_arrNonbonded#}}}
#########################################################


#--------- All-Atom uAA mutation and PEGylation ---------
################# Functions and Classes #################
# Functions#{{{
class Monomer:#{{{
# Class to represent a PEG monomer. Contains Atoms (of the monomer) and the
# Atoms removed both from this PEG monomer and from the previous monomer as
# this Monomer is bound to the preceding Monomer (or protein, in the case of
# the first Monomer). This information is used in backtracking when generating a
# PEG. Also keeps track of the number of times the following Monomer has been
# removed and a new one has been built from this one (also used in backtracking).
    def __init__(self, atoms=None, atomRemovedPrev=None, atomRemoved=None):
        self.atoms = atoms
        self.atomRemovedPrev = atomRemovedPrev
        self.atomRemoved = atomRemoved
        self.timesBacktracked = 0 #Number of times this Monomer has been backtracked to with these atoms}}}
    
def getAtom(atoms, atomName, resSeq=-1):#{{{
# Returns the atom with the matching name and resSeq from the Atoms specified.
# If resSeq==-1 or is not specified, resSeq is ignored when finding the atom.
    for atom in atoms:
        if (resSeq == -1 or atom.resSeq == resSeq) and atom.atomName == atomName:
            return atom#}}}

def removeAtom(atoms, atomName, resSeq=-1):#{{{
# Deletes the Atom with the matching name and resSeq from the Atoms specified.
# If resSeq==-1 or is not specified, resSeq is ignored when finding the atom.
    for i, atom in enumerate(atoms):    
        if (resSeq == -1 or atom.resSeq == resSeq) and atom.atomName == atomName:
            deleteIndex = i
            break
    del atoms[deleteIndex]#}}}

def getCoords(atoms, name, resSeq=-1):#{{{
# Returns the coordinates of the Atom with the matching name and resSeq from
# the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored when 
# finding the atom.
    for atom in atoms:
        if (resSeq == -1 or atom.resSeq == resSeq) and atom.atomName == name:
            return copy.deepcopy(atom.coordinates)#}}}

def getCCoords(atoms, resSeq=-1):#{{{
# Returns the coordinates of the Atom with name ' C  ' and matching resSeq 
# from the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored
# when finding the atom.
    return getCoords(atoms, ' C  ', resSeq)#}}}    

def getNCoords(atoms, resSeq=-1):#{{{
# Returns the coordinates of the Atom with name ' N  ' and matching resSeq 
# from the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored
# when finding the atom.
    return getCoords(atoms, ' N  ', resSeq)#}}}    

def getCACoords(atoms, resSeq=-1):#{{{
# Returns the coordinates of the Atom with name ' CA ' and matching resSeq 
# from the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored
# when finding the atom.
    return getCoords(atoms, ' CA ', resSeq)#}}}    

def getRadius(atom):#{{{
# Returns the radius based on the element of the specified Atom.
# If element is not specified, attempts to guess the element
# based on the atomName. This will not work perfectly for Atoms
# with names longer than 2 letters, but should work for Atoms
# typically found in organic compounds.
    r_C = .70
    r_H = .53
    r_N = .65
    r_O = .60
    r_S = .100
    r_other = 1.80
    if atom.element != '':
        element = atom.element
    else:
        element = ' ' + atom.atomName.lstrip()[0] 
    if element == ' C': return r_C
    elif element == ' H': return r_H
    elif element == ' N': return r_N
    elif element == ' O': return r_O
    elif element == ' S': return r_S
    else:
        return r_other#}}}

def getUnitVector(vec):#{{{
# Returns the vector of magnitude 1 in the same direction
# as the specified vector.
    mag = getMagVector(vec)
    u = vec/mag
    return u#}}}

def getMagVector(vec):#{{{
# Returns the magnitude of the specified vector.
    mag = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
    return mag#}}}

def getNtoCVector(atoms, resSeq=-1):#{{{
# Returns the vector between Atoms named ' N  ' and ' C  ' and matching 
# resSeq from the Atoms specified. If resSeq==-1 or is not specified, resSeq
# is ignored when finding the atom.
    N_coords = np.array(getNCoords(atoms, resSeq))
    C_coords = np.array(getCCoords(atoms, resSeq))
    vec = C_coords - N_coords
    return vec#}}}

def getNtoCUnitVector(atoms, resSeq=-1):#{{{
# Returns the vector of magnitude 1 and in the same direction as the 
# vector between Atoms named ' N  ' and ' C  ' and matching resSeq from
# the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored
# when finding the atom.
    vec = getNtoCVector(atoms, resSeq)
    return getUnitVector(vec)#}}}

def getCAtoXVector(atoms, name, resSeq=-1):#{{{
# Returns the vector between Atoms named ' CA ' and the specified name and 
# matching resSeq from the Atoms specified. If resSeq==-1 or is not specified,
# resSeq is ignored when finding the atom.
    CA_coords = np.array(getCACoords(atoms, resSeq))
    X_coords = np.array(getCoords(atoms, name, resSeq))
    vec = X_coords - CA_coords
    return vec#}}}

def getCAtoXUnitVector(atoms, name, resSeq=-1):#{{{
# Returns the vector of magnitude 1 and in the same direction as the 
# vector between Atoms named ' CA ' and the specified name and matching resSeq from
# the Atoms specified. If resSeq==-1 or is not specified, resSeq is ignored
# when finding the atom.
    vec = getCAtoXVector(atoms, name, resSeq)
    return getUnitVector(vec)#}}}

def getBondUnitVector(monAtom, pegs, numPEGadded, atomNames, protAtom=None, resSeq=-1):#{{{
# Returns the vector of magnitude 1 and in the same direction as the 
# bond between the monomer given (monAtom) and the previous monomer (or protein, in the case
# of the first monomer) from the Atoms specified.
    if numPEGadded == 0:
        oldBondAtom = getAtom(protAtom, atomNames["protBondAtom"], resSeq)
    else:
        oldBondAtom = getAtom(pegs[numPEGadded-1].atoms, atomNames["pegBondAtomBack"])
    newBondAtom = getAtom(monAtom, atomNames["pegBondAtomFront"])
    vect = newBondAtom.coordinates - oldBondAtom.coordinates
    return getUnitVector(vect)#}}}

def getAlignmentRotationMatrix(a, b):#{{{
# Returns rotation matrix to align unit vector a with unit vector b
# Matrix is designed to rotate a onto b, not vice versa
    I = np.identity(3)
    v = np.cross(a, b)
    v_x = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
    v_x_2 = np.dot(v_x, v_x)
    c = np.dot(a,b)
    R = I + v_x + v_x_2*(1/(1+c))
    return R#}}}

def shiftAtoms(atoms, shift):#{{{
# Shifts the specified atoms' coordinates by the shift vector specified
    for atom in atoms:
        atom.coordinates = atom.coordinates + shift#}}}

def scaleAtoms(atoms, scale):#{{{
# Multiplies the specified atoms' coordinates by the scale factor specified
    for atom in atoms:
        atom.coordinates = atom.coordinates * scale#}}}

def rotateAtoms(R, atoms):#{{{
# Rotates the specified atoms' coordinates using the rotation matrix specified
    for atom in atoms:
        atom.coordinates = np.dot(R, np.array(atom.coordinates))#}}}

def rotateResidue(resAtom, angle, u_vect):#{{{
# U_VECT MUST BE A UNIT VECTOR!!! 
# Rotates the specified atoms around unit vector u_vect by the specified angle.
# Treats atom with name ' N  ' as the origin when performing the rotation.
    ux = u_vect[0]
    uy = u_vect[1]
    uz = u_vect[2]
    cos_ang = np.cos(angle)
    sin_ang = np.sin(angle)
    rotationMatrix = np.array([[cos_ang+ux**2*(1-cos_ang), ux*uy*(1-cos_ang)-uz*sin_ang, ux*uz*(1-cos_ang)+uy*sin_ang],\
                        [uy*ux*(1-cos_ang)+uz*sin_ang, cos_ang+uy**2*(1-cos_ang), uy*uz*(1-cos_ang)-ux*sin_ang],\
                        [uz*ux*(1-cos_ang)-uy*sin_ang, uz*uy*(1-cos_ang)+ux*sin_ang, cos_ang+uz**2*(1-cos_ang)]])
    N_coords = getNCoords(resAtom)
    shiftAtoms(resAtom, -np.array(N_coords))
    for atom in resAtom:
        atom.coordinates = np.dot(atom.coordinates, rotationMatrix)
    shiftAtoms(resAtom, np.array(N_coords))#}}}

def removeResidue(protAtom, resSeq):#{{{
# Removes the residue with specified resSeq from the protein,
# excluding the atoms with the names listed below.
    atom_exclusions = {' N  ', ' C  ', ' HN ', ' O  '}
    removeAtomIndexes = []
    for i, atom in enumerate(protAtom):
        if atom.resSeq == resSeq and not atom.atomName in atom_exclusions:
            removeAtomIndexes.append(i)
    for index in reversed(removeAtomIndexes):
        protAtom.pop(index)#}}}

def alignResidue(resAtom, protAtom, resSeq):#{{{
# Aligns the atoms in resAtom with the residue they will be replacing
# in protAtom, specified byt resSeq. This includes rotating the atoms
# to align N and C in resAtom with the corresponding N and C in protAtom
# and shifting the atoms into the proper position. The residue is not
# inserted into the protein in this function, but its bacbone is positioned
# as if it were. There are likely clashes between the inserted residue
# and the rest of the protein though.
    a = getNtoCUnitVector(resAtom)  
    b = getNtoCUnitVector(protAtom, resSeq)
    R = getAlignmentRotationMatrix(a, b)
    mag_res = getMagVector(getNtoCVector(resAtom))
    mag_prot =  getMagVector(getNtoCVector(protAtom, resSeq))
    N_res_coords = getNCoords(resAtom)
    N_prot_coords = getNCoords(protAtom, resSeq)
    shiftAtoms(resAtom, -1*np.array(N_res_coords))
    scaleAtoms(resAtom, mag_prot/mag_res)
    rotateAtoms(R, resAtom)
    shiftAtoms(resAtom, np.array(N_prot_coords))#}}}

def fitResidue(allAtom, resSeq, secondAtomName=' CB '):#{{{
# Attempts to fit the residue specified by resSeq into the protein
# without clashes by calling rotateResidue and rotateSideChain. A
# printed warning statement is given if none of the combinations of
# rotating the entire residue and the sidechain avoid clashes. The 
# sidechain is rotated about the vector from the alpha carbon to
# the second atom in the residue, specified by secondAtomName.
    resRotations = 15
    sideChainRotations = 15
    resAngle = 2*np.pi/resRotations
    sideChainAngle = 2*np.pi/sideChainRotations
    resAtom = []
    protAtom = []
    for atom in allAtom:
        if atom.resSeq == resSeq:
            resAtom.append(atom)
        else:
            protAtom.append(atom)
    u_vect = getNtoCUnitVector(allAtom, resSeq)
    for i in range(resRotations):
        rotateResidue(resAtom, resAngle, u_vect)
        for j in range(sideChainRotations):
            u_vect2 = getCAtoXUnitVector(resAtom, secondAtomName)
            rotateSideChain(resAtom, sideChainAngle, u_vect2)
            if not resClashes(resAtom, protAtom):
                return
    print("WARNING: The inserted residue clashes with the rest of the protein. Output may be invalid.")#}}}

def insertResidue(resAtom, protAtom, resSeq):#{{{
# Inserts all Atoms in resAtom into protAtom except for those listed
# just below in atom_exclusions. No properties, including position, of
# the Atoms in resAtom are modified.
    atom_exclusions = {' N  ', ' C  ', ' HN ', ' O  '}
    insertion_atom_name = ' N  '
    for i, atom in enumerate(protAtom):
        if atom.resSeq == resSeq and atom.atomName == insertion_atom_name:
            insertIndex = i+1
            break
    for atom in reversed(resAtom):
        if not atom.atomName in atom_exclusions:
            protAtom.insert(insertIndex, atom)
    renumberAtoms(protAtom)#}}}

def rotateSideChain(resAtom, angle, u_vect):#{{{
# The non-backbone Atoms in resAtom are rotated by the specified angle
# about the vector specified (u_vect). This should be the vector from
# CA to the second atom in the sidechain, usually CB. 
    atom_exclusions = {' C  ', ' N  ', ' HN ', ' O  ', ' CA ', ' HA '}
    ux = u_vect[0]
    uy = u_vect[1]
    uz = u_vect[2]
    cos_ang = np.cos(angle)
    sin_ang = np.sin(angle)
    rotationMatrix = np.array([[cos_ang+ux**2*(1-cos_ang), ux*uy*(1-cos_ang)-uz*sin_ang, ux*uz*(1-cos_ang)+uy*sin_ang],\
                        [uy*ux*(1-cos_ang)+uz*sin_ang, cos_ang+uy**2*(1-cos_ang), uy*uz*(1-cos_ang)-ux*sin_ang],\
                        [uz*ux*(1-cos_ang)-uy*sin_ang, uz*uy*(1-cos_ang)+ux*sin_ang, cos_ang+uz**2*(1-cos_ang)]])
    CA_coords = getCACoords(resAtom)
    shiftAtoms(resAtom, -np.array(CA_coords))
    for atom in resAtom:
        if not atom.atomName in atom_exclusions:
            atom.coordinates = np.dot(atom.coordinates, rotationMatrix)
    shiftAtoms(resAtom, np.array(CA_coords))#}}}

def alignMonomer(protAtom, pegAtom, pegs, numPEGadded, atomNames, resSeq=-1):#{{{
# A new PEG Monomer is shifted and rotated so that its first carbon
# atom is set to be along the same vector as one of the hydrogens of
# the previous PEG Monomer. Additionally, one of the hydrogens of the
# new PEG Monomer is set along this same vector. The result is that
# the old PEG Monomer's C-H bond is now overlapped by an H-C bond of 
# the new PEG Monomer. The bond length between the last carbon of the
# previous PEG and the first carbon of the new PEG is set at 1.5. The
# PEG atoms are not added to the PEG array, byt they are aligned so
# they are ready to be.
    pegAtom = copy.deepcopy(pegAtom)
    newRemoveAtom = getAtom(pegAtom, atomNames["pegRemoveAtomFront"])
    newBondAtom = getAtom(pegAtom, atomNames["pegBondAtomFront"])
    if numPEGadded == 0:
        oldRemoveAtom = getAtom(protAtom, atomNames["protRemoveAtom"], resSeq)
        oldBondAtom = getAtom(protAtom, atomNames["protBondAtom"], resSeq)
    else:
        oldRemoveAtom = getAtom(pegs[numPEGadded-1].atoms, atomNames["pegRemoveAtomBack"])
        oldBondAtom = getAtom(pegs[numPEGadded-1].atoms, atomNames["pegBondAtomBack"])
    oldBondVector = np.array(oldBondAtom.coordinates) - np.array(oldRemoveAtom.coordinates) 
    oldBondLength = getMagVector(oldBondVector)
    newBondVector = np.array(newRemoveAtom.coordinates) - np.array(newBondAtom.coordinates)
    newBondLength = 1.5
    a = getUnitVector(newBondVector)
    b = getUnitVector(oldBondVector)
    R = getAlignmentRotationMatrix(a, b)
    shiftAtoms(pegAtom, -np.array(newRemoveAtom.coordinates))
    rotateAtoms(R, pegAtom)
    shiftAtoms(pegAtom, -np.array(newBondAtom.coordinates))
    newBondAtomCoords = (np.array(oldRemoveAtom.coordinates) - np.array(oldBondAtom.coordinates)) \
                         *(newBondLength/oldBondLength) + np.array(oldBondAtom.coordinates)
    shiftAtoms(pegAtom, newBondAtomCoords)
    return pegAtom#}}}

def fitMonomerAtom(monAtom, protAtom, pegs, numPEGadded, atomNames, resSeq=-1):#{{{
# Calls rotateMonomerAtom so that the atoms of monAtom are rotated about the vector
# between the previous PEG added and this one. This function should be called after
# alignMonomerAtom so that the rotation is done properly. The rotation is done 
# randomly and is attempted multiple times if there are clashes. Returns True if a 
# non-clashing rotation is found, and False if not. This leads to a randomized 
# direction of the PEG chain.
    monRotations = 15
    u_vect = getBondUnitVector(monAtom, pegs, numPEGadded, atomNames, protAtom, resSeq)
    for i in range(monRotations):
        angle = random.uniform(0, 2*np.pi)
        rotateMonomerAtom(monAtom, angle, u_vect, atomNames)
        if not pegClashes(protAtom, pegs, monAtom, numPEGadded, atomNames, resSeq):
            return True
    return False#}}}

def addMonomer(protAtom, pegAtom, pegs, numPEGadded, atomNames, resSeq=-1):#{{{
# Attemps to add another Monomer to the peg chain (pegs). It calls functions to 
# align the monomer and also to rotate it randomly until it does not clash. 
# Removes the atoms that are eliminated as the bond between PEG Monomers is
# formed from the main Atom array and stores them separately in the Monomer.
# Returns True if the Monomer is successfully added and False if not.
    monAtom = alignMonomer(protAtom, pegAtom, pegs, numPEGadded, atomNames, resSeq)
    if fitMonomerAtom(monAtom, protAtom, pegs, numPEGadded, atomNames, resSeq):
        atomRemoved = getAtom(monAtom, atomNames["pegRemoveAtomFront"])
        removeAtom(monAtom, atomNames["pegRemoveAtomFront"])
        if numPEGadded == 0:
            atomRemovedPrev = getAtom(protAtom, atomNames["protRemoveAtom"], resSeq)
            removeAtom(protAtom, atomNames["protRemoveAtom"], resSeq)
        else:
            atomRemovedPrev = getAtom(pegs[numPEGadded-1].atoms, atomNames["pegRemoveAtomBack"])
            removeAtom(pegs[numPEGadded-1].atoms, atomNames["pegRemoveAtomBack"])
        pegs[numPEGadded].atoms = monAtom
        pegs[numPEGadded].atomRemovedPrev = atomRemovedPrev
        pegs[numPEGadded].atomRemoved = atomRemoved
        return True
    else:
        return False#}}}

def rotateMonomerAtom(monAtom, angle, u_vect, atomNames):#{{{
# U_VECT MUST BE A UNIT VECTOR!!! 
# The atoms of monAtom are rotated about the specified unit vector (u_vect) by the
# specified angle. u_vect should be the unit bond vector between this Monomer and 
# the preceding one.
    ux = u_vect[0]
    uy = u_vect[1]
    uz = u_vect[2]
    cos_ang = np.cos(angle)
    sin_ang = np.sin(angle)
    rotationMatrix = np.array([[cos_ang+ux**2*(1-cos_ang), ux*uy*(1-cos_ang)-uz*sin_ang, ux*uz*(1-cos_ang)+uy*sin_ang],\
                        [uy*ux*(1-cos_ang)+uz*sin_ang, cos_ang+uy**2*(1-cos_ang), uy*uz*(1-cos_ang)-ux*sin_ang],\
                        [uz*ux*(1-cos_ang)-uy*sin_ang, uz*uy*(1-cos_ang)+ux*sin_ang, cos_ang+uz**2*(1-cos_ang)]])
    bondAtom = getAtom(monAtom, atomNames["pegBondAtomFront"])
    bondAtomCoords = copy.deepcopy(bondAtom.coordinates)
    shiftAtoms(monAtom, -bondAtomCoords)
    for atom in monAtom:
        atom.coordinates = np.dot(rotationMatrix, atom.coordinates)
    shiftAtoms(monAtom, bondAtomCoords)#}}}

def removeMonomer(pegs, numPEGadded, protAtom=None, resSeq=-1):#{{{
# Removes the atoms of the last Monomer added and updates the number
# of times the penultimate Monomer has been backtracked to.
    if numPEGadded == 1:
        protAtom.append(pegs[numPEGadded-1].atomRemovedPrev)
    else:
        pegs[numPEGadded-2].atoms.append(pegs[numPEGadded-1].atomRemovedPrev)
    pegs[numPEGadded-1].atoms = None
    pegs[numPEGadded-1].atomRemovedPrev = None
    pegs[numPEGadded-1].atomRemoved = None
    if numPEGadded > 1:
        pegs[numPEGadded-2].timesBacktracked+=1#}}}

def resClashes(resAtom, protAtom):#{{{
# Determines if the residue in its current position clashes with
# any other atoms in the protein, based on the atomic radius of 
# each atom.
    if selfClashes(resAtom):
        return True
    atom1_exclusions = {' C  ', ' N  ', ' CA '}
    atom2_exclusions = {''}
    for atom1 in resAtom:
        if not atom1.atomName in atom1_exclusions: 
            x1 = atom1.coordinates[0] 
            y1 = atom1.coordinates[1] 
            z1 = atom1.coordinates[2]
            r1 = getRadius(atom1)
            for atom2 in protAtom:
                if not atom2.atomName in atom2_exclusions:
                    x2 = atom2.coordinates[0] 
                    y2 = atom2.coordinates[1] 
                    z2 = atom2.coordinates[2]
                    r2 = getRadius(atom2)
                    distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                    minDistance = r1 + r2
                    if distance < minDistance:
                        print(atom1.atomName, atom1.resSeq, atom1.coordinates, atom2.atomName, atom2.resSeq, atom2.coordinates)
                        return True
    return False#}}}

def selfClashes(resAtom):#{{{
# Determines if the sidechain of the residue in its current 
# position clashes with the backbone of the same residue, based 
# on the atomic radius of each atom. This function is likely not
# extremely important anymore due to changes in the implementation
# of the uAA mutation functions, but it does no harm to the 
# mutation process.
    sidechain_names = {' C  ', ' O  ', ' N  ', ' HN ', ' HA '} 
    backbone_exclusions = {' CA '}
    sidechainAtom = []
    backboneAtom = []
    for atom in resAtom:
        if atom.atomName in sidechain_names:
            sidechainAtom.append(atom)
        elif not atom.atomName in backbone_exclusions:
            backboneAtom.append(atom)
    for atom1 in backboneAtom:
        x1 = atom1.coordinates[0] 
        y1 = atom1.coordinates[1] 
        z1 = atom1.coordinates[2]
        r1 = getRadius(atom1)
        for atom2 in sidechainAtom:
            x2 = atom2.coordinates[0] 
            y2 = atom2.coordinates[1] 
            z2 = atom2.coordinates[2]
            r2 = getRadius(atom2)
            distance = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            minDistance = r1 + r2
            if distance < minDistance:
                print(atom1.atomName, atom1.resSeq, atom1.coordinates, atom2.atomName, atom2.resSeq, atom2.coordinates)
                return True
    return False#}}}

def pegClashes(protAtom, pegs, testAtom, numPEGadded, atomNames, resSeq=-1):#{{{
# Determines if the peg (testAtom) in its current position clashes with
# any other atoms in the protein, based on the atomic radius of each atom.
    atom1_exclusions = {atomNames["pegBondAtomFront"], atomNames["pegRemoveAtomFront"]}
    if numPEGadded == 0:
        atom2_exclusions = {atomNames["protBondAtom"], atomNames["protRemoveAtom"]}
    else:
        atom2_exclusions = {atomNames["pegBondAtomBack"], atomNames["pegRemoveAtomBack"]}
    for atom1 in testAtom:
        if not atom1.atomName in atom1_exclusions:
            r1 = getRadius(atom1)
            for i in range(numPEGadded):
                for atom2 in pegs[i].atoms:
                    if not (atom2.atomName in atom2_exclusions and i == numPEGadded-1):
                        r2 = getRadius(atom2)
                        vect = atom1.coordinates - atom2.coordinates
                        distance = getMagVector(vect)
                        minDistance = r1 + r2
                        if distance < minDistance:
                            return True
            for atom2 in protAtom:
                if not (atom2.atomName in atom2_exclusions and resSeq==atom2.resSeq):
                    r2 = getRadius(atom2)
                    vect = atom1.coordinates - atom2.coordinates
                    distance = getMagVector(vect)
                    minDistance = r1 + r2
                    if distance < minDistance:
                        return True
    return False#}}}

def modifyAtomProperties(allAtom, resAtom, resSeq, resName):#{{{
# Changes resSeq, resType, chainType, resNum, o, and tf of atoms
# in the protein and in the mutated residue appropriately to 
# enable creation of a new PDB file. Should be called after the
# residue is positioned correctly.
    template_atom_name = ' N  '
    for atom in resAtom:
        atom.resSeq = resSeq
    for atom in allAtom:
       if atom.resSeq == resSeq and atom.atomName == template_atom_name:
           template_atom = atom
    chainID = template_atom.chainType
    resNumID = template_atom.resNum
    for atom in allAtom:
        if atom.resSeq == resSeq:
            atom.resType = resName
            atom.chainType = chainID
            atom.resSeq = resSeq
            atom.resNum = resNumID
            atom.o = 1.0
            atom.tf = 0.0#}}}

def renumberAtoms(atoms):#{{{
# All atoms in the protein are renumbered from 1 to account for
# the mutation of a residue in the protein.
    for i, atom in enumerate(atoms):
        atom.atomNum = i+1
        atom.atomSeq = i+1#}}}

def insertPEGAtoms(pegs, protAtom):#{{{
# Properties resSeq, chainType, and atomSeq are added to all of the PEG
# Atoms, and then the atoms are added into the protein in preparation
# for creating a new PDB file.
    allAtom = copy.copy(protAtom)
    maxResSeq = 0
    maxChainID = 'A'
    maxAtomSeq = 0
    for atom in protAtom:
        if atom.chainType > maxChainID:
            maxChainID = atom.chainType
        if atom.atomSeq > maxAtomSeq:
            maxAtomSeq = atom.atomSeq
    
    resSeq = maxResSeq + 1
    chainID = chr(ord(maxChainID)+1)
    atomSeq = maxAtomSeq + 1

    for peg in pegs:
        for atom in peg.atoms:
            atom.resSeq = resSeq
            atom.chainType = chainID
            atom.atomSeq = atomSeq
            allAtom.append(atom)
            atomSeq += 1
        resSeq += 1
    return allAtom#}}}
#}}}
#########################################################


################## mutateResidue #########################
def mutateResidue(pAtom, rAtom, resSeq, resName, secAtomName=' CB '):#{{{
    # Description: 
    # Mutate an amino acid (typically uAA) into an all-atom protein at a specified position.
    #     WARNING: The structure should be minimized after using this script!!!!!
    #     WARNING: This function does not provide input paramter error checking!!!!!

    # Inputs:
    # pAtom        This is the array of KPI.Atoms for the protein
    # rAtom        This is the array of KPI.Atoms for the amino acid that will be mutated in.
    # resSeq       This is the resSeq of the amino acid in the protein that you would like to
    #                  replace with the amino acid represented by rAtom.
    # resName      This is the name of the residue that is being mutated into the protein, 
    #                  typically a three character residue name, as in a PDB. We often use
    #                  'PAZ'.
    # secAtomName  This is the name of the second atom of the residue that is being mutated
    #                  into the protein. The amino acid is rotated about the bond from CA to
    #                  this atom to find a good fit when mutating the residue in. Typically,
    #                  and by default, this is ' CB '.

    # Actions:
    # Given the proper inputs as mentioned above, will remove the residue specified.
    # Will then insert the residue given (typically an uAA) at the same location.
    # The residue will then be rotated so that it fits without clashing with any other atoms
    #     in the protein Atom array. If there is absolutely no position and rotation that will 
    #     avoid clashes, an error message will be printed, and a protein with the new residue
    #     in a poor position will be generated. It may be possible to use the protein with 
    #     minimization, but more likely is that there was an error in the inputs.
    # Atom properties of both the protein and new residue are modified to create a useable
    #     PDB file.
    # Finally, an array of Atoms with the new residue replacing the old is returned, ready
    # to be output as a PDB file.

    # Outputs:
    # An array of "KPI.Atoms" representing the mutated protein.

    protAtom = copy.deepcopy(pAtom)
    resAtom = copy.deepcopy(rAtom)
    alignResidue(resAtom, protAtom, resSeq)
    removeResidue(protAtom, resSeq)
    insertResidue(resAtom, protAtom, resSeq)
    modifyAtomProperties(protAtom, resAtom, resSeq, resName)
    fitResidue(protAtom, resSeq, secAtomName)
    return protAtom#}}}
#########################################################


################## pegGenAllAtom #########################
def pegGenAllAtom(protAtom, pegAtom, resSeq, atomNames, numPEG):#{{{
    # Description: 
    # Generate a PEG of the specified number of monomers at the specified residue
    #     amd atom of a protein.
    #     WARNING: The structure should be minimized after using this script!!!!!
    #     WARNING: This function does not provide input paramter error checking!!!!!

    # Inputs:
    # protAtom     This is the array of KPI.Atoms for the protein
    # pegAtom      This is the array of KPI.Atoms used as a template for a PEG monomer. 
    # resSeq       This is the resSeq of the amino acid in the protein where you would like
    #                  to add PEG.
    # atomNames    This is a dictionary that gives the names of all the atoms that will be 
    #                  removed or bonded to other atoms as a result of PEGylation.
    #              Needed entries:
    #                  protRemoveAtom         Atom in protein that is removed when first PEG
    #                                             is bonded to it. (Usually ' H* ')
    #                  protBondAtom           Atom in protein that first PEG is bonded to.
    #                                             (Usually ' C* ')
    #                  pegRemoveAtomFront     Atom in PEG that is removed when PEG is bonded
    #                                             to a previous atom. (Usually ' H1C')
    #                  pegRemoveAtomBacki     Atom in PEG that is removed when PEG is bonded
    #                                             to a following atom. (Usually ' H2C')
    #                  pegBondAtomFront       Atom in PEG that is bonded to when PEG is bonded
    #                                             to a previous atom. (Usually ' C1 ')
    #                  pegBondAtomBack        Atom in PEG that is bonded to when PEG is bonded
    #                                             to a following atom. (Usually ' C2 ')
    #                 
    # numPEG       This is the desired length of the PEG chain (number of PEG monomers).

    # Actions:
    # Given the proper inputs, a PEG will be attached at the specified residue.
    # The polymer is created one monomer at a time until the desired length has been achieved.
    # Each monomer is added in a random rotation around the bond formed between it and the 
    #     previous monomer.
    # Numbering and properties of the protein and PEG are modified as needed to produce a valid
    #     PDB format.
    # If a monomer cannot be placed after trying 15 rotations, the previous monomer is removed
    #     and placed in a new rotation, after which PEG extension is again attempted.
    # If a certain monomer is backtracked to 5 times, it is removed and PEG extension begins at
    #     the monomer before it.
    # If there are 5 times as many removals of a placed PEG as there are PEGs to be added, the
    #     script is ended. The PEGylated protein is still returned, but an error message is 
    #     generated, and the PEG will be shorter than specified.
    # Finally, an array of Atoms with the protein and PEG is returned, ready
    #     to be output as a PDB file.

    # Outputs:
    # An array of "KPI.Atoms" representing the PEGylated protein.

    protAtom = copy.deepcopy(protAtom)
    pegAtom = copy.deepcopy(pegAtom)
    numPEGadded = 0
    myPEGs = []
    removals = 0
    maxRemovals = 5*numPEG
    for i in range(numPEG):
        myPEGs.append(Monomer())
    while numPEGadded < numPEG:
        if numPEGadded % 50 == 0:
            print(numPEGadded, "pegs added.")
        if numPEGadded == 0:
            if addMonomer(protAtom, pegAtom, myPEGs, numPEGadded, atomNames, resSeq):
                numPEGadded+=1
            else:
                print("Unable to place first PEG.",
                    "Make sure given residue, atom name, and file are correct and that the residue has surface accessibility")
                quit()
        else:
            if addMonomer(protAtom, pegAtom, myPEGs, numPEGadded, atomNames):
                numPEGadded+=1
            else:
                if removals > maxRemovals:
                    print("Unable to create a sufficiently long polymer.")
                    print("Created a chain ", numPEGadded, " PEGS long.")
                    myPEGs = myPEGs[0:numPEGadded]
                    break
                else:
                    while True:
                        removeMonomer(myPEGs, numPEGadded, protAtom, resSeq)
                        removals += 1
                        numPEGadded -= 1
                        if myPEGs[numPEGadded-1].timesBacktracked < 5 or numPEGadded == 0:
                            break
                        elif myPEGs[numPEGadded-1].timesBacktracked == 5:
                            print("PEG", numPEGadded, "has been backtracked to 5 times.")
                            myPEGs[numPEGadded-1].timesBacktracked = 0
                    
    return insertPEGAtoms(myPEGs, protAtom)#}}}
##########################################################


#-------------- CHARMM FF and data file -----------------
###################### catPRM ###########################
def catPRM(fileNameNew, arrPRM) :#{{{
    # Inputs

    # Actions

    # Outputs

    # Set-up and Pre-define variables
    arrHeaders = ["ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER" ,"CMAP", \
                  "NONBONDED", "NBFIX", "HBOND", "END"] # Not currently being used...
    arrAtom = []
    arrBond = []
    arrAngle = []
    arrDihedral = []
    arrImproper = []
    arrCMAP = []
    arrNonbonded = []
    nonbondedHeader = []
    
    def readPRM(fileName) :#{{{
        # fields: Atom   Bond   Angle  Dihedral Improper CMAP   Nonbonded
        active = [False, False, False, False,   False,   False, False]
        readFile = open(fileName, 'r')
        data = readFile.readlines()
        readFile.close()
        for i, line in enumerate(data) :
            if re.match("ATOMS", line) :
                active = [True, False, False, False, False, False, False]
                arrAtom.append("!atom data from %s\n" % (fileName))
                continue
            if re.match("BONDS", line) :
                active = [False, True, False, False, False, False, False]
                arrBond.append("!bond data from %s\n" % (fileName))
                continue
            if re.match("ANGLES", line) :
                active = [False, False, True, False, False, False, False]
                arrAngle.append("!angle data from %s\n" % (fileName))
                continue
            if re.match("DIHEDRALS", line) :
                active = [False, False, False, True, False, False, False]
                arrDihedral.append("!dihedral data from %s\n" % (fileName))
                continue
            if re.match("IMPROPER", line) :
                active = [False, False, False, False, True, False, False]
                arrImproper.append("!improper data from %s\n" % (fileName))
                continue
            if re.match("CMAP", line) :
                active = [False, False, False, False, False, True, False]
                arrCMAP.append("!CMAP data from %s\n" % (fileName))
                continue
            if re.match("NONBONDED", line) :
                active = [False, False, False, False, False, False, True]
                arrNonbonded.append("!Nonbonded data from %s\n" % (fileName))
                nonbondedHeader.append(data[i] + data[i+1])
                continue
            if re.match("END", line) or re.match("NBFIX", line) or re.match("HBOND", line):
                active = [False, False, False, False, False, False, False]
                continue
            # ADD NBFIX functionality

            # Add lines to appropriate fields
            if active[0] == True :
                arrAtom.append(line)
            if active[1] == True :
                arrBond.append(line)
            if active[2] == True :
                arrAngle.append(line)
            if active[3] == True :
                arrDihedral.append(line)
            if active[4] == True :
                arrImproper.append(line)
            if active[5] == True :
                arrCMAP.append(line)
            if active[6] == True and not re.match("cutnb", line):
                arrNonbonded.append(line)
    # End Function#}}}

    def writePRM() :#{{{
        writeFile = open(fileNameNew, 'w')
        writeFile.write("*>>>>>>>>>>>>> Concatenated Parameter File <<<<<<<<<<<<<\n")
        writeFile.write("*>>> The .prm files that are contained in this file: <<<\n")
        for file in arrPRM : writeFile.write("*>>> %s\n" % (file))
        
        writeFile.write("\nATOMS\n")
        for line in arrAtom : writeFile.write(line)

        writeFile.write("BONDS\n")
        for line in arrBond : writeFile.write(line)

        writeFile.write("ANGLES\n")
        for line in arrAngle : writeFile.write(line)

        writeFile.write("DIHEDRALS\n")
        for line in arrDihedral : writeFile.write(line)

        writeFile.write("IMPROPER\n")
        for line in arrImproper : writeFile.write(line)

        writeFile.write("CMAP\n")
        for line in arrCMAP : writeFile.write(line)

        writeFile.write(nonbondedHeader[0] + "\n")
        for line in arrNonbonded : writeFile.write(line)
        
        writeFile.write("END")
        writeFile.close()
    #End Function#}}}

    def nonbondedHeaderDiff() :#{{{
        check = False
        for header in nonbondedHeader :
            if set(nonbondedHeader[0]) - set(header) : check=True

        if check == True :
            print(">>>>>>>> CAUTION <<<<<<<<<")
            print("Differences in NONBONDED headers between the .prm files exist!")
            print("The header in the first given .prm file (%s) will be used!!" % (arrPRM[0]))
            print("\nThe list of headers:")
            for i, header in enumerate(nonbondedHeader) :
                print("\nFrom file: %s\n%s" % (arrPRM[i], header.rstrip()))
            print(">>>>>>>>>>>>><<<<<<<<<<<<<")
    #End Function#}}}

    # Read-in data from all .prm files
    print("Reading .prm files...")
    for file in arrPRM : readPRM(file)

    # Check data for errors
    print("Checking for consistency...")
    nonbondedHeaderDiff()
    # ADD checks to make sure there aren't any doubles in the data

    # Write data to fileNameNew
    print("Writing data to %s\n\n" % (fileNameNew))
    writePRM()

    return fileNameNew#}}}
#########################################################


###################### catRTF ###########################
def catRTF(fileNameNew, arrRTF, Other=False, ResiCaution=300) :#{{{
    # Inputs

    # Actions
    #-This script assumes that once you encounter the first RESI the rest of the
    #.rtf file only contains residues. It will throw a caution if the residue
    #exceeds 100 lines.
    #-Other stuff (non-mass and non-Resi data) by default will not be included in the 
    #new concatenated file. By default it will be printed to screen.
    #-Note: This this script will not re-number your atom type MASS number. If you have
    #repeats or a value of "-1" it will remain as such. Such .rtf files are not
    #compatible with charmm2lammps

    # Outputs

    #Set-up and Pre-define variables
    arrNewRESI = ["RESI", "!RESI", "PRES"] # These key words indicate a new residue is being encountered
    arrMasses = []  # arrMasses includes lines that start with MASS
    arrRESI = []    # arrRESI includes residues with the headers: RESI, !RESI and PRES
    arrOther = []   # arrOther records all other lines, but ignores lines that start with !, [0-9], END, or *
    
    def readRTF(fileName) :#{{{
        arrRESI.append("!Residue data from %s\n" % (fileName))
        arrMasses.append("!Mass data from %s" % (fileName)) # DO NOT add newline! "\n" removes compatibility with charmm2lammps.pl
        arrOther.append("\n!Other data from %s\n" % (fileName))
        active = False # Activates RESI read
        oneRESI = []   # Holds the lines of the current RESI being read, used as a check to make sure 
                       # only residue data is being read.

        # Open file and read data
        readFile = open(fileName, 'r')
        data = readFile.readlines()
        readFile.close()

        # Parce data into correct fields
        for line in data :
            line = line.rstrip()
            if active == True and line == '' : 
                arrRESI.append("\n")
                continue
            elif line == '' : continue

            if re.match("END", line) : 
                if len(oneRESI) > ResiCaution : print("CAUTION: %s exceeds %i lines and could contain errors" % \
                                                     (oneRESI[0].split()[1], ResiCaution))
                oneRESI = []
                return
            elif line.split()[0] == "MASS" : 
                arrMasses.append(line)
                continue
            elif active == False and line.split()[0] in arrNewRESI : 
                active = True
                arrRESI.append(line)
                oneRESI.append(line)
                continue
            elif active == False and not re.match("!", line) and not re.match("[0-9]", line) and not re.match("\*",line) : 
                arrOther.append(line)
                continue
            elif active == True and line.split()[0] in arrNewRESI :
                if len(oneRESI) > ResiCaution : print("CAUTION: RESI %s exceeds %i lines and could contain errors" % \
                                                     (oneRESI[0].split()[1], ResiCaution))
                oneRESI = []
                arrRESI.append(line)
                oneRESI.append(line)
                continue
            elif active == True :
                arrRESI.append(line)
                oneRESI.append(line)
    #End Function#}}}

    def writeRTF():#{{{
        writeFile = open(fileNameNew, 'w')
        writeFile.write("*>>>>>>>>>>>>> Concatenated Parameter File <<<<<<<<<<<<<\n")
        writeFile.write("*>>> The .prm files that are contained in this file: <<<\n")
        for file in arrRTF : writeFile.write("*>>> %s\n" % (file))
        writeFile.write("36  1\n\n")
        
        for line in arrMasses : writeFile.write(line + "\n")
        if Other == True : 
            for line in arrOther : writeFile.write(line + "\n")
        for line in arrRESI : writeFile.write(line + "\n")
        writeFile.write("END")
        writeFile.close()
        return
    #End Function#}}}

    # Read-in data from all .rtf files
    print("Reading .rtf files...")
    for file in arrRTF : readRTF(file)

    # Check data and for errors
    
    # Print Other data
    if Other == False :
        print(">>>>>>>>>>>>><<<<<<<<<<<<<")
        print("The following data was not included:")
        for line in arrOther :
            print(line)
        print(">>>>>>>>>>>>><<<<<<<<<<<<<")

    # Write data to fileNameNew
    print("Writing data to %s\n\n" % (fileNameNew))
    writeRTF()

    return fileNameNew#}}}
#########################################################


################## charmmPDBprep ########################
def charmmPDBprep(fileNameRead, fileNameNew, HIS="HSD") :#{{{
    # Inputs:
    # fileNameRead   This is the .pdb file that you got from the protein databank website
    # fileNameNew    This is the name of the CHARMM-compatiable file you want to write

    # Actions:
    # Will read the .pdb and then will write a .pdb file that has removed all the experimental comments
    # and has reset the numbering. 
    # CAUTION: When publishing make sure you're published atoms and residues are appropriate

    # Outputs:
    # There are no outputs, though fileNameNew is created

    # Read pdb
    arrAtom = KPR.read_pdb(fileNameRead)

    for atom in arrAtom :
        # Adjust all the HIS to HSD
        if atom.resType == "HIS" and HIS == "HSD" : atom.resType = "HSD"

        # Sometimes ILE's CD has improper nameing
        if atom.resType == "ILE" and atom.atomName == " CD1" : atom.atomName = " CD "


    errorCheck_SEQREStoATOM(fileNameRead)
    errorCheck_ATOM(fileNameRead)
 
    # Write CHARMM-compatibale pdb file
    KPW.write_pdb(fileNameNew, arrAtom, reset=True)

    return#}}}
#########################################################


############### lammpsRepX_inputWorld ###################
def lammpsRepX_inputWorld(startTemp, endTemp, delta) :#{{{
    fileName = "RepX_worldVarInputs.txt"
    diff = endTemp - startTemp
    print ("\nCalculating appropriate box range...")
    if diff % delta != 0 : 
        print("Changing your max temp to one that equally divides")
        endTemp = endTemp - (diff % delta)
        diff = endTemp - startTemp
        numBoxes = int(diff / delta) + 1
        print("New max temp: %i" % endTemp)
    else : numBoxes = int(diff / delta) + 1
    
    tempArr = np.zeros(numBoxes)
    boxArr  = np.zeros(numBoxes)
    
    for i in range(numBoxes) :
        tempArr[i] = startTemp + i * delta
        boxArr[i]  = i
    
    file = open(fileName, "w")
    pntStr =  'variable          temp    world &\n'
    pntStr += '                                '
    for i, temp in enumerate(tempArr) :
        pntStr += str(int(temp)).ljust(3) + ' '
        if i+1 == len(tempArr) : pntStr += '\n'
        elif (i+1)%10 == 0 : pntStr += '&\n                                '
    file.write(pntStr)
    
    pntStr =  'variable          box     world &\n'
    pntStr += '                                '
    for i, box in enumerate(boxArr) :
        pntStr += str(int(box)).ljust(3) + ' ' 
        if i+1 == len(tempArr) : pntStr += '\n'
        elif (i+1)%10 == 0 : pntStr += '&\n                                '
    file.write(pntStr)
    file.close()
    print("The file %s was created\n" % fileName)
#}}}
#########################################################


#------------------- Error Checks -----------------------
############## errorCheck_SEQREStoATOM ##################
def errorCheck_SEQREStoATOM(fileName) : #{{{
    # This check is to compare a .pdb file that has been downloaded from the internet. It will check sequence consistiency between the
    # the SEQRES section to the ATOM section. The ATOM section is considered more correct.
    # Inputs:
    # fileName        This is the file name of you .pdb file you got from the internet

    # Actions :
    # Will read the pdb and check if there are insertions to homology sequence (ex. 49 50 50A 51...)
    # Will directly compre the sequence in ATOM and the sequence in SEQRES to check if they are identical

    # Outputs:
    # Nothing is outputted, will print to screen if anything is wrong

    arrSEQRES = [] # residue sequence as devined in SEQRES section
    arrPDBRES = [] # residue type in the pdb
    arrPDBCNT = [] # residue count in the PDB

    PDBfile = open(fileName, 'r')
    for line in PDBfile :
        comment = line[0:6]
        if comment == "SEQRES" : 
            line = line[19:].split()
            for seq in line : 
                arrSEQRES.append(seq)

    print("SEQRES section contains %i residues" % len(arrSEQRES))
    PDBfile.close()

    PDBfile = open(fileName, 'r')
    resCount = 0
    additionCheck = ""  # All addition variables are related to insertions into the primary sequence standard (fairly rare)
    additionStr = ""
    for line in PDBfile :
        comment = line[0:6]
        if comment == "ATOM  " :
            resType = line[17:20]
            resSeq  = int(line[22:26])
            addition = line[26]  # If there is a sequence addition into the PDB standard, this will be a letter A-Z
            if addition in KPI.alphabet : additionStr = str(resSeq) + addition
            else : additionStr = ""
            if additionStr == "" and resCount < resSeq :
                arrPDBRES.append(resType)
                arrPDBCNT.append(resSeq)
                resCount = resSeq
            elif additionStr != additionCheck and additionStr != "" and resCount == resSeq : 
                arrPDBRES.append(resType)
                arrPDBCNT.append(resSeq)
                additionCheck = additionStr
                print("\nWARNING: Duplicate residue sequence numbers indicated!!!")
                print("Error occured on sequence number: %s\n" % (additionStr))

    print("ATOM section contains %i residues" % len(arrPDBRES))
    PDBfile.close()

    if   len(arrSEQRES) > len(arrPDBRES) : maxRES = len(arrSEQRES)
    elif len(arrPDBRES) > len(arrSEQRES) : maxRES = len(arrPDBRES)
    else : 
        print("The lengths are equal.")
        return

    for i in range(maxRES) :
        if arrSEQRES[i].strip() != arrPDBRES[i].strip() :
            print("\nERROR on residue %i" % arrPDBCNT[i])
            print("SEQRES: %s;      the 5-AA sequence series: %s %s %s %s %s" % 
                  (arrSEQRES[i], arrSEQRES[i-2], arrSEQRES[i-1], arrSEQRES[i], arrSEQRES[i+1], arrSEQRES[i+2]))
            print("PDB residue: %s; the 5-AA sequence series: %s %s %s %s %s" % 
                  (arrPDBRES[i], arrPDBRES[i-2], arrPDBRES[i-1], arrPDBRES[i], arrPDBRES[i+1], arrPDBRES[i+2]))
            quit()
        if i == maxRES-1 : print("Got to the end of the list...wellll now thats a pickle")
#}}}
#########################################################


################## errorCheck_ATOM ######################
def errorCheck_ATOM(fileName) : #{{{
    pdbSeq = 1

    PDBfile = open(fileName, 'r')
    for line in PDBfile :
        comment = line[0:6]
        if comment == "ATOM  " :
            resSeq = int(line[22:26])
            if pdbSeq == resSeq : continue
            elif resSeq == pdbSeq + 1 : pdbSeq += 1
            elif resSeq > pdbSeq + 1 :
                print("Numbering skip was found on residue %i" % resSeq)
                pdbSeq = pdbSeq + (resSeq-pdbSeq)
            else : 
                print("ERROR: numbering is off count on residue %i" % resSeq)
                pdbSeq = pdbSeq + (resSeq-pdbSeq)
    PDBfile.close()
#}}}
#########################################################


#--------- Prep Scripts for Post-processing -------------
################# mbar prep fuctions ####################
# Functions#{{{
def glob_traj() :#{{{
    # This function globs the trajectory files in the directory
    # It ignores files ending in "vid.dcd" because these are used
    # for videos, rather than nc analysis.
    files = gb.glob('*.dcd')
    correctFiles = []
    
    for file in files:
        if not re.search('.*.vid.dcd', file):
            correctFiles.append(file)
    if(len(correctFiles) == 0):
        print("No valid dcd files in this directory.")
        quit()
    return sorted(correctFiles, key=penult_elem)#}}}

def glob_NCmod() :#{{{
    # This function finds the native contact-specifying file(s) to use.
    # It sorts multiple files alphabetically.
    mod_glob = sorted(gb.glob('../*.nc.mod'))
    if(len(mod_glob) == 0):
        print("No valid mod files in this directory.")
        quit()
    return mod_glob #}}}

def glob_screen() :#{{{
    # This function globs all screen files in the directory
    screen_glob = sorted(gb.glob('screen.*'), key=second_elem)
    if(len(screen_glob) == 0):
        print("No valid screen files in this directory.")
        quit()
    return screen_glob#}}}

def glob_box() :#{{{
    # This function globs all box files in the directory
    box_glob = sorted(gb.glob('box.*'), key=second_elem)
    if(len(box_glob) == 0):
        print("No valid box files in this directory.")
        quit()
    return box_glob#}}}

def glob_nc() :#{{{
    # This function globs all nc files in the directory
    nc_glob = sorted(gb.glob('nc.*'), key=second_elem)
    if(len(nc_glob) == 0):
        print("No valid nc files in this directory.")
        quit()
    return nc_glob#}}}

def read_native_contacts(nc_files, nc_array, distance_multiplier) :#{{{
    # This function reads the native contact data from the nc mod file(s) into an array.
    # Instead of storing sigma, (distance_multiplier*sigma)**2 is stored for
    # efficiency. The order of the contacts is the order found within the files, and the files
    # are ordered alphabetically.
    linePattern = "pair_coeff\s*" + KPR.intPat + KPR.intPat + KPR.floatPat + KPR.floatPat
    lines = []
    for filename in nc_files:
        file = open(filename)
        for line in file:
            if re.match(linePattern, line):
                lines.append(line)
        file.close()

    for line in lines:
        #Minus 1 on this line and the next to account for starting at index 0 instead of 1
        atom1 = int(line.split()[1]) - 1
        atom2 = int(line.split()[2]) - 1
        if atom1 > atom2 :
            atom1, atom2 = atom2, atom1
        extended_sigma_sq = (distance_multiplier*float(line.split()[4]))**2
        nc_array.append((atom1, atom2, extended_sigma_sq))#}}}

def contact_made(u, contact):#{{{
    # This function determines whether the supplied contact is made for the current timestep
    # by comparing the square of the distance between the two atoms with the square of the 
    # stored extended distance.
    atom1 = contact[0]
    atom2 = contact[1]
    extended_sigma_sq = contact[2]
    atom1_pos = u.atoms[atom1].position
    atom2_pos = u.atoms[atom2].position
    x_distance = atom1_pos[0] - atom2_pos[0]
    y_distance = atom1_pos[1] - atom2_pos[1]
    z_distance = atom1_pos[2] - atom2_pos[2]

    distance_sq = x_distance*x_distance + y_distance*y_distance + z_distance*z_distance 
    
    if distance_sq <= extended_sigma_sq:
        return True
    else:
        return False#}}}

def read_configuration_data(config_file) :#{{{
    # Function used to read data from the file that contains information about when the boxes switch temperatures
    file = open(config_file)
    for i, line in enumerate(file):
        if re.search('Step', line):
            headerRows=i+1
            break
    file.close()
    configurations = np.loadtxt(config_file, dtype=int, skiprows=headerRows) 
    return configurations #}}}

def check_user_inputs(top_file, config_file, sim_timestep) :#{{{
    # Function used to check validity of user input files
      if not os.path.isfile(top_file) :
          print("User specified top_file does not exist.")
          quit()
      if not os.path.isfile(config_file) :
          print("User specified config_file does not exist.")
          quit()
      if sim_timestep <= 0 :
          sim_timestep = 1
          print("User specified sim_timestep not possible.")
          print("Changing value to default: sim_timestep = 0")#}}}

def filter_lines(currFileName, beginPattern, endPattern, nthMatch) :#{{{
    # This function provides the regular expression (re) file search
    # It returns a sub-list containing lines from currFileName between beginPattern and endPattern
    # It uses the 'nthMatch' line that matches beginPattern as the starting point
    lines = []
    file = open(currFileName)
    foundBeginPattern = 0
    
    for line in file:
        if not foundBeginPattern < nthMatch  and endPattern != '' and re.search(endPattern, line):
            break
        elif not foundBeginPattern < nthMatch:
            lines.append(line)
        elif foundBeginPattern < nthMatch  and re.search(beginPattern, line):
            foundBeginPattern+=1
    
    file.close()
    return lines#}}}

def limData(maxRows, arrVal) :#{{{
    while len(arrVal[0]) > maxRows :
        check = len(arrVal[0])
        if check%1000 == 0 : print("Timesteps still in excess; %i points long" % (check))
        removeLine = np.random.randint(len(arrVal[0])-1)
        for i in range(len(arrVal)) : arrVal[i] = np.delete(arrVal[i], removeLine, 0)
    return arrVal#}}}
#}}}
#########################################################


############## mbarNC_screenDCD_to_nc ###################
def mbarNC_screenDCD_to_nc (top_file, sim_timestep, config_file, distance_multiplier=1.2) :#{{{
    # Typical values for inputs: 
    # top_file = ../protein.data
    # sim_timestep = 3
    # config_file = screen
    # distance_multiplier = 1.2 # Multiplier used to create the neighborhood list of "made" native contacts (1.2 is std)
    
    # Output:
    # Text files containing a string of 1's and 0's for each timestep.
    # The files will be called nc.[temperature]
    # 1 means the contact was made, 0 means it wasn't.
    # Order of contacts is the same as that of the param file in the directory.
    
    import MDAnalysis as mda
    
    print("\nBegining native contact analysis using .dcd files...")
    ############## Vars and Checks #################
    # Set-up and Check inputs
    nc_array = []                 # Array where native contact info will be stored
    screen_array = []             # Array where the native contact data is stored before being sorted
    out_array = []                # Array where the output is stored before writing it to a file
    stepNums = []                 # Array that holds the steps where nc data was obtained
    
    check_user_inputs(top_file, config_file, sim_timestep)
    
    # Define the trajectory, topology, and native contact files
    traj_files = glob_traj()            # Glob of all dcd files in the directory
    nc_files = glob_NCmod()             # Glob of nc mod file(s)
    screen_files = glob_screen()        # Glob of all screen files in the directory
    numBoxes = len(traj_files)          # Number of temperature boxes used
    if numBoxes == 0 : 
        print("There are no .dcd files in your directory")
        quit()
    
    ################### Main ########################
    # Fill native contact array and create a Universe from the given files
    print("Reading configuration data...")
    configurations = read_configuration_data(config_file)
    print("Reading native contact data...")
    read_native_contacts(nc_files, nc_array, distance_multiplier)
    
    # Initializing
    ts_conversion = KPI.tsConv * sim_timestep #Factor to convert time stored in dcd files to the timesteps used in the simulation
    tempRange = np.zeros(numBoxes)
    startProduction = 0.0 
    foundInitTemp = False
    foundProdStep = False
    
    # Defining simulation parameters
    print("Defining simulation parameters...")
    for i, file in enumerate(screen_files):
        foundInitTemp = False
        open_file = open(file)
        for line in open_file:
            if re.match('Initial screen temperature', line) :
                tempRange[i] = int(line.split()[3])
                foundInitTemp = True
            if re.match('Start of production phase on step', line) :
                startProduction = int(line.split()[6])
                foundProdStep = True
            if foundInitTemp == True and foundProdStep == True : break
         
        # Checks to see if constants were labeled in the simulation...
        if foundInitTemp == False :
            tempRange[i] = i 
            print('No initial temperature found. Numbering box numerically. WARNING: temperature mbar input file incorrect.')
        if foundProdStep == False :
            startProduction = 0 
            print ('WARNING: No indication of production step mentioned, assuming equilibration phase of 0 steps.')
        open_file.close()
    
    # Run MDAnalysis
    print("Running MDAnalysis...")
    # for i, t in enumerate(traj_files):
      #  unis.append(mda.Universe(top_file, t))
      #  screen_array.append([])
      #  out_array.append('')
    
    # For each timestep, test whether each native contacts is made. Store data to out_array.
    print("Computing number of made contacts at each temperature...")
    for i, t in enumerate(traj_files):
        print("Working on file %s..." % traj_files[i])
        screen_array.append([])
        out_array.append('')
        u = mda.Universe(top_file, t) 
        for ts in u.trajectory:
            # Read in the timestep numbers from the first universe
            if i == 0:
                stepNums.append(int(round(ts.time/ts_conversion)))
            out_str = str(int(round(ts.time/ts_conversion))).ljust(14) # Starts out string with step number
            nc_str = ""
            contacts_made = 0
            for j in range (len(nc_array)):
                if contact_made(u, nc_array[j]):
                    nc_str += "1"
                    contacts_made += 1
                else:
                    nc_str += "0"
            out_str += str(contacts_made).ljust(11) 
            out_str += nc_str
            out_str += "\n"
            screen_array[i].append(out_str)
    
    # Find  number of steps betwen swaps
    swapEvery = configurations[1][0]
    
    # Write data from out_array to file            
    # Append lines from each screen to the proper string of boxLines for each swap,
    # and continue until the last line of configurations
    print("Sorting data...")
    stepIndex=0
    while stepNums[stepIndex]<startProduction:
        stepIndex+=1 
    stop=len(configurations)-1
    if startProduction == 0:
        i=0 
    else:
        i=int((startProduction-1)/swapEvery)   #The -1 is because of the point when LAMMPS swaps boxes

    while i<stop:
        if i%1000==0 : print('Sorting step ', i*swapEvery)
        for j in range(numBoxes):
            screenIndex=configurations[i][j+1]
            k=stepIndex
            while k<len(stepNums) and stepNums[k]<=(i+1)*swapEvery:
                out_array[j]+=screen_array[screenIndex][k]
                k+=1
        i+=1    
        stepIndex=k
    
    print("Writing output files...")
    for i in range(numBoxes):
        with open("nc." + str(int(tempRange[i])), 'w') as writeFile:
            header_string = "# File \"" + "nc." + str(int(tempRange[i])) + "\". "
            header_string +=  "Created " + str(datetime.now()) + ".\n"
            header_string += "# Total native contacts: " + str(len(nc_array))  + "\n"
            writeFile.write(header_string)
            writeFile.write(out_array[i])
    
    for i in range(numBoxes):
        with open("nc_screen." + str(i), 'w') as writeFile:
            header_string = "# File \"" + "nc_screen." + str(i) + "\". "
            header_string +=  "Created " + str(datetime.now()) + ".\n"
            header_string += "# Total native contacts: " + str(len(nc_array))  + "\n"
            writeFile.write(header_string)
            for line in screen_array[i]:
                    writeFile.write(line)
    print("Native contact analysis complete!\n")#}}}
#########################################################


################ mbarNC_screen_to_box ###################
def mbarEner_screen_to_box() :#{{{
    # Function to filter and sort screen file information into box files
    # Read-in the data lines from "log.lammps" or "screen" if "log.lammps" doesn't exist
    # and create a sorted glob to manage "screen.#" files 
    if   os.path.isfile('./log.lammps') :       currFileName='log.lammps'
    elif os.path.isfile('./log.setup.lammps') : currFileName='log.setup.lammps'
    elif os.path.isfile('screen') :             currFileName='screen'
    else: print('No Replica Exchange key-file was found')
    configurations = read_configuration_data(currFileName)
    screenFiles = glob_screen()
    numBoxes = len(screenFiles)

    #Extract header line from screen.0, and create formatted header line to be added to all box files
    file = open(screenFiles[0])
    for line in file:
        if re.search('TotEng', line) and re.search('Temp', line) and re.search('Step', line):
            headerLine=line
            formattedHeaderLine=line.split()[0].rjust(8)
            for i in range(len(line.split())-1):
                formattedHeaderLine= ''.join((formattedHeaderLine, line.split()[i+1].rjust(13)))
            formattedHeaderLine+='\n'
    headerMatches = 0 #how many lines contain header info. Only data beneath the last header will be used
    file.seek(0)
    for line in file:
        if re.search(headerLine, line):
            headerMatches+=1
    file.close()

    # Read-in the data lines from "screen.#" files and store them as a list of lists
    screenLines=[]
    tempRange = np.zeros(numBoxes)
    startProduction = 0.0
    foundInitTemp = False
    foundProdStep = False
    for i, glob in enumerate(screenFiles):
        foundInitTemp = False
        file = open(glob)
        for line in file:
            if re.search('Initial screen temperature', line) :
                tempRange[i] = int(line.split()[3])
                foundInitTemp = True
            if re.search('Start of production phase on step', line) :
                startProduction = int(line.split()[6])
                foundProdStep = True
            if foundInitTemp == True and foundProdStep == True : break
        
        # Checks to see if constants were labeled in the simulation...
        if foundInitTemp == False :
            tempRange[i] = i
            print('No initial temperature found. Numbering box numerically. WARNING: temperature mbar input file incorrect.')
        if foundProdStep == False :
            startProduction = 0
            print ('WARNING: No indication of production step mentioned, assuming equilibration phase of 0 steps.')

        screenLines.append(filter_lines(glob, headerLine, 'Loop time of', headerMatches))
        file.close()

    # Read-in the thermo step numbers from screenLines[0] (which comes from "screen.0") and store in a list
    thermoSteps=[]
    for i, line in enumerate(screenLines[0]):
        thermoSteps.append(int(screenLines[0][i].split()[0]))


    # Find number of boxes and  number of steps between swaps
    swapEvery = configurations[1][0] #step number of second configuration line


    # Create a list to store the lines for each box file and write the header line to each
    boxLines=[]
    for i in range(numBoxes):
        boxLines.append(formattedHeaderLine)


    print('\nReading from \"screen\", sorting data and filtering into box files...')

    # Append lines from each screen to the proper string of boxLines for each swap,
    # and continue until the last line of configurations
    thermoIndex=0
    while thermoSteps[thermoIndex]<startProduction:
        thermoIndex+=1 #Added this line and above, may need error checking
    stop=len(configurations)-1
    if startProduction == 0:
        i=0
    else:
        i=int((startProduction-1)/swapEvery)   #The -1 is because of the point when LAMMPS swaps boxes
    while i<stop:
        if i%1000==0 : print ('Sorting step ', i*swapEvery)
        for j in range(numBoxes):
            screenIndex=configurations[i][j+1]
            k=thermoIndex
            while k<len(thermoSteps) and thermoSteps[k]<=(i+1)*swapEvery:
                boxLines[j]+=screenLines[screenIndex][k]
                k+=1
        i+=1    
        thermoIndex=k

    # Create "box.#" files and write data from boxLines to them        
    for i in range(numBoxes):
        with open('box.' + str(int(tempRange[i])), 'w') as currentWriteFile:
            currentWriteFile.write(boxLines[i])#}}}
#########################################################


################ mbarPrep_box_to_mbar ###################
def mbarPrep_box_to_mbar(maxRows):#{{{
    # Read-in data from box.# files and put them in mbar-readable format
    files = glob_box()
    numBoxes   = len(files)

    arrTemp = np.zeros(numBoxes)

    # Grab relevent data from box.* files
    for i, file in enumerate(files):
        posStep   = -1
        posTemp   = -1
        posPotEng = -1
        posTotEng = -1
        boxTemp = int(file.split(".")[1])

        # Print status
        print ("Reading file: box", boxTemp)

        # Search the header to find col positions for important data
        fileTMP = open(file)
        line = fileTMP.readline()
        count = 0
        for j in range(0, len(line.split())) :
            if line.split()[j] == "Step" : posStep = count
            if line.split()[j] == "Temp" : posTemp = count
            if line.split()[j] == "PotEng" : posPotEng = count
            if line.split()[j] == "TotEng" : posTotEng = count
            count += 1
        fileTMP.close()

        # Check to make sure col positions are defined
        if posStep < 0 : print("\"Step\" column undefined")
        if posTemp < 0 : print("\"Temp\" column undefined")
        if posPotEng < 0 : print("\"PotEng\" column undefined")
        if posTotEng < 0 : print("\"TotEng\" column undefined")

        # Read-in and save Temperature, Total Energy and Potential energy
        data = np.loadtxt(file,skiprows=1)
        arrTemp[i] = boxTemp
        col = data[:,posTotEng]
        if i == 0 : arrTE = col.copy()
        if i > 0  : arrTE = np.vstack([arrTE,col])
        col = data[:,posPotEng]
        if i == 0 : arrPE = col.copy()
        if i > 0  : arrPE = np.vstack([arrPE,col])
    #------

    # Limit total number of data points
    print("Decreasing the number of data proints to %i..." % (maxRows))
    arrPET = arrPE.T
    arrTET = arrTE.T
    if len(arrPET) != len(arrTET) : 
        print("\nERROR: Your PE and TE data do not match!!!\n")
        quit()
    arrVal = [arrPET, arrTET]
    limData(maxRows, arrVal)
    arrPE = arrVal[0].T
    arrTE = arrVal[1].T

    currDirectory = "mbar_input"

    # Write Directory
    if not os.path.exists(currDirectory) : os.makedirs(currDirectory)
    # Write Temperatures
    pntFile = open(currDirectory + '/temperatures', 'w')
    pntStr = ''.join('{:7}'.format(val) for val in arrTemp)
    pntFile.write(pntStr)
    pntFile.close()
    # Write Total Energies
    pntFile = open(currDirectory + '/total-energies', 'w')
    pntStr = ''
    for row in arrTE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()
    # Write Total Energies
    pntFile = open(currDirectory + '/potential-energies', 'w')
    pntStr = ''
    for row in arrPE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()

    # Print status
    for j in range(numBoxes) :
        print ("Directory" , currDirectory, "contains PE, TE and Temperature mbar inputs for box", int(arrTemp[j]))
    print("MBAR preparations are finished. Ready for MBAR analysis.\n")#}}}
#########################################################


############ mbarPrep_box_and_nc_to_mbar ################
def mbarPrep_box_and_nc_to_mbar(maxRows):#{{{
    # Read-in data from box.# and nc.# files and put them in mbar-readable format
    boxFiles = glob_box()
    ncFiles = glob_nc()
    numBoxes = len(boxFiles)
    if len(ncFiles) != numBoxes : 
        print("Number of native contact files do not match number of box files!")
        quit()

    arrTemp = np.zeros(numBoxes)

    # Grab relevent data from box.* files
    for i, file in enumerate(boxFiles):
        posStep   = -1
        posTemp   = -1
        posPotEng = -1
        posTotEng = -1
        boxTemp = int(file.split(".")[1])

        # Print status
        print ("Reading file: box", boxTemp)

        # Search the header to find col positions for important data
        fileTMP = open(file)
        line = fileTMP.readline()
        count = 0
        for j in range(0, len(line.split())) :
            if line.split()[j] == "Step" : posStep = count
            if line.split()[j] == "Temp" : posTemp = count
            if line.split()[j] == "PotEng" : posPotEng = count
            if line.split()[j] == "TotEng" : posTotEng = count
            count += 1
        fileTMP.close()

        # Check to make sure col positions are defined
        if posStep < 0 : print("\"Step\" column undefined")
        if posTemp < 0 : print("\"Temp\" column undefined")
        if posPotEng < 0 : print("\"PotEng\" column undefined")
        if posTotEng < 0 : print("\"TotEng\" column undefined")

        # Read-in and save Temperature, Total Energy and Potential energy
        data = np.loadtxt(file,skiprows=1)
        arrTemp[i] = boxTemp
        col = data[:,posTotEng]
        if i == 0 : arrTE = col.copy()
        if i > 0  : arrTE = np.vstack([arrTE,col])
        col = data[:,posPotEng]
        if i == 0 : arrPE = col.copy()
        if i > 0  : arrPE = np.vstack([arrPE,col])
    #------

    # Grab relevent data from nc.* files
    for i, file in enumerate(ncFiles) :
        boxTemp = int(file.split(".")[1])
        
        # Print status
        print ("Reading file: nc", boxTemp)

        # Read-in and save total contats made column
        data = np.loadtxt(file, skiprows=2)
        col = data[:, 1]
        if i == 0 : arrNC = col.copy()
        if i > 0  : arrNC = np.vstack([arrNC,col])
    #------

    # Limit total number of data points
    print("\nDecreasing the number of timesteps to %i..." % (maxRows))
    arrPET = arrPE.T
    arrTET = arrTE.T
    arrNCT = arrNC.T
    if len(arrPET) != len(arrTET) : 
        print("\nERROR: Your PE and TE timesteps do not match!!!\n")
        quit()
    elif len(arrPET) != len(arrNCT) : 
        print("\nWARNING: Your Energetic and NC timesteps do not match.\n")
        arrVal1 = [arrNCT]
        limData(maxRows, arrVal1)
        arrVal2 = [arrPET, arrTET]
        limData(maxRows, arrVal2)
        arrVal = [arrVal2[0], arrVal2[1], arrVal1[0]]
    else:
        print("Number of PE, TE and NC timesteps are equal.")
        arrVal = [arrPET, arrTET, arrNCT]
        limData(maxRows, arrVal)
    arrPE = arrVal[0].T
    arrTE = arrVal[1].T
    arrNC = arrVal[2].T

    currDirectory = "mbar_input"

    # Write Directory
    if not os.path.exists(currDirectory) : os.makedirs(currDirectory)
    # Write Temperatures
    pntFile = open(currDirectory + '/temperatures', 'w')
    pntStr = ''.join('{:7}'.format(val) for val in arrTemp)
    pntFile.write(pntStr)
    pntFile.close()
    # Write Total Energies
    pntFile = open(currDirectory + '/total-energies', 'w')
    pntStr = ''
    for row in arrTE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()
    # Write Potential Energies
    pntFile = open(currDirectory + '/potential-energies', 'w')
    pntStr = ''
    for row in arrPE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()
    # Write Native Contacts
    pntFile = open(currDirectory + '/contact_info', 'w')
    pntStr = ''
    for row in arrNC.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()

    # Print status
    for j in range(numBoxes) :
        print ("Directory" , currDirectory, "contains PE, TE, contact and Temperature mbar inputs for box", int(arrTemp[j]))
    print("MBAR preparations are finished. Ready for MBAR analysis.\n")#}}}
#########################################################
    

######## mbarPrep_box_and_nc_to_mbar_domain #############
def mbarPrep_box_and_nc_to_mbar_domain(maxRows,domainAA,NCmodFile):#{{{
    # Read-in data from box.# and nc.# files and put them in mbar-readable format
    boxFiles = glob_box()
    ncFiles = glob_nc()
    numBoxes = len(boxFiles)
    if len(ncFiles) != numBoxes : 
        print("Number of native contact files do not match number of box files!")
        quit()

    arrTemp = np.zeros(numBoxes)

    # Grab relevent data from box.* files
    for i, file in enumerate(boxFiles):
        posStep   = -1
        posTemp   = -1
        posPotEng = -1
        posTotEng = -1
        boxTemp = int(file.split(".")[1])

        # Print status
        print ("Reading file: box", boxTemp)

        # Search the header to find col positions for important data
        fileTMP = open(file)
        line = fileTMP.readline()
        count = 0
        for j in range(0, len(line.split())) :
            if line.split()[j] == "Step" : posStep = count
            if line.split()[j] == "Temp" : posTemp = count
            if line.split()[j] == "PotEng" : posPotEng = count
            if line.split()[j] == "TotEng" : posTotEng = count
            count += 1
        fileTMP.close()

        # Check to make sure col positions are defined
        if posStep < 0 : print("\"Step\" column undefined")
        if posTemp < 0 : print("\"Temp\" column undefined")
        if posPotEng < 0 : print("\"PotEng\" column undefined")
        if posTotEng < 0 : print("\"TotEng\" column undefined")

        # Read-in and save Temperature, Total Energy and Potential energy
        data = np.loadtxt(file,skiprows=1)
        arrTemp[i] = boxTemp
        col = data[:,posTotEng]
        if i == 0 : arrTE = col.copy()
        if i > 0  : arrTE = np.vstack([arrTE,col])
        col = data[:,posPotEng]
        if i == 0 : arrPE = col.copy()
        if i > 0  : arrPE = np.vstack([arrPE,col])
    #------

    # Grab relevent data from nc.* files 
    # Section defines 'domain' NC data
    def referenceDomains(NCmodFile) :
        readFile = open(NCmodFile)
        refDomainAA = []
        for i,line in enumerate(readFile):
            line = line.split()
            if int(line[1]) in domainAA or int(line[2]) in domainAA :
                refDomainAA.append(i)
        readFile.close()
        return refDomainAA
    def domainPercentNC(refDomaAA,ncFile) :
        readFile = open(ncFile)
        colNC = []
        for line in readFile :
            if line[0] == "#" : continue
            line = line.split()
            line = line[2]
            numerator = 0
            denominator = 0
            for item in refDomainAA :
                refVal = line[item]
                if   refVal == '0' : 
                    denominator += 1
                elif refVal == '1' : 
                    denominator += 1
                    numerator   += 1
                else :
                    print("string in %s contains something other than a 1 or 0: %s" % (ncFile,refVal))
            colNC.append(numerator / denominator)
        readFile.close()
        return colNC

    for i, file in enumerate(ncFiles) :
        boxTemp = int(file.split(".")[1])
        
        # Print status
        print ("Reading file: nc", boxTemp)

        # Read-in and save total contats made column
        #data = np.loadtxt(file, skiprows=2)
        refDomainAA = referenceDomains(NCmodFile)
        col = domainPercentNC(refDomainAA,file)
        if i == 0 : arrNC = col.copy()
        if i > 0  : arrNC = np.vstack([arrNC,col])
    #------

    # Limit total number of data points
    print("\nDecreasing the number of timesteps to %i..." % (maxRows))
    arrPET = arrPE.T
    arrTET = arrTE.T
    arrNCT = arrNC.T
    if len(arrPET) != len(arrTET) : 
        print("\nERROR: Your PE and TE timesteps do not match!!!\n")
        quit()
    elif len(arrPET) != len(arrNCT) : 
        print("\nWARNING: Your Energetic and NC timesteps do not match.\n")
        arrVal1 = [arrNCT]
        limData(maxRows, arrVal1)
        arrVal2 = [arrPET, arrTET]
        limData(maxRows, arrVal2)
        arrVal = [arrVal2[0], arrVal2[1], arrVal1[0]]
    else:
        print("Number of PE, TE and NC timesteps are equal.")
        arrVal = [arrPET, arrTET, arrNCT]
        limData(maxRows, arrVal)
    arrPE = arrVal[0].T
    arrTE = arrVal[1].T
    arrNC = arrVal[2].T

    currDirectory = "mbar_input"

    # Write Directory
    if not os.path.exists(currDirectory) : os.makedirs(currDirectory)
    # Write Temperatures
    pntFile = open(currDirectory + '/temperatures', 'w')
    pntStr = ''.join('{:7}'.format(val) for val in arrTemp)
    pntFile.write(pntStr)
    pntFile.close()
    # Write Total Energies
    pntFile = open(currDirectory + '/total-energies', 'w')
    pntStr = ''
    for row in arrTE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()
    # Write Potential Energies
    pntFile = open(currDirectory + '/potential-energies', 'w')
    pntStr = ''
    for row in arrPE.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()
    # Write Native Contacts
    pntFile = open(currDirectory + '/contact_info', 'w')
    pntStr = ''
    for row in arrNC.T :
        pntRow = ' '.join('{:12}'.format(val) for val in row)
        pntStr += pntRow + '\n'
    pntFile.write(pntStr)
    pntFile.close()

    # Print status
    for j in range(numBoxes) :
        print ("Directory" , currDirectory, "contains PE, TE, contact and Temperature mbar inputs for box", int(arrTemp[j]))
    print("MBAR preparations are finished. Ready for MBAR analysis.\n")#}}}
#########################################################
 

####################### shiftDCD ########################
def shiftDCD(top_file, traj_file, outname, box_dims, x_atom, y_atom, z_atom):#{{{
    #Import statements for MDAnalysis
    import MDAnalysis as mda 
    from MDAnalysis import transformations
    from MDAnalysis.coordinates.LAMMPS import DCDWriter

    #Box dimensions
    x_lo = box_dims[0]
    x_hi = box_dims[1]
    y_lo = box_dims[2]
    y_hi = box_dims[3]
    z_lo = box_dims[4]
    z_hi = box_dims[5]

    #X atom vals
    move_x = x_atom[0]
    x_segname = x_atom[1]
    x_resid = x_atom[2]
    x_atom_name = x_atom[3]
    x_orig = x_atom[4]

    #Y atom vals
    move_y = y_atom[0]
    y_segname = y_atom[1]
    y_resid = y_atom[2]
    y_atom_name = y_atom[3]
    y_orig = y_atom[4]

    #Z atom vals
    move_z = z_atom[0]
    z_segname = z_atom[1]
    z_resid = z_atom[2]
    z_atom_name = z_atom[3]
    z_orig = z_atom[4]

    u = mda.Universe(top_file, traj_file)

    new_traj = []

    for ts in u.trajectory:
        x_shift = 0
        y_shift = 0
        z_shift = 0
        if(move_x):
            x_atom = u.select_atoms("atom" + " " + x_segname + " " + str(x_resid) + " " + x_atom_name)
            x_num = x_atom[0].position[0]
            x_shift = x_orig - x_num
        if(move_y):
            y_atom = u.select_atoms("atom" + " " + y_segname + " " + str(y_resid) + " " + y_atom_name)
            y_num = y_atom[0].position[1]
            y_shift = y_orig - y_num
        if(move_z):
            z_atom = u.select_atoms("atom" + " " + z_segname + " " + str(z_resid) + " " + z_atom_name)
            z_num = z_atom[0].position[2]
            z_shift = z_orig - z_num
        new_traj.append(transformations.translate([x_shift, y_shift, z_shift])(ts))

    for i, ts in enumerate(new_traj):
        if i%25 == 0:
            print("Timestep", i)
        for i in range(len(ts.positions)):
            position = ts.positions[i]
            x = position[0]
            y = position[1]
            z = position[2]
            if x < x_lo:
                x = x_hi + (x-x_lo)
            elif x > x_hi:
                x = x_lo + (x-x_hi)
            if y < y_lo:
                y = y_hi + (y-y_lo)
            elif y > y_hi:
                y = y_lo + (y-y_hi)
            if z < z_lo:
                z = z_hi + (z-z_lo)
            elif z > z_hi:
                z = z_lo + (z-z_hi)
            position = [x, y, z]
            ts.positions[i] = position

    with DCDWriter(outname+".dcd", len(u.atoms)) as W:
        for ts in new_traj:
            W.write(ts)#}}}
#########################################################


######################## catDCD #########################
def catDCD(top_file, basename, numFiles, outname):#{{{
    import MDAnalysis as mda 
    from MDAnalysis.coordinates.LAMMPS import DCDWriter

    traj_files = []

    for i in range(numFiles):
        traj_files.append(basename + str(i+1) + ".dcd")

    u = mda.Universe(top_file, traj_files)

    with DCDWriter(outname, len(u.atoms)) as W:
        for ts in u.trajectory:
            W.write(ts)#}}}
#########################################################


######################## catLog #########################
def catLog(basename, numFiles, outname, startLineFirst, endLineFirst, startLine, endLine):#{{{
    log_files = []
    lines = []

    for i in range(numFiles):
        log_files.append("log." + basename + str(i+1) + ".lammps")

    with open(log_files[0]) as file:
        for i in range(startLineFirst-1):
            file.readline()
        for i in range(endLineFirst-startLineFirst + 1): 
            line = file.readline()
            lines.append(line.lstrip())
            if i == 2:
                stepsBetween = int(line[:8])
            if i == endLineFirst-startLineFirst:
                lastNum = int(line[:8])

    for filename in log_files[1:]:
        with open(filename) as file:
            for i in range(startLine-1):
                file.readline()
            for i in range(endLine-startLine + 1): 
                line = file.readline()
                lastNum = lastNum + stepsBetween
                line = str(lastNum) + line[8:]
                lines.append(line)

    with open(outname, 'w') as file:
        for line in lines:
            file.write(line)#}}}
#########################################################


####################### excelPrep #######################
def excelPrep(startLine, endLine, selectionCriteria):#{{{
    filenames = gb.glob(selectionCriteria)

    for filename in filenames:
        with open(filename, 'r') as file:
            lines = file.readlines()
        newFilename = filename[:-7] + ".txt"
        with open(newFilename, 'w') as newFile:
            for line in lines[startLine-1:endLine]:
                newFile.write(line.lstrip())#}}}
#########################################################
