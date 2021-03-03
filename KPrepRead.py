#
# KPrepRead.py
# --------------------
# Authors: Addison Smith and Derek Bush
# Contributors: Anthony Gillepsie
# Written to be run on python 2 or 3

#########################################################
#--------------- Index of Function ----------------------
# read_pdb            Reads .pdb files directly from the protein data bank website
# readDiSulfides_pdb  Reads .pdb files directly from the protein data bank website for disulfide information
# readGO_pdb          Reads .pdb files produced by the mmtsb web server
# read_rtf            Reads .rtf files to find an Atom's atomType (can be accomplished with readpsf)
# readpsfprm          Reads the .psf file and relevent .prm files to define system parameters
# readparam_GO        Reads .param files produced by the mmtsb web server
# editGO_DiSulfides   Takes a list of residues that participate in disulfide bonds and adds it to the list of bonds
# readLAMMPS_data     Reads .data files used in LAMMPS
# readLAMMPS_nc       Reads .mod files that contain Native Contact (nc) (ie. NBFIX) information
# readLAMMPS_aa       Reads .mod files that contain Amino Acid (aa) information
# readKPREP_kp        Reads .kp files that contain KPI data-dump of variables
#--------------------------------------------------------
# Note: All function in KPrepRead.py create new class variables.
#       To modify an exitsing class variable see KPrepFunctions.py
#--------------------------------------------------------
# Note: to toggle VIM folds: put cursor at your location -> command mode -> "za"
#       to create a fold: command mode -> ex. ":12,34fold" the numbers are line numbers
#--------------------------------------------------------
#########################################################


#########################################################
import KPrepInfo as KPI
import re                   # re is short for regular expression
import copy
import numpy as np
#########################################################


########### Class Variables and Functions ###############
#-------- Universal Variables and Functions -------------#{{{
#------- specifically written for internal use ----------
floatPat = '[0-9.-]+\s*'    # "float pattern"   Looks for any float number
                            # a re that defines a any (digit, decimal or negative) + any space
intPat = '\d+\s*'           # "integer pattern" Loods for any integer number
                            # a re that defines a any positive number + any space
resGoPat = 'G\d+\s*'        # "residue Go-model pattern" Looks for residues; ex. G123
                            # a re that defines a "G" + digit(s) + any space
capPat = '[^A-Z]*[A-Z]{3}[^A-Z]*'  # "capital pattern" 
                            # Looks for 3+ consecutive upper-case letters ex. BOND
spaceFloatPat = '\s+[0-9.-]+\s*'    # "space + float + space pattern"   Looks for any float number
#--------------------------------------------------------#}}}
#########################################################


##################### read_pdb ##########################
def read_pdb(fileName, HETATM=True):#{{{
    # PDB file format explained on the RCSB PDF literature: PDB File Format v. 3.3 page: 187
    # Inputs:
    # FileName      This is a PDB file
    # HETATM        Binary True/False. If True, will read HETATMs from the PDB. Some models or simulators don't 
    #                   need HETATMs and so define as False.

    # Actions:
    # Will read the pdb and extract the following information for each atom (see KPI.Atom class):
    #      type, coordinates, atomSeq, atomName, atomNameAlt, resSeq, resName, chainID, q
    # Will sequentially define the follow variables based on sequencial order:
    #      atomNum, resNum
    # Will extract information from ATOM and HETATM types

    # Outputs:
    # An array of "KPI.Atom" class variables that contains all reported atoms (exclude ANISOU atoms)

    # Set-up and Pre-define variables
    arrAtom = []      # Array of all the HETATM and ATOM in the pdb
    i = 0               # Atom counter for each ATOM and HETATM encountered
    j = 0               # Residue counter for each ATOM and HETATM encountered
    additionAdjust = 0  # addition (insertion) offset if established PDB sequence was added to

    # Begin reading the pdb
    print("Reading file: %s" % (fileName))
    PDBfile = open(fileName, 'r')
    for line in PDBfile :
        comment = line[0:6]

        if HETATM == True and (comment == "ATOM  " or comment == "HETATM") :
            # Read line and extract information
            arrAtom.append(KPI.Atom())                      # Spacing
            arrAtom[i].type        = comment                # 6
            arrAtom[i].atomSeq     = int(line[6:11])        # 5 (single blank space after)
            arrAtom[i].atomName    = line[12:16]            # 4
            arrAtom[i].atomNameAlt = line[16]               # 1
            arrAtom[i].resType     = line[17:20]            # 3 (single blank space after)
            arrAtom[i].chainType   = line[21]               # 1
            arrAtom[i].resSeq      = int(line[22:26])       # 4 (4 blank spaces after)
            addition               = line[26]               # 1 A-Z indicating an addition to the PDB standard (rare occurance)
            x = float(line[30:38])                          # 8
            y = float(line[38:46])                          # 8
            z = float(line[46:54])                          # 8
            arrAtom[i].coordinates = [x, y, z]
            if line[54:60] == '      ' : arrAtom[i].o = ''  # 6
            else :                       arrAtom[i].o = float(line[54:60])
            if line[60:66] == '      ' : arrAtom[i].tf = '' # 6 (11 blank spaces after)
            else :                       arrAtom[i].tf = float(line[60:66])
            arrAtom[i].element     = line[76:78].rstrip()   # 2
            arrAtom[i].q           = line[78:80].rstrip()   # 2

            # Sequential Numbering as seen by the read
            arrAtom[i].atomNum = i+1

            if j == 0 : 
                arrAtom[i].resNum = j+1
                j += 1
            elif addition in KPI.alphabet and arrAtom[i-1].resSeq == arrAtom[i].resSeq and additionCheck == False : 
                arrAtom[i].resNum = j+1
                additionCheck = True
                j += 1
            elif arrAtom[i-1].resSeq == arrAtom[i].resSeq : arrAtom[i].resNum = j
            else :
                arrAtom[i].resNum = j+1
                additionCheck = False
                j += 1
            
            # Increment Atom counter
            i += 1


        elif HETATM == False and comment == "ATOM  " :
            # Read line and extract information
            arrAtom.append(KPI.Atom())                      # Spacing
            arrAtom[i].type        = comment                # 6
            arrAtom[i].atomSeq     = int(line[6:11])        # 5 (single blank space after)
            arrAtom[i].atomName    = line[12:16]            # 4
            arrAtom[i].atomNameAlt = line[16]               # 1
            arrAtom[i].resType     = line[17:20]            # 3 (single blank space after)
            arrAtom[i].chainType   = line[21]               # 1
            arrAtom[i].resSeq      = int(line[22:26])       # 4 (4 blank spaces after)
            addition               = line[26]               # 1 A-Z indicating an addition to the PDB standard (rare occurance)
            x = float(line[30:38])                          # 8
            y = float(line[38:46])                          # 8
            z = float(line[46:54])                          # 8
            arrAtom[i].coordinates = [x, y, z]
            if line[54:60] == '      ' : arrAtom[i].o = ''  # 6
            else :                       arrAtom[i].o = float(line[54:60])
            if line[60:66] == '      ' : arrAtom[i].tf = '' # 6 (11 blank spaces after)
            else :                       arrAtom[i].tf = float(line[60:66])
            arrAtom[i].element     = line[76:78].rstrip()   # 2
            arrAtom[i].q           = line[78:80].rstrip()   # 2

            # Sequential Numbering as seen by the read
            arrAtom[i].atomNum = i+1

            if j == 0 : 
                arrAtom[i].resNum = j+1
                j += 1
            elif addition in KPI.alphabet and arrAtom[i-1].resSeq == arrAtom[i].resSeq and additionCheck == False : 
                arrAtom[i].resNum = j+1
                additionCheck = True
                j += 1
            elif arrAtom[i-1].resSeq == arrAtom[i].resSeq : arrAtom[i].resNum = j
            else :
                arrAtom[i].resNum = j+1
                additionCheck = False
                j += 1

            # Increment Atom counter
            i += 1
    # End For

    # Define resType using KPI.dictAA dictionary
    dictUAA = {}
    uAAresID = 30 # reserving 1-29 for natural and published AA values
    for atom in arrAtom :
        dictEntry = atom.resType
        if dictEntry in KPI.dictAA : atom.resID = KPI.dictAA[dictEntry]
        if dictEntry in dictUAA : atom.resID = dictUAA[dictEntry]
        if dictEntry not in KPI.dictAA and dictEntry not in dictUAA : 
            print("\nWARNING: %s is not in the KPI dictionary!!" % (dictEntry))
            print("Setting %s to the unnatural AA resID: %i\n" % (dictEntry, uAAresID))
            newEntry = {dictEntry : uAAresID}
            dictUAA.update(newEntry)
            uAAresID += 1

    # This function cannont and does not define:
    #atomType
    #atomID
    #mw

    PDBfile.close()
    
    return arrAtom#}}}
#########################################################


################# readDiSulfides_pdb ####################
def readDiSulfides_pdb(fileName) :#{{{
    # PDB file format explained on the RCSB PDF literature: PDB File Format v. 3.3 page: 187
    # Inputs:
    # FileName      This is a PDB file

    # Actions:
    # Will read the pdb and extract a list of all residues that participate in disulfide bonds.

    # Outputs:
    # A list of all residues that participate in disulfide bonds. This list is meant to be used
    #   by other process or functions. Doesn't output anything to a file.

    # Set-up and Pre-define variables
    arrDiS = []         # Array of all the disulfide residue pairs
    i = 0               # Number of disulfide bonds in protein
    
    # Begin reading the pdb
    print("Reading file %s for disulfide analysis" % (fileName))
    PDBfile = open(fileName, 'r')
    for line in PDBfile :
        comment = line[0:6]
        if comment == "SSBOND" :
            i += 1
            line = line.split()
            arrDiS.append([int(line[4]), int(line[7])])
    PDBfile.close()

    # Put variable to be returned in martix format
    arrDiS = np.matrix(arrDiS)

    if i > 0 :
        print("Found %i disulfide bond(s) in %s" % (i, fileName))
    else :
        print("Did not find any disulfide bonds in %s" % (fileName))
        return []

    return arrDiS#}}}
#########################################################


##################### readGO_pdb ########################
def readGO_pdb (fileName) :#{{{
    # NOTE: the function readpdb can also read GO_model pdb files!
    # Inputs:
    # FileName   This must be a Go model style PDB file

    # Actions:
    # Will read the pdb created by mmtsb web server and extract the following information for each atom (see KPI.Atom class):
    #      type, coordinates, atomID, atomNum, atomName, resID, resNum, resName
    # At time of writing, mmtsb automatically defines resName as ALA.
    # At time of writing, mmtsb assumes atomID and resID both start at 1 therefor, atomID=atomNum & resID=resNum
    # At time of writing, mmtsb automatically defines atomName as CA
    # At time of writing, mmtsb does not provide information on multiple chains
    # Will extract information from ATOM and HETATM types

    # Outputs:
    # An array of "KPI.Atom" class variables that contains all ATOM 

    # Set-up and Pre-define variables
    arrAtom = []      # Array of all the ATOM in the pdb
    i = 0               # Atom counter for each ATOM encountered

    # Begin reading the PDB
    PDBfile = open(fileName, 'r')
    for line in PDBfile :
        comment = line[0:6]

        if comment == "ATOM  " :
            # Read line and extract information
            arrAtom.append(KPI.Atom())
            arrAtom[i].type        = comment
            arrAtom[i].atomSeq     = int(line[6:11])
            arrAtom[i].atomNum     = int(line[6:11])
            arrAtom[i].atomID      = int(line[6:11])
            arrAtom[i].atomName    = line[12:16]
            arrAtom[i].resType     = line[17:20]
            arrAtom[i].resSeq      = int(line[22:26])
            arrAtom[i].resNum      = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            arrAtom[i].coordinates = [x, y, z]
        
            # Increment Atom counter
            i += 1
    # End For
    PDBfile.close()

    return arrAtom#}}}
#########################################################


###################### read_rtf ##########################
def read_rtf(arrFileName, arrAtom) :#{{{
    # Note: This read ONLY adds the additional atomType field to your atomArray. You must have an appropriate atomArray
    #       prior calling this function
    # Inputs:
    # arrFileName  This is an array of .rtf files the function will look in. If you are only looking in only one file, the
    #              entry MUST still be in the form of an array. For example: [theFile.rtf]. Will not work if you pass-in a string.
    # atomArray    An array of Atoms that are missing atomType fields

    # Actions:
    # a standard NTER and CTER will ALWAYS be assumed (to be edited)

    # Outputs:
    # An array of "KPI.Atoms" class variables with their associated atomType

    # Set-up and Pre-define variables
    totAtom = len(arrAtom)
    firstLine = False

    for i in range(totAtom) :
        resType = arrAtom[i].resType
        atomName = arrAtom[i].atomName + arrAtom[i].atomNameAlt
        for file in arrFileName :
            readFile = open(file, 'r')
            lastLine = 'END'
            
            for line in readFile:
                # find resType; make sure another residue doesn't come-up, omit termini from check
                if re.search('RESI ' + resType, line) : 
                    firstLine = True
                    continue
                if arrAtom[i].resSeq != arrAtom[0].resSeq and arrAtom[i].resSeq != arrAtom[-1].resSeq and re.search('RESI', line) :
                    firstLine = False

                # if residue has been found, search for the next atomName
                if firstLine == True and re.search(atomName, line[0:13]) : 
                    arrAtom[i].atomType = line.split()[2]
                    break

                # check to see if an atomName wasn't found
                if re.search(lastLine, line) and file == arrFileName[-1] :
                    print('CATASTROPHIC ERROR: Atom', atomRes, atomName, 'was not found! Quitting Routine!!!')
                    quit()
            
            # reset for next residue
            if firstLine == True :
                firstLine = False
                break

            # find line within that region with unique atomName
            # parce line and define atomType
            # break
            readFile.close()
    return arrAtom#}}}
#########################################################


##################### readpsfprm ########################
#Functions #{{{
def modifyATOM(arrline, i, arrAtom, dictAtom) :#{{{
    # Define initial atom values
    atomNumLine    = int(arrline[0])
    if atomNumLine != arrAtom[i].atomNum :
        print("CATASTROPHIC ERROR: atomNum in psf does not equal atomNum in array")
        print("Error occured on line %i of the NATOM section in psf" % (i))
        quit()
    arrAtom[i].atomType = arrline[5]
    arrAtom[i].q        = float(arrline[6])
    arrAtom[i].mw       = float(arrline[7])

    # Identify atomID
    dictEntry = "%s" % (arrAtom[i].atomType)
    if dictEntry in dictAtom : arrAtom[i].atomID = int(dictAtom[dictEntry])
    if dictEntry not in dictAtom :
        arrAtom[i].atomID = len(dictAtom)+1
        newEntry = {dictEntry : str(arrAtom[i].atomID)}
        dictAtom.update(newEntry)

    return arrAtom, dictAtom#}}}

def createBOND(arrPair, i, arrAtom, dictBond, arrFileNamePRM) :#{{{
    xbond = KPI.Bond()

    # Define initial bond values
    xbond.bondNum = i+1
    xbond.bond[0] = arrPair[0]
    xbond.bond[1] = arrPair[1]
    xbond.bondType[0] = str(arrAtom[arrPair[0]-1].atomType)
    xbond.bondType[1] = str(arrAtom[arrPair[1]-1].atomType)

    # Identify bondID
    dictEntry0 = "%s %s" % (xbond.bondType[0], xbond.bondType[1])
    dictEntry1 = "%s %s" % (xbond.bondType[1], xbond.bondType[0])
    if dictEntry0 in dictBond : xbond.bondID = int(dictBond[dictEntry0])
    if dictEntry1 in dictBond : xbond.bondID = int(dictBond[dictEntry1])
    if dictEntry0 not in dictBond and dictEntry1 not in dictBond :
        xbond.bondID = len(dictBond)+1
        newEntry = {dictEntry0 : str(xbond.bondID)}
        dictBond.update(newEntry)

    # Find epsilon and sigma
    fileTot = len(arrFileNamePRM)
    bondPattern0 = '%s\s+%s\s+[0-9.-]+\s+[0-9.-]+' % (xbond.bondType[0], xbond.bondType[1])
    bondPattern1 = '%s\s+%s\s+[0-9.-]+\s+[0-9.-]+' % (xbond.bondType[1], xbond.bondType[0])
    for iFile in arrFileNamePRM :
        readFile = open(iFile, 'r')
        data = readFile.readlines()
        readFile.close()
        for line in data :
            if re.match(bondPattern0, line) :
                xbond.epsil = float(line.split()[2])
                xbond.sigma = float(line.split()[3])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return xbond, dictBond
            if re.match(bondPattern1, line) :
                xbond.epsil = float(line.split()[2])
                xbond.sigma = float(line.split()[3])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return xbond, dictBond
            # Check to see if you never defined epsil and sigma
            if line == data[-1] and iFile == arrFileNamePRM[-1] :
                print("CATASTROPHIC ERROR: could not find the LJ parameters for : %s" % (dictEntry0))
                print("Setting epsilon and sigma to 0 and 0, respectively.")
                xbond.epsil = 0.0
                xbond.sigma = 0.0
                return xbond, dictBond

    print("Error in creatBOND: should not have ever gotten here!!!")
    quit()#}}}

def createANGLE(arrPair, i, arrAtom, dictAngle, arrFileNamePRM) :#{{{
    xangle = KPI.Angle()

    # Define initial angle values
    xangle.angleNum = i+1
    xangle.angle[0] = arrPair[0]
    xangle.angle[1] = arrPair[1]
    xangle.angle[2] = arrPair[2]
    xangle.angleType[0] = str(arrAtom[arrPair[0]-1].atomType)
    xangle.angleType[1] = str(arrAtom[arrPair[1]-1].atomType)
    xangle.angleType[2] = str(arrAtom[arrPair[2]-1].atomType)

    # Identify angleID
    dictEntry0 = "%s %s %s" % (xangle.angleType[0], xangle.angleType[1], xangle.angleType[2])
    dictEntry1 = "%s %s %s" % (xangle.angleType[2], xangle.angleType[1], xangle.angleType[0])
    if dictEntry0 in dictAngle : xangle.angleID = int(dictAngle[dictEntry0])
    if dictEntry1 in dictAngle : xangle.angleID = int(dictAngle[dictEntry1])
    if dictEntry0 not in dictAngle and dictEntry1 not in dictAngle :
        xangle.angleID = len(dictAngle)+1
        newEntry = {dictEntry0 : str(xangle.angleID)}
        dictAngle.update(newEntry)

    # Find epsilon and sigma
    anglePattern0 = '%s\s+%s\s+%s\s+[0-9.-]+\s+[0-9.-]+' % \
                     (xangle.angleType[0], xangle.angleType[1], xangle.angleType[2])
    anglePattern1 = '%s\s+%s\s+%s\s+[0-9.-]+\s+[0-9.-]+' % \
                     (xangle.angleType[2], xangle.angleType[1], xangle.angleType[0])
    for iFile in arrFileNamePRM :
        readFile = open(iFile, 'r')
        data = readFile.readlines()
        readFile.close()
        for line in data :
            if re.match(anglePattern0, line) :
                xangle.epsil = float(line.split()[3])
                xangle.sigma = float(line.split()[4])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return xangle, dictAngle
            if re.match(anglePattern1, line) :
                xangle.epsil = float(line.split()[3])
                xangle.sigma = float(line.split()[4])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return xangle, dictAngle
            # Check to see if you never defined epsil and sigma
            if line == data[-1] and iFile == arrFileNamePRM[-1] :
                print("CATASTROPHIC ERROR: could not find the LJ parameters for : %s" % (dictEntry0))
                print("Setting epsilon and sigma to 0 and 0, respectively.")
                xangle.epsil = 0.0
                xangle.sigma = 0.0
                return xangle, dictAngle

    print("Error in creatAngle: should not have ever gotten here!!!")
    quit()#}}}

def createDIHEDRAL(arrPair, i, arrAtom, dictDihedral, arrFileNamePRM) :#{{{
    xarray = []
    xdihedral = KPI.Dihedral()

    # Define initial dihedral values
    saveLine = [] # Variable holds exact matches
    wildLine = [] # Variable holds all wildcard matches
    xdihedral.dihedralNum = i+1
    xdihedral.dihedral[0] = arrPair[0]
    xdihedral.dihedral[1] = arrPair[1]
    xdihedral.dihedral[2] = arrPair[2]
    xdihedral.dihedral[3] = arrPair[3]
    xdihedral.dihedralType[0] = str(arrAtom[arrPair[0]-1].atomType)
    xdihedral.dihedralType[1] = str(arrAtom[arrPair[1]-1].atomType)
    xdihedral.dihedralType[2] = str(arrAtom[arrPair[2]-1].atomType)
    xdihedral.dihedralType[3] = str(arrAtom[arrPair[3]-1].atomType)

    # Identify dihedralID
    def dihedralID(xdihedral, dictDihedral) :
        dictEntry0 = "%s %s %s %s %s" % (xdihedral.dihedralType[0], xdihedral.dihedralType[1], \
                                      xdihedral.dihedralType[2], xdihedral.dihedralType[3], xdihedral.n)
        dictEntry1 = "%s %s %s %s %s" % (xdihedral.dihedralType[3], xdihedral.dihedralType[2], \
                                      xdihedral.dihedralType[1], xdihedral.dihedralType[0], xdihedral.n)
        if dictEntry0 in dictDihedral : dihID = int(dictDihedral[dictEntry0])
        if dictEntry1 in dictDihedral : dihID = int(dictDihedral[dictEntry1])
        if dictEntry0 not in dictDihedral and dictEntry1 not in dictDihedral :
            dihID = len(dictDihedral)+1
            newEntry = {dictEntry0 : str(dihID)}
            dictDihedral.update(newEntry)
        return dihID, dictDihedral

    # Exact matches
    dihedralPattern0 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (xdihedral.dihedralType[0], xdihedral.dihedralType[1], \
                      xdihedral.dihedralType[2], xdihedral.dihedralType[3])
    dihedralPattern1 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (xdihedral.dihedralType[3], xdihedral.dihedralType[2], \
                      xdihedral.dihedralType[1], xdihedral.dihedralType[0])
    # Wildcard X matches
    dihedralPattern00 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     ("X", xdihedral.dihedralType[1], \
                      xdihedral.dihedralType[2], "X")
    dihedralPattern11 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     ("X", xdihedral.dihedralType[2], \
                      xdihedral.dihedralType[1], "X")

    # Search files and find matches
    for iFile in arrFileNamePRM :
        readFile = open(iFile, 'r')
        data = readFile.readlines()
        readFile.close()
        for line in data :
            if re.match(dihedralPattern0, line) or re.match(dihedralPattern1,line) :
                saveLine.append(line.rstrip())
                #print(line.rstrip()) # Print out entire line for troubleshooting

            # Wildcard dihedral matches
            elif re.match(dihedralPattern00, line) or re.match(dihedralPattern11,line) :
                wildLine.append(line.rstrip())
                #print(line.rstrip()) # Print out entire line for troubleshooting
            
            # Check to see if you never defined epsil and sigma
            if line == data[-1] and iFile == arrFileNamePRM[-1] :
                if wildLine == '' and saveLine == '' :
                    print("\nCAUTION: could not find an exact match for improper pattern : %s" % (dictEntry0))
                    print("reference dihedralNum: ", xdihedral.dihedralNum)
                    print("CATASTROPHIC ERROR: Could not find any related matches either...")
                    print("Setting epsilon and sigma to 0 and 0, respectively.\n")
                    xdihedral.K = 0.0
                    xdihedral.n = 0
                    xdihedral.d = 0.0
                    xarray.append(copy.copy(xdihedral))
                    return xarray, dictDihedral
    
    # If alternate options exist, define parameters accordingly
    if saveLine != [] :
        for line in saveLine :
            xdihedral.K = float(line.split()[4])
            xdihedral.n = int(line.split()[5])
            xdihedral.d = float(line.split()[6])
            dihID, dictDihedral = dihedralID(xdihedral, dictDihedral)
            xdihedral.dihedralID = dihID
            xarray.append(copy.copy(xdihedral))
        return xarray, dictDihedral
    if wildLine != [] :
        for line in wildLine :
            xdihedral.K = float(line.split()[4])
            xdihedral.n = int(line.split()[5])
            xdihedral.d = float(line.split()[6])
            dihID, dictDihedral = dihedralID(xdihedral, dictDihedral)
            xdihedral.dihedralID = dihID
            xarray.append(copy.copy(xdihedral))
        return xarray, dictDihedral

    print("Error in creatDihedral: should not have ever gotten here!!!")
    quit()#}}}

def createIMPROPER(arrPair, i, arrAtom, dictImproper, arrFileNamePRM) :#{{{
    ximproper = KPI.Improper()

    # Define initial improper values
    wildLine = '' # Variable holds all wildcard matches
    saveLine = '' # Variable holds similar-type matches
    ximproper.improperNum = i+1
    ximproper.improper[0] = arrPair[0]
    ximproper.improper[1] = arrPair[1]
    ximproper.improper[2] = arrPair[2]
    ximproper.improper[3] = arrPair[3]
    ximproper.improperType[0] = str(arrAtom[arrPair[0]-1].atomType)
    ximproper.improperType[1] = str(arrAtom[arrPair[1]-1].atomType)
    ximproper.improperType[2] = str(arrAtom[arrPair[2]-1].atomType)
    ximproper.improperType[3] = str(arrAtom[arrPair[3]-1].atomType)

    # Identify improperID
    dictEntry0 = "%s %s %s %s" % (ximproper.improperType[0], ximproper.improperType[1], \
                                  ximproper.improperType[2], ximproper.improperType[3])
    dictEntry1 = "%s %s %s %s" % (ximproper.improperType[3], ximproper.improperType[2], \
                                  ximproper.improperType[1], ximproper.improperType[0])
    if dictEntry0 in dictImproper : ximproper.improperID = int(dictImproper[dictEntry0])
    if dictEntry1 in dictImproper : ximproper.improperID = int(dictImproper[dictEntry1])
    if dictEntry0 not in dictImproper and dictEntry1 not in dictImproper :
        ximproper.improperID = len(dictImproper)+1
        newEntry = {dictEntry0 : str(ximproper.improperID)}
        dictImproper.update(newEntry)

    # Exact matches
    improperPattern0 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[0], ximproper.improperType[1], \
                      ximproper.improperType[2], ximproper.improperType[3])
    improperPattern1 = '%s\s+%s\s+%s\s+%s\s+[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[3], ximproper.improperType[2], \
                      ximproper.improperType[1], ximproper.improperType[0])
    # Wildcard X matches
    improperPattern00 = '.*%s.*%s.*%s.*%s.*[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[0], "X", \
                      "X", ximproper.improperType[3])
    improperPattern11 = '.*%s.*%s.*%s.*%s.*[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[3], "X", \
                      "X", ximproper.improperType[0])
    # Similar-type matches
    improperPattern000 = '.*%s.*%s.*%s.*%s.*[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[0], ximproper.improperType[1], \
                      ximproper.improperType[2], ximproper.improperType[3])
    improperPattern111 = '.*%s.*%s.*%s.*%s.*[0-9.-]+\s+\d\s+[0-9.-]+' % \
                     (ximproper.improperType[3], ximproper.improperType[2], \
                      ximproper.improperType[1], ximproper.improperType[0])
    # Search files and find matches
    for iFile in arrFileNamePRM :
        readFile = open(iFile, 'r')
        data = readFile.readlines()
        readFile.close()
        for line in data :
            if re.match(improperPattern0, line) or re.match(improperPattern1,line) :
                ximproper.K = float(line.split()[4])
                ximproper.chi = float(line.split()[6])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return ximproper, dictImproper

            # Wildcard improper matches
            elif re.match(improperPattern00, line) or re.match(improperPattern11,line) :
                wildLine += line
                #print(line.rstrip()) # Print out entire line for troubleshooting
            
            # Similar-type matches
            elif re.match(improperPattern000,line) or re.match(improperPattern111,line) : 
                saveLine += line
                #print(line.rstrip()) # Print out entire line for troubleshooting

            # Check to see if you never defined epsil and sigma
            if line == data[-1] and iFile == arrFileNamePRM[-1] :
                if wildLine == '' and saveLine != '' :
                    print("\nCAUTION: could not find an exact match for improper pattern : %s" % (dictEntry0))
                    print("reference improperNum: ", ximproper.improperNum)
                if wildLine == '' and saveLine == '' :
                    print("\nCAUTION: could not find an exact match for improper pattern : %s" % (dictEntry0))
                    print("reference improperNum: ", ximproper.improperNum)
                    print("CATASTROPHIC ERROR: Could not find any related matches either...")
                    print("Setting epsilon and sigma to 0 and 0, respectively.\n")
                    ximproper.epsil = 0.0
                    ximproper.sigma = 0.0
                    return ximproper, dictImproper
    
    # If alternate options exist, define parameters accordingly
    if wildLine != '' :
        ximproper.K   = float(wildLine.split()[4])
        ximproper.chi = float(wildLine.split()[6])
        return ximproper, dictImproper
    if saveLine != '' :
        print("WARNING: related match(es) found, saving LJ parameters of first encounter:\n", saveLine.rstrip())
        ximproper.K   = float(saveLine.split()[4])
        ximproper.chi = float(saveLine.split()[6])
        return ximproper, dictImproper

    print("Error in creatImproper: should not have ever gotten here!!!")
    quit()#}}}

def createNONBONDED(xatom, arrFileNamePRM) :#{{{
    xnonbonded = KPI.Nonbonded()

    # Define initial bond values
    xnonbonded.atomNum = xatom.atomNum
    xnonbonded.atomID = xatom.atomID
    xnonbonded.atomType = xatom.atomType

    # Find Rmin/2 and epsilon
    fileTot = len(arrFileNamePRM)
    nonbondedPattern = '%s\s+[0-9.-]+\s+[0-9.-]+\s+[0-9.-]+' % (xnonbonded.atomType)
    for iFile in arrFileNamePRM :
        readFile = open(iFile, 'r')
        data = readFile.readlines()
        readFile.close()
        for line in data :
            if re.match(nonbondedPattern, line) :
                xnonbonded.epsil = float(line.split()[2])
                xnonbonded.rmin2 = float(line.split()[3])
                #print(line.rstrip()) # Print out entire line for troubleshooting
                return xnonbonded
            # Check to see if you never defined epsil and rmin2
            if line == data[-1] and iFile == arrFileNamePRM[-1] :
                print("\nCATASTROPHIC ERROR: could not find the LJ parameters for : %s" % \
                        (xnonbonded.atomType))
                print("Setting Rmin/2 and epsilon to 0 and 0, respectively.\n")
                xnonbonded.epsil = 0.0
                xnonbonded.rmin2 = 0.0
                return xnonbonded

    print("Error in creatNONBONDED: should not have ever gotten here!!!")
    quit()#}}}
#}}}

def readpsfprm(fileNamePSF, arrFileNamePRM, arrAtom) :#{{{
    
    
    # Actions:

    # Outputs:

    # Set-up and Pre-define variables
    orig_arrAtom = arrAtom.copy()
    totAtom = len(arrAtom)
    dictAtom = {}
    arrBond = []
    dictBond = {}
    arrAngle = []
    dictAngle = {}
    arrDihedral = []
    dictDihedral = {}
    arrImproper = []
    dictImproper = {}
    arrNonbonded = []

    # Read & parce data from .psf
    readFile = open(fileNamePSF, 'r')#{{{
    data = readFile.readlines()
    readFile.close()

    # Parce data into each sections
    for i, line in enumerate(data) :
        if re.search('NATOM', line) :
            nAtom = int(line.split()[0])
            atomdata = data[i+1 : i+nAtom+1]
            #print(atomdata[-1])
        if re.search('NBOND', line) :
            nBond = int(line.split()[0])
            rows = int(nBond/4)   # Assume there are 4 bond pairs in each line
            if nBond%4 != 0 : rows += 1
            bonddata = data[i+1 : i+rows+1]
            #print(bonddata[-1])
        if re.search('NTHETA', line) :
            nTheta = int(line.split()[0])
            rows = int(nTheta/3) # Assume there are 3 angle pairs in each line
            if nTheta%3 != 0 : rows += 1
            angledata = data[i+1 : i+rows+1]
            #print(angledata[-1])
        if re.search('NPHI', line) :
            nPhi = int(line.split()[0])
            rows = int(nPhi/2)   # Assume there are 2 dihedral pairs in each line
            if nPhi%2 != 0 : rows += 1
            dihedraldata = data[i+1 : i+rows+1]
            #print(dihedraldata[-1])
        if re.search('NIMPHI', line) :
            nImPhi = int(line.split()[0])
            rows = int(nImPhi/2) # Assume there are 2 improper pairs in each line
            if nImPhi%2 != 0 : rows += 1
            improperdata = data[i+1 : i+rows+1]
            #print(improperdata[-1])
        if re.search('NDON', line) :
            nDon = int(line.split()[0])
            rows = int(nDon/4)   # Assume there are 4 donor pairs in each line
            if nDon%4 != 0 : rows += 1
            donordata = data[i+1 : i+rows+1]
            #print(donordata[-1])
        if re.search('NACC', line) :
            nAcc = int(line.split()[0])
            rows = int(nAcc/4)   # Assume there are 4 acceptor pairs in each line
            if nAcc%4 != 0 : rows += 1
            acceptordata = data[i+1 : i+rows+1]
            #print(acceptordata[-1])
        # Skipping NNB section
        # Skipping NGRP section
        # Skipping MOLNT section
        # Skipping NUMLP section
        if re.search('NCRTERM', line) :
            nCrTerm = int(line.split()[0])
            crosstermdata = data[i+1 : i+nCrTerm+1]
            #print(crosstermdata[-1])
    data = 0 # to save on memory#}}}

    # Define Atom parameters: atomType, mw, q
    print("Defining atom type, MW and charge...")#{{{
    for i, line in enumerate(atomdata) :
        arrline = line.split()
        arrAtom, dictAtom = modifyATOM(arrline, i, arrAtom, dictAtom)#}}}

    # Define Bond parameters
    print("Defining Bond parameters...")#{{{
    for i, line in enumerate(bonddata) :
        line = line.split()
        pairTot = int(len(line)/2)
        for j in range(pairTot) :
            pair = [int(line[2*j]), int(line[2*j+1])]
            bond, dictBond = createBOND(pair, 4*i+j, arrAtom, dictBond, arrFileNamePRM)
            arrBond.append(bond) # For some reason bondType info doesn't transfer properly
            arrBond[-1].bondType = [bond.bondType[0], bond.bondType[1]]
    
    # For Troubleshooting
    #print(dictBond)
    #ii = 3
    #print(arrBond[ii].bond, arrBond[ii].bondType, arrBond[ii].bondID, arrBond[ii].bondNum, arrBond[ii].epsil)
    #ii = -1
    #print(arrBond[ii].bond, arrBond[ii].bondType, arrBond[ii].bondID, arrBond[ii].bondNum, arrBond[ii].epsil)#}}}

    # Define Angle parameters
    print("Defining Angle parameters...")#{{{
    for i, line in enumerate(angledata) :
        line = line.split()
        pairTot = int(len(line)/3)
        for j in range(pairTot) :
            pair = [int(line[3*j]), int(line[3*j+1]), int(line[3*j+2])]
            angle, dictAngle = createANGLE(pair, 3*i+j, arrAtom, dictAngle, arrFileNamePRM)
            arrAngle.append(angle) # For some reason angleType info doesn't transfer properly
            arrAngle[-1].angleType = [angle.angleType[0], angle.angleType[1], angle.angleType[2]]
    
    # For Troubleshooting
    #print(dictAngle)
    #ii = 3
    #print(arrAngle[ii].angle, arrAngle[ii].angleType, arrAngle[ii].angleID, arrAngle[ii].angleNum, arrAngle[ii].epsil)
    #ii = -1
    #print(arrAngle[ii].angle, arrAngle[ii].angleType, arrAngle[ii].angleID, arrAngle[ii].angleNum, arrAngle[ii].epsil)#}}}
    
    # Define Dihedral parameters
    print("Defining Dihedral parameters...")#{{{
    for i, line in enumerate(dihedraldata) :
        line = line.split()
        pairTot = int(len(line)/4)
        for j in range(pairTot) :
            pair = [int(line[4*j]), int(line[4*j+1]), int(line[4*j+2]), int(line[4*j+3])]
            arrpair, dictDihedral = createDIHEDRAL(pair, 2*i+j, arrAtom, dictDihedral, arrFileNamePRM)
            for dihedral in arrpair :
                arrDihedral.append(copy.copy(dihedral)) # For some reason dihedralType doesn't copy properly
                arrDihedral[-1].dihedralType = [dihedral.dihedralType[0], dihedral.dihedralType[1], \
                                                dihedral.dihedralType[2], dihedral.dihedralType[3]]
    
    # For Troubleshooting
    #print(dictDihedral)
    #ii = 3
    #print(arrDihedral[ii].dihedral, arrDihedral[ii].dihedralType, arrDihedral[ii].dihedralID, arrDihedral[ii].dihedralNum, arrDihedral[ii].K)
    #ii = -1
    #print(arrDihedral[ii].dihedral, arrDihedral[ii].dihedralType, arrDihedral[ii].dihedralID, arrDihedral[ii].dihedralNum, arrDihedral[ii].K)#}}}
 
    # Define Improper parameters
    print("Defining Improper parameters...")#{{{
    for i, line in enumerate(improperdata) :
        line = line.split()
        pairTot = int(len(line)/4)
        for j in range(pairTot) :
            pair = [int(line[4*j]), int(line[4*j+1]), int(line[4*j+2]), int(line[4*j+3])]
            improper, dictImproper = createIMPROPER(pair, 2*i+j, arrAtom, dictImproper, arrFileNamePRM)
            arrImproper.append(improper) # For some reason improperType info doesn't transfer properly
            arrImproper[-1].improperType = [improper.improperType[0], improper.improperType[1], \
                                               improper.improperType[2], improper.improperType[3]]
    # For Troubleshooting
    #print(dictImproper)
    #ii = 3
    #print(arrImproper[ii].improper, arrImproper[ii].improperType, arrImproper[ii].improperID, arrImproper[ii].improperNum, arrImproper[ii].chi)
    #ii = -1
    #print(arrImproper[ii].improper, arrImproper[ii].improperType, arrImproper[ii].improperID, arrImproper[ii].improperNum, arrImproper[ii].K)#}}}

    # Define Nonbonded parameters: atomType, mw, q
    print("Defining Nonbonded parameters...")#{{{
    for i, atom in enumerate(arrAtom) :
        nonbonded = createNONBONDED(atom, arrFileNamePRM)
        arrNonbonded.append(nonbonded)#}}}

    return arrBond, arrAngle, arrDihedral, arrImproper, arrNonbonded#}}}
#########################################################


#################### readGO_param #######################
#Functions #{{{
def createGoBOND(strLine, i) :#{{{
    # Inputs:
    # strLine   This is a string variable that contains pertinant BOND information from a param file in the format:
    #             resGoPat + resGoPat + floatPat + floatPat
    # i         This is the consecutive occurance of this parameter ie. the "Num"

    # Actions:
    # Will interpret the line, create a Bond class and appropriatley define class information:
    #      bondID, bondNum, bond, sigma, epsil
    # At time of writing, Go-model is written where ID parameters = Num parameters

    # Outputs:
    # A single Bond class variable

    xbond = KPI.Bond()

    # Split line on all "G" and spaces. [1:-1] removes the leading and trailing strings
    xbond.bondID  = i
    xbond.bondNum = i
    strArray = re.split('[G\s]+', strLine)[1:-1]
    site1 = int(strArray[0])
    site2 = int(strArray[1])
    xbond.bond = [site1, site2]
    xbond.epsil = float(strArray[2])
    xbond.sigma = float(strArray[3])
    return xbond#}}}

def createGoANGLE(strLine, i) :#{{{
    # Inputs:
    # strLine   This is a string variable that contains pertinant ANGLE information from a param file
    #             resGoPat + resGoPat + resGoPat + floatPat + floatPat
    # i         This is the consecutive occurance of this parameter ie. the "Num"

    # Actions:
    # Will interpret the line, create a Angle class and appropriatley define class information:
    #      angleID, angleNum, angle, sigma, epsil
    # At time of writing, Go-model is written where ID parameters = Num parameters

    # Outputs:
    # A single Angle class variable

    xangle = KPI.Angle()

    # Split line on all "G" and spaces. [1:-1] removes the leading and trailing strings
    xangle.angleID  = i
    xangle.angleNum = i
    strArray = re.split('[G\s]+', strLine)[1:-1]
    site1 = int(strArray[0])
    site2 = int(strArray[1])
    site3 = int(strArray[2])
    xangle.angle = [site1, site2, site3]
    xangle.epsil = float(strArray[3])
    xangle.sigma = float(strArray[4])
    return xangle#}}}

def createGoDIHEDRAL(strLine, i) :#{{{
    # Inputs:
    # strLine   This is a string variable that contains pertinant DIHEDRAL information from a param file
    #             resGoPat + resGoPat + resGoPat + resGoPat + floatPat + floatPat
    # i         This is the consecutive occurance of this parameter ie. the "Num"

    # Actions:
    # Will interpret the line, create a Dihedral class and appropriatley define class information:
    #      dihedral, K, n, d
    # At time of writing, Go-model is written where ID parameters = Num parameters

    # Outputs:
    # A single Dihedral class variable

    xdihedral = KPI.Dihedral()

    # Split line on all "G" and spaces. [1:-1] removes the leading and trailing strings
    xdihedral.dihedralID  = i
    xdihedral.dihedralNum = i
    strArray = re.split('[G\s]+', strLine)[1:-1]
    site1 = int(strArray[0])
    site2 = int(strArray[1])
    site3 = int(strArray[2])
    site4 = int(strArray[3])
    xdihedral.dihedral = [site1, site2, site3, site4]
    xdihedral.K = float(strArray[4])
    xdihedral.n = int(strArray[5])
    xdihedral.d = float(strArray[6])
    return xdihedral#}}}

def createGoNONBONDED(strLine) :#{{{
    # Inputs:
    # strLine   This is a string variable that contains pertinant NONBONDED information from a param file
    #             resGoPat + floatPat + floatPat

    # Actions:
    # Will interpret the line, create a Nonbonded class and appropriatley define class information:
    #      atomID, rmin2, epsil

    # Outputs:
    # A single Nonbonded class variable

    xnonbonded = KPI.Nonbonded()

    # Split line on all "G" and spaces. [1:-1] removes the leading and trailing strings
    strArray = re.split('[G\s]+', strLine)[1:-1]
    site1 = int(strArray[0])
    xnonbonded.atomID  = site1
    xnonbonded.atomNum = site1
    xepsil = float(strArray[1])
    if xepsil >= 0 : xnonbonded.epsil = xepsil
    else :           xnonbonded.epsil = -xepsil
    xnonbonded.rmin2 = float(strArray[2])
    return xnonbonded#}}}

def createGoNBFIX(strLine, i) :#{{{
    # Inputs:
    # strLine   This is a string variable that contains pertinant NBFIX information from a param file
    #             resGoPat + resGoPat + floatPat + floatPat
    # i         This is the consecutive occurance of this parameter

    # Actions:
    # Will interpret the line, create a NBfix class and appropriatley define class information:
    #      nbfixID, nbfix, epsil, sigma

    # Outputs:
    # A single NBfix class variable

    xnbfix = KPI.NBfix()

    # Split line on all "G" and spaces. [1:-1] removes the leading and trailing strings
    xnbfix.nbfixID = i
    strArray = re.split('[G\s]+', strLine)[1:-1]
    site1 = int(strArray[0])
    site2 = int(strArray[1])
    xnbfix.nbfix = [site1, site2]
    xepsil = float(strArray[2])
    if xepsil >= 0 : xnbfix.epsil = xepsil
    else :           xnbfix.epsil = -xepsil
    xnbfix.sigma = float(strArray[3])
    return xnbfix#}}}
#}}}

def readGO_param (fileName, BOND=True, ANGLE=True, DIHEDRAL=True, NONBONDED=True, NBFIX=True) :#{{{
    # Inputs:
    # FileileName_var   This must be a GO model .param file
    # BOND, ANGLE, DIHEDRAL, NONBONDED, NBFIX   These bool variables tell the script if you want to
    #                   actually search for these values. You can do a specific read if you turn off
    #                   unwanted parameters

    # Actions:
    # Will read the param created by mmtsb web server and extract the following information:
    #      Bond, Angle, Dihedral, Nonbonded, and NBfix. 
    # The function will then create arrays of class variables that correspond to the proper data type.

    # Outputs:
    # The function will output a series of arrays in the order:
    #      arrBond, arrAngle, arrDihedral, arrNonbonded, arrNBfix. 

    # Warnings:
    # If there is a break/newline in the middle of a field in your go_param file, the script may skip-over that
    #      line and may skip all following lines too
    # If there is any sort of deviation from the standard pattern in the line, the script may skip that line

    # Set-up and Pre-defined variables
    arrBond = []
    arrAngle = []
    arrDihedral = []
    arrNonbonded = []
    arrNBfix = []
    activeB  = False    # "active Bond field" These bool variables inform if the field is actively reading or not
    activeA  = False    # "active Angle field"
    activeD  = False    # "active Dihedral field"
    activeN  = False    # "active Nonbonded field"
    activeNf = False    # "active NBfix field"
    linePat = ""        # "line pattern" An additional check to ensure what we are reading is in the right form
    modLine = ""        # If a line isn't in the correct pattern for the function, create a string that is
    i = 0               # A counter the defines how many units have been added to a given field
    arrReturn = []

    # Begin reading the PARAM
    PARAMfile = open(fileName, 'r')
    for line in PARAMfile :
        
        # BOND field reads
        if BOND == True and re.search('BOND\s', line) :
            # Activate field and define line pattern struture
            activeB = True
            linePat = resGoPat + resGoPat + floatPat + floatPat
            continue
        if activeB == True and re.search(linePat, line):
            # Define a new BOND
            arrBond.append(createGoBOND(line, i+1))
            i += 1
            continue
        if activeB == True and re.search(capPat, line) :
            # Close field and reset read variables
            print("Read-in %i BOND parameters" % i)
            activeB = False
            linePat = ""
            i = 0

        # ANGLE field reads
        if ANGLE == True and re.search('ANGLE\s', line) :
            # Activate field and define line pattern struture
            activeA = True
            linePat = resGoPat + resGoPat + resGoPat + floatPat + floatPat
            continue
        if activeA == True and re.search(linePat, line):
            # Define a new ANGLE
            arrAngle.append(createGoANGLE(line, i+1))
            i += 1
            continue
        if activeA == True and re.search(capPat, line) :
            # Close field and reset read variables
            print("Read-in %i ANGLE parameters" % i)
            activeA = False
            linePat = ""
            i = 0

        # DIHEDRAL field reads
        if DIHEDRAL == True and re.search('DIHEDRAL\s', line) :
            # Activate field and define line pattern struture
            activeD = True
            linePat = resGoPat + resGoPat + resGoPat + resGoPat + floatPat + "\d\s*" + floatPat
            continue
        if activeD == True and re.search(linePat, line):
            # Define a new DIHEDRAL
            arrDihedral.append(createGoDIHEDRAL(line, i+1))
            i += 1
            continue
        if activeD == True and re.search(capPat, line) :
            # Close field and reset read variables
            print("Read-in %i DIHEDRAL parameters" % i)
            activeD = False
            linePat = ""
            i = 0

        # NONBONDED field reads
        if NONBONDED == True and re.search('NONBONDED\s', line) :
            # Activate field and define line pattern struture
            activeN = True
            linePat = resGoPat + floatPat + floatPat + floatPat
            continue
        if activeN == True and re.search(linePat, line):
            modLine = line.split()[0] + " " + line.split()[2] + " " + line.split()[3] + " "
            # Define a new NONBONDED
            arrNonbonded.append(createGoNONBONDED(modLine))
            i += 1
            continue
        if activeN == True and re.search(capPat, line) and i > 1 :
            # Close field and reset read variables
            print("Read-in %i NONBONDED parameters" % i)
            activeN = False
            linePat = ""
            modLine = ""
            i = 0

        # NBfix field reads
        if NBFIX == True and re.search('NBFIX\s', line) :
            # Activate field and define line pattern struture
            activeNf = True
            linePat = resGoPat + resGoPat + floatPat + floatPat
            continue
        if activeNf == True and re.search(linePat, line):
            # Define a new NBFIX
            arrNBfix.append(createGoNBFIX(line, i+1))
            i += 1
            continue
        if activeNf == True and re.search(capPat, line) :
            # Close field and reset read variables
            print("Read-in %i NBFIX parameters" % i)
            activeNf = False
            linePat = ""
            i = 0
    # End For
    PARAMfile.close()

    # Compile all the arrays of data
    if BOND == True : arrReturn.append(arrBond)
    if ANGLE == True : arrReturn.append(arrAngle)
    if DIHEDRAL == True : arrReturn.append(arrDihedral)
    if NONBONDED == True: arrReturn.append(arrNonbonded)
    if NBFIX == True : arrReturn.append(arrNBfix)

    return arrReturn#}}}
#########################################################


################# editGO_DiSulfides #####################
# Functions#{{{
def distanceTwoPoints(X, Y, frst, scnd) :
    dist = np.sqrt((X[0]-Y[0])**2 + (X[1]-Y[1])**2 + (X[2]-Y[2])**2)
    if dist > 7 :
        print("\nCAUTION: Distance between residues is large for a disulfide bond\n\
                Residues: %i and %i; Distance = %.4f" % (frst, scnd, dist))
    return dist


def disulfideChecks(i, atom) :
    # Begin checks
    boolFail = False
    if i != atom.atomNum : print("missmatch in reference")
    if atom.resType != "CYS" : boolFail = True
    # Print error message
    if boolFail == True :
        print("\nWARNING: Some data might be incorrect!\nResidue Number = %i; Residue type = %s" \
             % (atom.atomNum, atom.resType))
#}}}

def editGO_DiSulfides(arrAtom, arrBond, arrDiSulfides) :#{{{
    # Note: This is an edit that creates additional disulfide bondes for a Go model. It is assumed
    #  that the user has already determined the disulfide bonds (see readDiSulfides_pdb) and transformed the
    #  the list of participating amino acids to comply with the Go model pdb residue IDs

    # Inputs: 
    #arrAtom is the array of Go model amino acids and is used for reference only.
    #arrBond is the array of already-created bonds from the Go model (see readGO_params)
    #arrDiSulfides contains the information about the participating disulfide bonds in the protein.
    #  assumed is that these numbers have already been transformed into the appropriate Go model residue IDs
    
    # Actions:
    #First the distance between the amino acids is determined.
    #Lastly the new bond is added to the system using the creatGoBOND function from readGo_param

    # Outputs:
    #arrBondNew is the updated list of bonds in the system that contains the new disulfide bonds

    # Set-up and Pre-define variables
    DiS_bondEnergy = 378.00 
    arrBondNew = arrBond.copy()
    length = len(arrDiSulfides)
    numBonds = len(arrBond)
    
    if length == 0 : return arrBond

    # Determine residues, check information, calculate distance and create new bond
    for i in range(length) :
        frst = arrDiSulfides[i,0]
        scnd = arrDiSulfides[i,1]
        frstAA = arrAtom[frst-1]
        scndAA = arrAtom[scnd-1]
        disulfideChecks(frst, frstAA)
        disulfideChecks(scnd, scndAA)
        frstCoord = frstAA.coordinates
        scndCoord = scndAA.coordinates
        dist = distanceTwoPoints(frstCoord, scndCoord, frst, scnd)
        strBond = " %i  %i  %.3f  %.4f " %(frst, scnd, DiS_bondEnergy, dist)
        bondNum = numBonds + i + 1
        arrBondNew.append(createGoBOND(strBond, bondNum))

    return arrBondNew
#}}}
#########################################################


#################### readKPREP_kp #######################
#Functions #{{{
def createKPIAtom(line) :#{{{
    xatom = KPI.Atom()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = ''
    xatom.type = line[0]
    if line[1] == "N/A" : line[1] = 0
    xatom.atomNum = int(line[1])
    if line[2] == "N/A" : line[2] = 0
    xatom.atomSeq = int(line[2])
    if line[3] == "N/A" : line[3] = ''
    xatom.atomName = line[3]
    if line[4] == "N/A" : line[4] = ''
    xatom.atomNameAlt = line[4]
    if line[5] == "N/A" : line[5] = 0
    xatom.resID = int(line[5])
    if line[6] == "N/A" : line[6] = ''
    xatom.resType = line[6]
    if line[7] == "N/A" : line[7] = ''
    xatom.chainType = line[7]
    if line[8] == "N/A" : line[8] = 0
    xatom.resNum = int(line[8])
    if line[9] == "N/A" : line[9] = 0
    xatom.resSeq = int(line[9])
    if line[10] == "N/A" or line[11] == "N/A" or line[12] == "N/A" : 
        line[10] = 0.0
        line[11] = 0.0
        line[12] = 0.0
    xatom.coordinates = [float(line[10]), float(line[11]), float(line[12])]
    if line[13] == "N/A" : line[13] = 0.0
    xatom.o = float(line[13])
    if line[14] == "N/A" : line[14] = 0.0
    xatom.tf = float(line[14])
    if line[15] == "N/A" : line[15] = ''
    xatom.element = line[15]
    if line[16] == "N/A" : line[16] = 0.0
    xatom.q = float(line[16])
    if line[17] == "N/A" : line[17] = 0.0
    xatom.mw = float(line[17])
    if line[18] == "N/A" : line[18] = 0
    xatom.atomID = int(line[18])
    if line[19] == "N/A" : line[19] = ''
    xatom.atomType = line[19]

    return xatom#}}}

def createKPIBond(line) :#{{{
    xbond = KPI.Bond()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = 0
    xbond.bondID = int(line[0])
    if line[1] == "N/A" : line[1] = 0
    xbond.bondNum = int(line[1])
    if line[2] == "N/A" or line[3] == "N/A" : 
        line[2] = 0
        line[3] = 0
    xbond.bond = [int(line[2]), int(line[3])]
    if line[4] == "N/A" or line[5] == "N/A" : 
        line[4] = ''
        line[5] = ''
    xbond.bondType = [line[4], line[5]]
    if line[6] == "N/A" : line[6] = 0.0
    xbond.sigma = float(line[6])
    if line[7] == "N/A" : line[7] = 0.0
    xbond.epsil = float(line[7])

    return xbond#}}}

def createKPIAngle(line) :#{{{
    xangle = KPI.Angle()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = 0
    xangle.angleID = int(line[0])
    if line[1] == "N/A" : line[1] = 0
    xangle.angleNum = int(line[1])
    if line[2] == "N/A" or line[3] == "N/A" or line[4] == "N/A" : 
        line[2] = 0
        line[3] = 0
        line[4] = 0
    xangle.angle = [int(line[2]), int(line[3]), int(line[4])]
    if line[5] == "N/A" or line[6] == "N/A" or line[7] == "N/A" : 
        line[5] = ''
        line[6] = ''
        line[7] = ''
    xangle.angleType = [line[5], line[6], line[7]]
    if line[8] == "N/A" : line[8] = 0.0
    xangle.sigma = float(line[8])
    if line[9] == "N/A" : line[9] = 0.0
    xangle.epsil = float(line[9])

    return xangle#}}}

def createKPIDihedral(line) :#{{{
    xdihedral = KPI.Dihedral()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = 0
    xdihedral.dihedralID = int(line[0])
    if line[1] == "N/A" : line[1] = 0
    xdihedral.dihedralNum = int(line[1])
    if line[2] == "N/A" or line[3] == "N/A" or line[4] == "N/A" or line[5] == "N/A" : 
        line[2] = 0
        line[3] = 0
        line[4] = 0
        line[5] = 0
    xdihedral.dihedral = [int(line[2]), int(line[3]), int(line[4]), int(line[5])]
    if line[6] == "N/A" or line[7] == "N/A" or line[8] == "N/A" or line[9] == "N/A" : 
        line[6] = ''
        line[7] = ''
        line[8] = ''
        line[9] = ''
    xdihedral.dihedralType = [line[6], line[7], line[8], line[9]]
    if line[10] == "N/A" : line[10] = 0.0
    xdihedral.K = float(line[10])
    if line[11] == "N/A" : line[11] = 0
    xdihedral.n = int(line[11])
    if line[12] == "N/A" : line[12] = 0.0
    xdihedral.d = float(line[12])
    if line[13] == "N/A" : line[13] = 0
    xdihedral.w = int(line[13])

    return xdihedral#}}}

def createKPINonbonded(line) :#{{{
    xnonbonded = KPI.Nonbonded()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = 0
    xnonbonded.atomID = int(line[0])
    if line[1] == "N/A" : line[1] = 0
    xnonbonded.atomNum = int(line[1])
    if line[2] == "N/A" : line[2] = '' 
    xnonbonded.atomType = line[2]
    if line[3] == "N/A" : line[3] = 0.0
    xnonbonded.rmin2 = float(line[3])
    if line[4] == "N/A" : line[4] = 0.0
    xnonbonded.epsil = float(line[4])

    return xnonbonded#}}}

def createKPIImproper(line) :#{{{
    ximproper = KPI.Improper()

    # Split line and define atom values
    # all "N/A" entries will define the value as the default value defined in KPI
    line = line.split()
    if line[0] == "N/A" : line[0] = 0
    ximproper.improperID = int(line[0])
    if line[1] == "N/A" : line[1] = 0
    ximproper.improperNum = int(line[1])
    if line[2] == "N/A" or line[3] == "N/A" or line[4] == "N/A" or line[5] == "N/A" : 
        line[2] = 0
        line[3] = 0
        line[4] = 0
        line[5] = 0
    ximproper.improper = [int(line[2]), int(line[3]), int(line[4]), int(line[5])]
    if line[6] == "N/A" or line[7] == "N/A" or line[8] == "N/A" or line[9] == "N/A" : 
        line[6] = ''
        line[7] = ''
        line[8] = ''
        line[9] = ''
    ximproper.improperType = [line[6], line[7], line[8], line[9]]
    if line[10] == "N/A" : line[10] = 0.0
    ximproper.K = float(line[10])
    if line[11] == "N/A" : line[11] = 0.0
    ximproper.chi = float(line[11])

    return ximproper#}}}
#}}}

def readKPREP_kp(fileName) :#{{{
    # Inputs
    
    # Actions

    # Outputs

    # Set-up and Pre-defined variables
    arrAtom = []
    arrBond = []
    arrAngle = []
    arrDihedral = []
    arrNonbonded = []
    arrImproper = []
    arrNBfix = []
    atomHeader = "type\s+atomNum\s+atomSeq\s+atomName\s+atomNameAlt\s+"
    bondHeader = "bondID\s+bondNum\s+bond\s+bondType"
    angleHeader = "angleID\s+angleNum\s+angle\s+angleType"
    dihedralHeader = "dihedralID\s+dihedralNum\s+dihedral\s+dihedralType"
    nonbondedHeader = "atomID\s+atomNum\s+atomType"
    improperHeader = "improperID\s+improperNum\s+improper\s+improperType"
    activeAt = False    # "active Atom fiels" These bool variables inform if the field is actively reading or not
    activeB  = False    # "active Bond field"
    activeA  = False    # "active Angle field"
    activeD  = False    # "active Dihedral field"
    activeN  = False    # "active Nonbonded field"
    activeI  = False    # "active Improper field"
    activeNf = False    # "active NBfix field"

    readFile = open(fileName, 'r')
    for line in readFile :
        if re.match(atomHeader, line) :
            activeAt = True
            continue
        if re.match(bondHeader, line) :
            activeB = True
            continue
        if re.match(angleHeader, line) :
            activeA = True
            continue
        if re.match(dihedralHeader, line) :
            activeD = True
            continue
        if re.match(nonbondedHeader, line) :
            activeN = True
            continue
        if re.match(improperHeader, line) :
            activeI = True
            continue
        if re.match("END", line) :
            activeAt = False
            activeB  = False
            activeA  = False
            activeD  = False
            activeN  = False
            activeI  = False
            activeNf = False
            continue
            
        if activeAt == True :
            arrAtom.append(createKPIAtom(line))
        if activeB == True :
            arrBond.append(createKPIBond(line))
        if activeA == True :
            arrAngle.append(createKPIAngle(line))
        if activeD == True :
            arrDihedral.append(createKPIDihedral(line))
        if activeN == True :
            arrNonbonded.append(createKPINonbonded(line))
        if activeI == True :
            arrImproper.append(createKPIImproper(line))
        #Still need to add NBfix
    
    readFile.close()
    
    return arrAtom, arrBond, arrAngle, arrDihedral, arrNonbonded, arrImproper#}}}
#########################################################


################## readLAMMPS_data ######################
#Functions #{{{
def readLAMMPS_data_init(fileName, printExcess=False) : #{{{
    # Inputs
    # This is where you would add additional sections to read-in. If you do add sections make
    #sure to append what is already coded; do not change the order of returStr
    
    # Actions

    # Outputs

    # Set-up and Pre-defined variables
    #sectionStatus = #{{{
    sectionStatus = {"Header"           : True,
                     "Masses"           : False,
                     "Pair Coeffs"      : False,
                     "Atoms"            : False,
                     "Bond Coeffs"      : False,
                     "Bonds"            : False,
                     "Angle Coeffs"     : False,
                     "Angles"           : False,
                     "Dihedral Coeffs"  : False,
                     "Dihedrals"        : False,
                     "Improper Coeffs"  : False,
                     "Impropers"        : False} #}}}
    
    headerStr = []

    massStr = []
    atomsStr = []

    bondCoeffsStr = []
    bondsStr = []
    
    angleCoeffsStr = []
    anglesStr = []
    
    dihedralCoeffsStr = []
    dihedralsStr = []

    pairCoeffsStr = []
    
    improperCoeffsStr = []
    impropersStr = []

    excessStr = ""

    # Resets the sectionStatus variable and returns any un-used lines
    def sectionStatusReset(line, sectionStatus) :
        change = False
        for index in sectionStatus :
            if re.match(index, line) : 
                for i in sectionStatus : sectionStatus[i] = False
                sectionStatus[index] = True
                change = True
        if   change == False : return line.rstrip()
        elif change == True  : return ''
        
    # Sort data into appropriate section or reset the active section
    readFile = open(fileName, 'r')
    for i, line in enumerate(readFile) :
        if line == '' : continue
        if re.match(spaceFloatPat, line) :
            if   sectionStatus["Header"] == True : headerStr.append(line.rstrip())
            elif sectionStatus["Masses"] == True : massStr.append(line.rstrip())
            elif sectionStatus["Pair Coeffs"]   == True : pairCoeffsStr.append(line.rstrip())
            elif sectionStatus["Atoms"]         == True : atomsStr.append(line.rstrip())
            elif sectionStatus["Bond Coeffs"]   == True : bondCoeffsStr.append(line.rstrip())
            elif sectionStatus["Bonds"]         == True : bondsStr.append(line.rstrip())
            elif sectionStatus["Angle Coeffs"]  == True : angleCoeffsStr.append(line.rstrip())
            elif sectionStatus["Angles"]        == True : anglesStr.append(line.rstrip())
            elif sectionStatus["Dihedral Coeffs"]  == True : dihedralCoeffsStr.append(line.rstrip())
            elif sectionStatus["Dihedrals"]        == True : dihedralsStr.append(line.rstrip())
            elif sectionStatus["Improper Coeffs"]  == True : improperCoeffsStr.append(line.rstrip())
            elif sectionStatus["Impropers"]        == True : impropersStr.append(line.rstrip())
        else :
            excessStr += "line" + str(i+1) + ": " + sectionStatusReset(line, sectionStatus) + "\n"
    readFile.close()

    # Print the lines that did not get sorted into a header or data
    if printExcess == True : print("The following lines were not recorded :\n%s" % excessStr)

    # Return arrays of sorted data
    returnStr =  [headerStr, massStr, pairCoeffsStr, atomsStr, bondCoeffsStr, bondsStr, \
                  angleCoeffsStr, anglesStr, dihedralCoeffsStr, dihedralsStr, \
                  improperCoeffsStr, impropersStr]
    return returnStr #}}}

def readLAMMPS_data_header_charmm(headerStr) : #{{{
    headerData    = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    headerData[0] = headerStr[0].split()[0] #totAtomNum    
    headerData[1] = headerStr[1].split()[0] #totBondNum    
    headerData[2] = headerStr[2].split()[0] #totAngleNum   
    headerData[3] = headerStr[3].split()[0] #totDihedralNum
    headerData[4] = headerStr[4].split()[0] #totImproperNum
    headerData[5] = headerStr[5].split()[0] #totAtomID     
    headerData[6] = headerStr[6].split()[0] #totBondID     
    headerData[7] = headerStr[7].split()[0] #totAngleID    
    headerData[8] = headerStr[8].split()[0] #totDihedralID 
    headerData[9] = headerStr[9].split()[0] #totImproperID 
    headerData[10] = headerStr[10].split()[0] #box dim xlo
    headerData[11] = headerStr[10].split()[1] #box dim xhi 
    headerData[12] = headerStr[11].split()[0] #box dim ylo 
    headerData[13] = headerStr[11].split()[1] #box dim yhi 
    headerData[14] = headerStr[12].split()[0] #box dim zlo 
    headerData[15] = headerStr[12].split()[1] #box dim zhi 
    return headerData 
#}}}

def readLAMMPS_data_parce_atom_molecular(arrAtom, atomsStr, massStr) : #{{{
    if len(arrAtom) != len(atomsStr) : 
        #print(len(arrAtom,), len(atomsStr))
        print("ERROR: number of atoms in header do not match number of atoms in \"Atoms\" field!")
        quit()

    for i, mass in enumerate(massStr) : 
        massStr[i] = mass.split()
        if len(massStr[i]) != 4 : massStr[i] = [massStr[i][0], massStr[i][1], "#", "UnkType"]
    
    for i, atom in enumerate(atomsStr) : atomsStr[i] = atom.split()

    for i, atom in enumerate(arrAtom) :
        coordinates = [float(atomsStr[i][3]), float(atomsStr[i][4]), float(atomsStr[i][5])]
        atomID = int(atomsStr[i][2])
        arrAtom[i].type = "UNDEF"
        arrAtom[i].atomNum = int(atomsStr[i][0])
        arrAtom[i].atomSeq = int(atomsStr[i][0])
        arrAtom[i].resType = "UNK"
        arrAtom[i].coordinates = coordinates
        arrAtom[i].mw = massStr[atomID-1][1]
        arrAtom[i].atomID = atomID
        arrAtom[i].atomType = massStr[atomID-1][3]
#}}}

def readLAMMPS_data_parce_atom_molecular_GO(arrAtom, atomsStr, massStr) : #{{{
    if len(arrAtom) != len(atomsStr) : 
        #print(len(arrAtom,), len(atomsStr))
        print("ERROR: number of atoms in header do not match number of atoms in \"Atoms\" field!")
        quit()

    for i, mass in enumerate(massStr) : 
        massStr[i] = mass.split()
        if len(massStr[i]) != 4 : massStr[i] = [massStr[i][0], massStr[i][1], "#", "UnkType"]
    
    for i, atom in enumerate(atomsStr) : atomsStr[i] = atom.split()

    for i, atom in enumerate(arrAtom) :
        coordinates = [float(atomsStr[i][3]), float(atomsStr[i][4]), float(atomsStr[i][5])]
        atomID = int(atomsStr[i][2])
        arrAtom[i].type = "UNDEF"
        arrAtom[i].atomNum = int(atomsStr[i][0])
        arrAtom[i].atomSeq = int(atomsStr[i][0])
        arrAtom[i].coordinates = coordinates
        arrAtom[i].mw = massStr[atomID-1][1]
        arrAtom[i].atomID = atomID
        arrAtom[i].atomType = massStr[atomID-1][3]
        arrAtom[i].resType = massStr[atomID-1][3]
#}}}

def readLAMMPS_data_parce_bond_harmonic(arrBond, bondsStr, bondCoeffsStr) : #{{{
    if len(arrBond) != len(bondsStr) : 
        print(len(arrBond,), len(bondsStr))
        print("ERROR: number of bonds in header do not match number of atoms in \"Bonds\" field!")
        quit()

    for i, bondCoeff in enumerate(bondCoeffsStr) :
        bondCoeffsStr[i] = bondCoeff.split()
        if len(bondCoeffsStr[i]) != 6 : 
            bondCoeffsStr[i] = [bondCoeffsStr[i][0], bondCoeffsStr[i][1], bondCoeffsStr[i][2], "#", "UnkType", "UnkType"]

    for i, bond in enumerate(bondsStr) : bondsStr[i] = bond.split()

    for i, bond in enumerate(arrBond) :
        bondID = int(bondsStr[i][1])
        bonded = [int(bondsStr[i][2]), int(bondsStr[i][3])]
        bondedType = [bondCoeffsStr[bondID-1][4], bondCoeffsStr[bondID-1][5]]
        arrBond[i].bondID = bondID
        arrBond[i].bondNum = int(bondsStr[i][0])
        arrBond[i].bond = bonded
        arrBond[i].bondType = bondedType
        arrBond[i].sigma = float(bondCoeffsStr[bondID-1][2])
        arrBond[i].epsil = float(bondCoeffsStr[bondID-1][1]) 
#}}}

def readLAMMPS_data_parce_angle_charmm(arrAngle, anglesStr, angleCoeffsStr) : #{{{
    if len(arrAngle) != len(anglesStr) : 
        print(len(arrAngle,), len(anglesStr))
        print("ERROR: number of angles in header do not match number of atoms in \"Angles\" field!")
        quit()

    for i, angleCoeff in enumerate(angleCoeffsStr) :
        angleCoeffsStr[i] = angleCoeff.split()
        if len(angleCoeffsStr[i]) != 9 : 
            angleCoeffsStr[i] = [angleCoeffsStr[i][0], angleCoeffsStr[i][1], angleCoeffsStr[i][2], \
                                 angleCoeffsStr[i][3], angleCoeffsStr[i][4], "#", "UnkType", "UnkType", "UnkType"]

    for i, angle in enumerate(anglesStr) : anglesStr[i] = angle.split()

    for i, angle in enumerate(arrAngle) :
        angleID = int(anglesStr[i][1])
        angled = [int(anglesStr[i][2]), int(anglesStr[i][3]), int(anglesStr[i][4])]
        angledType = [angleCoeffsStr[angleID-1][6], angleCoeffsStr[angleID-1][7], angleCoeffsStr[angleID-1][8]]
        arrAngle[i].angleID = angleID
        arrAngle[i].angleNum = int(anglesStr[i][0])
        arrAngle[i].angle = angled
        arrAngle[i].angleType = angledType
        arrAngle[i].sigma = float(angleCoeffsStr[angleID-1][2])
        arrAngle[i].epsil = float(angleCoeffsStr[angleID-1][1])
#}}}

def readLAMMPS_data_parce_dihedral_charmm(arrDihedral, dihedralsStr, dihedralCoeffsStr) : #{{{
    if len(arrDihedral) != len(dihedralsStr) : 
        print(len(arrDihedral), len(dihedralsStr))
        print("ERROR: number of dihedrals in header do not match number of atoms in \"Dihedrals\" field!")
        quit()

    for i, dihedralCoeff in enumerate(dihedralCoeffsStr) :
        dihedralCoeffsStr[i] = dihedralCoeff.split()
        if len(dihedralCoeffsStr[i]) != 10 : 
            dihedralCoeffsStr[i] = [dihedralCoeffsStr[i][0], dihedralCoeffsStr[i][1], dihedralCoeffsStr[i][2], \
                                    dihedralCoeffsStr[i][3], dihedralCoeffsStr[i][4], \
                                    "#", "UnkType", "UnkType", "UnkType", "UnkType"]

    for i, dihedral in enumerate(dihedralsStr) : dihedralsStr[i] = dihedral.split()

    for i, dihedral in enumerate(arrDihedral) :
        dihedralID = int(dihedralsStr[i][1])
        dihedraled = [int(dihedralsStr[i][2]), int(dihedralsStr[i][3]), int(dihedralsStr[i][4]), int(dihedralsStr[i][5])]
        dihedraledType = [dihedralCoeffsStr[dihedralID-1][6], dihedralCoeffsStr[dihedralID-1][7], \
                          dihedralCoeffsStr[dihedralID-1][8], dihedralCoeffsStr[dihedralID-1][9]]
        arrDihedral[i].dihedralID = dihedralID
        arrDihedral[i].dihedralNum = int(dihedralsStr[i][0])
        arrDihedral[i].dihedral = dihedraled
        arrDihedral[i].dihedralType = dihedraledType
        arrDihedral[i].K  = float(dihedralCoeffsStr[dihedralID-1][1])
        arrDihedral[i].n  = float(dihedralCoeffsStr[dihedralID-1][2])
        arrDihedral[i].d  = float(dihedralCoeffsStr[dihedralID-1][3])
        arrDihedral[i].wt = float(dihedralCoeffsStr[dihedralID-1][4]) 
#}}}

def readLAMMPS_data_parce_dihedral_eten(arrDihedral, dihedralsStr, dihedralCoeffsStr) : #{{{
    if len(arrDihedral) != len(dihedralsStr) : 
        print(len(arrDihedral), len(dihedralsStr))
        print("ERROR: number of dihedrals in header do not match number of atoms in \"Dihedrals\" field!")
        quit()

    for i, dihedralCoeff in enumerate(dihedralCoeffsStr) :
        dihedralCoeffsStr[i] = dihedralCoeff.split()
        if len(dihedralCoeffsStr[i]) != 10 : 
            dihedralCoeffsStr[i] = [dihedralCoeffsStr[i][0], dihedralCoeffsStr[i][1], dihedralCoeffsStr[i][2], \
                                    dihedralCoeffsStr[i][3], dihedralCoeffsStr[i][4], \
                                    "#", "UnkType", "UnkType", "UnkType", "UnkType"]

    for i, dihedral in enumerate(dihedralsStr) : dihedralsStr[i] = dihedral.split()

    for i, dihedral in enumerate(arrDihedral) :
        dihedralID = int(dihedralsStr[i][1])
        dihedraled = [int(dihedralsStr[i][2]), int(dihedralsStr[i][3]), int(dihedralsStr[i][4]), int(dihedralsStr[i][5])]
        dihedraledType = [dihedralCoeffsStr[dihedralID-1][6], dihedralCoeffsStr[dihedralID-1][7], \
                          dihedralCoeffsStr[dihedralID-1][8], dihedralCoeffsStr[dihedralID-1][9]]
        arrDihedral[i].dihedralID = dihedralID
        arrDihedral[i].dihedralNum = int(dihedralsStr[i][0])
        arrDihedral[i].dihedral = dihedraled
        arrDihedral[i].dihedralType = dihedraledType
        arrDihedral[i].K  = float(dihedralCoeffsStr[dihedralID-1][1])
        arrDihedral[i].n  = float(dihedralCoeffsStr[dihedralID-1][2])
        arrDihedral[i].d  = float(dihedralCoeffsStr[dihedralID-1][3])
        arrDihedral[i].wt = float(dihedralCoeffsStr[dihedralID-1][4]) 
#}}}

def readLAMMPS_data_parce_improper_harmonic(arrImproper, impropersStr, improperCoeffsStr) : #{{{
    if len(arrImproper) != len(impropersStr) : 
        print(len(arrImproper), len(impropersStr))
        print("ERROR: number of impropers in header do not match number of atoms in \"Impropers\" field!")
        quit()

    for i, improperCoeff in enumerate(improperCoeffsStr) :
        improperCoeffsStr[i] = improperCoeff.split()
        if len(improperCoeffsStr[i]) != 8 : 
            improperCoeffsStr[i] = [improperCoeffsStr[i][0], improperCoeffsStr[i][1], improperCoeffsStr[i][2], \
                                    "#", "UnkType", "UnkType", "UnkType", "UnkType"]

    for i, improper in enumerate(impropersStr) : impropersStr[i] = improper.split()

    for i, improper in enumerate(arrImproper) :
        improperID = int(impropersStr[i][1])
        impropered = [int(impropersStr[i][2]), int(impropersStr[i][3]), int(impropersStr[i][4]), int(impropersStr[i][5])]
        improperedType = [improperCoeffsStr[improperID-1][4], improperCoeffsStr[improperID-1][5], \
                          improperCoeffsStr[improperID-1][6], improperCoeffsStr[improperID-1][7]]
        arrImproper[i].improperID = improperID
        arrImproper[i].improperNum = int(impropersStr[i][0])
        arrImproper[i].improper = impropered
        arrImproper[i].improperType = improperedType
        arrImproper[i].K   = float(improperCoeffsStr[improperID-1][1])
        arrImproper[i].chi = float(improperCoeffsStr[improperID-1][2])
#}}}

def readLAMMPS_data_parce_pair_lj(arrNonbonded, arrAtom, pairCoeffsStr) : #{{{
    for i, pair in enumerate(pairCoeffsStr) :
        pairCoeffsStr[i] = pair.split()
        if len(pairCoeffsStr[i]) != 5 : 
            pairCoeffsStr[i] = [pairCoeffsStr[i][0], pairCoeffsStr[i][1], pairCoeffsStr[i][2], "#", "UnkType"]

    for i, atom in enumerate(arrNonbonded) :
        atomID   = arrAtom[i].atomID
        atomType = pairCoeffsStr[atomID-1][4]
        atomTypeMasses = arrAtom[i].atomType
        if atomType != atomTypeMasses and atomType != "UnkType" and atomTypeMasses != "UnkType" :
            print("Warning: atom type labeling mis-match between \"Masses\" section and \"Pair Coeff\" section")
            print("Atom ID %i \t Masses: %s \t Pair Coeff: %s" % (atomID, atomTypeMasses, atomType))
        arrNonbonded[i].atomID  = atomID
        arrNonbonded[i].atomNum = arrAtom[i].atomNum
        arrNonbonded[i].atomType = atomType
        arrNonbonded[i].rmin2 = float(pairCoeffsStr[atomID-1][2]) #!!!!!! Do I need to divide by the rmin2 value???
        arrNonbonded[i].epsil = float(pairCoeffsStr[atomID-1][2]) 
#}}}

def readLAMMPS_data_parce_pair_eten(arrNonbonded, arrAtom, pairCoeffsStr) : #{{{
    for i, pair in enumerate(pairCoeffsStr) :
        pairCoeffsStr[i] = pair.split()
        if len(pairCoeffsStr[i]) != 5 : 
            pairCoeffsStr[i] = [pairCoeffsStr[i][0], pairCoeffsStr[i][1], pairCoeffsStr[i][2], "#", "UnkType"]

    for i, atom in enumerate(arrNonbonded) :
        atomID   = arrAtom[i].atomID
        atomType = pairCoeffsStr[atomID-1][4]
        atomTypeMasses = arrAtom[i].atomType
        if atomType != atomTypeMasses and atomType != "UnkType" and atomTypeMasses != "UnkType" :
            print("Warning: atom type labeling mis-match between \"Masses\" section and \"Pair Coeff\" section")
            print("Atom ID %i \t Masses: %s \t Pair Coeff: %s" % (atomID, atomTypeMasses, atomType))
        arrNonbonded[i].atomID  = atomID
        arrNonbonded[i].atomNum = arrAtom[i].atomNum
        arrNonbonded[i].atomType = atomType
        arrNonbonded[i].rmin2 = float(pairCoeffsStr[atomID-1][2]) #!!!!!Do I need to divide by some constant??
        arrNonbonded[i].epsil = float(pairCoeffsStr[atomID-1][2]) 
#}}}
#}}}

def readLAMMPS_data(fileName) : #{{{
    # Inputs
    
    # Actions

    # Outputs

    # Set-up and Pre-defined variables
    arrAtom = []
    arrBond = []
    arrAngle = []
    arrDihedral = []
    arrNonbonded = []
    arrImproper = []
    
    dataStr = readLAMMPS_data_init(fileName, printExcess=False)
    headerStr      = dataStr[0]
    massStr        = dataStr[1]
    pairCoeffsStr  = dataStr[2]
    atomsStr       = dataStr[3]
    bondCoeffsStr  = dataStr[4]
    bondsStr       = dataStr[5]
    angleCoeffsStr = dataStr[6]
    anglesStr      = dataStr[7]
    dihedralCoeffsStr = dataStr[8]
    dihedralsStr      = dataStr[9]
    improperCoeffsStr = dataStr[10]
    impropersStr      = dataStr[11]
    # Check to see if massStr = len(pairCoeffsStr)
    if len(massStr) != len(pairCoeffsStr) :
        print("Error: Number of atoms does not match number or nonbonded interactions!")
        print("Atoms %i \nPair Coeffs %i " % (len(massStr), len(pairCoeffsStr)))
        quit()

    # Define Header information
    headerData = readLAMMPS_data_header_charmm(headerStr)

    # Create Atoms and sort data into proper locations
    for i in range(int(headerData[0])) : arrAtom.append(KPI.Atom())
    readLAMMPS_data_parce_atom_molecular(arrAtom, atomsStr, massStr)

    # Create Bonds and sort data into proper locations
    for i in range(int(headerData[1])) : arrBond.append(KPI.Bond())
    readLAMMPS_data_parce_bond_harmonic(arrBond, bondsStr, bondCoeffsStr)

    # Create Angles and sort data into proper locations
    for i in range(int(headerData[2])) : arrAngle.append(KPI.Angle())
    readLAMMPS_data_parce_angle_charmm(arrAngle, anglesStr, angleCoeffsStr)

    # Create Dihedrals and sort data into proper locations
    for i in range(int(headerData[3])) : arrDihedral.append(KPI.Dihedral())
    readLAMMPS_data_parce_dihedral_charmm(arrDihedral, dihedralsStr, dihedralCoeffsStr)

    # Create Impropers and sort data into proper locations
    for i in range(int(headerData[4])) : arrImproper.append(KPI.Improper())
    readLAMMPS_data_parce_improper_harmonic(arrImproper, impropersStr, improperCoeffsStr)

    # Create Nonbonded and sort pair_coeff data into proper locations
    for i in range(int(headerData[0])) : arrNonbonded.append(KPI.Nonbonded())
    readLAMMPS_data_parce_pair_lj(arrNonbonded, arrAtom, pairCoeffsStr)

    return arrAtom, arrBond, arrAngle, arrDihedral, arrNonbonded, arrImproper
#}}}

def readLAMMPS_data_GO(fileName) : #{{{
    # Inputs
    
    # Actions

    # Outputs

    # Set-up and Pre-defined variables
    arrAtom = []
    arrBond = []
    arrAngle = []
    arrDihedral = []
    arrNonbonded = []
    arrImproper = []
    
    dataStr = readLAMMPS_data_init(fileName, printExcess=False)
    headerStr      = dataStr[0]
    massStr        = dataStr[1]
    pairCoeffsStr  = dataStr[2]
    atomsStr       = dataStr[3]
    bondCoeffsStr  = dataStr[4]
    bondsStr       = dataStr[5]
    angleCoeffsStr = dataStr[6]
    anglesStr      = dataStr[7]
    dihedralCoeffsStr = dataStr[8]
    dihedralsStr      = dataStr[9]
    improperCoeffsStr = dataStr[10]
    impropersStr      = dataStr[11]
    # Check to see if massStr = len(pairCoeffsStr)
    if len(massStr) != len(pairCoeffsStr) :
        print("Error: Number of atoms does not match number or nonbonded interactions!")
        print("Atoms %i \nPair Coeffs %i " % (len(massStr), len(pairCoeffsStr)))
        quit()

    # Define Header information
    headerData = readLAMMPS_data_header_charmm(headerStr)

    # Create Atoms and sort data into proper locations
    for i in range(int(headerData[0])) : arrAtom.append(KPI.Atom())
    readLAMMPS_data_parce_atom_molecular_GO(arrAtom, atomsStr, massStr)

    # Create Bonds and sort data into proper locations
    for i in range(int(headerData[1])) : arrBond.append(KPI.Bond())
    readLAMMPS_data_parce_bond_harmonic(arrBond, bondsStr, bondCoeffsStr)

    # Create Angles and sort data into proper locations
    for i in range(int(headerData[2])) : arrAngle.append(KPI.Angle())
    readLAMMPS_data_parce_angle_charmm(arrAngle, anglesStr, angleCoeffsStr)

    # Create Dihedrals and sort data into proper locations
    for i in range(int(headerData[3])) : arrDihedral.append(KPI.Dihedral())
    readLAMMPS_data_parce_dihedral_eten(arrDihedral, dihedralsStr, dihedralCoeffsStr)

    # Create Impropers and sort data into proper locations
    for i in range(int(headerData[4])) : arrImproper.append(KPI.Improper())
    readLAMMPS_data_parce_improper_harmonic(arrImproper, impropersStr, improperCoeffsStr)

    # Create Nonbonded and sort pair_coeff data into proper locations
    for i in range(int(headerData[0])) : arrNonbonded.append(KPI.Nonbonded())
    readLAMMPS_data_parce_pair_eten(arrNonbonded, arrAtom, pairCoeffsStr)

    return arrAtom, arrBond, arrAngle, arrDihedral, arrNonbonded, arrImproper
#}}}
#########################################################


################### readLAMMPS_nc #######################
def readLAMMPS_nc(fileName, arrAtom) : #{{{
    arrNBfix = []
    readFile = open(fileName, 'r')
    for i, line in enumerate(readFile) :
        line = line.split()
        if len(line) != 5 or line[0] != "pair_coeff" :
            print("Error: This is not set-up like a proper .nc.mod file")
            print("Issue first encountered on line %i" % (i+1))
            quit()
        p1 = int(line[1]) # first atom in pair
        p2 = int(line[2]) # second atom in pair
        if p1 > p2 : p1, p2 = p2, p1
        p1type = arrAtom[p1-1].resType
        if p1type == "" : p1type = "UnkType"
        p2type = arrAtom[p2-1].resType
        if p2type == "" : p2type = "UnkType"

        arrNBfix.append(KPI.NBfix())
        arrNBfix[i].nbfixID = i+1
        arrNBfix[i].nbfix = [p1, p2]
        arrNBfix[i].nbfixType=[p1type, p2type]
        arrNBfix[i].sigma = float(line[4])
        arrNBfix[i].epsil = float(line[3])
    readFile.close()
    return arrNBfix
#}}}
#########################################################


