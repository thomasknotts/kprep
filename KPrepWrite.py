#
# KPrepWrite.py
# --------------------
# Authors: Addison Smith and Derek Bush
# Contributors: Anthony Gillepsie
# Written to be run on python 2 or 3

#########################################################
#--------------- Index of Function ----------------------
# write_pdb            Writes .pdb files directly from the protein data bank website
# writeGO_pdb           Writes .pdb files directly from the protein data bank website
# writeLAMMPS_data     Writes .data files used in LAMMPS for all-atom simulations
# writeLAMMPS_GO_data  Writes .data files used in LAMMPS for G0-model simulations
# writeLAMMPS_GO_nc    Writes .mod files that contain Native Contact (nc) (ie. NBFIX) information
# writeLAMMPS_GO_aa    Writes .mod files that contain Amino Acid (aa) information
# writeLAMMPS_GO_wc    Writes .mod files that contain HISAA Wall Coefficient (wc) information
# writeKPREP_kp        Writes a .kp file, this is just a data-dump of all KPI variables
#--------------------------------------------------------
#########################################################


#########################################################
import KPrepInfo as KPI
import re
#########################################################


##################### write_pdb #########################
def write_pdb(fileNameNew, arrAtom, reset=False) :#{{{
    # Inputs:
    # fileNameNew          This is the name of the file that will be created
    # arrAtom              This is an arrays of class variables. Make sure that you have all needed
    #                         data for the write
    # reset                Binary True/False. If True, the PDB will reset the numbering so that atom and
    #                         residue values start at 1 (commonly needed for CHARMM, KSS and other programs)

    # Actons:
    # Will produce a PDB file.

    # Outputs:
    # Nothing is returned, but the fileNameNew file is created and prepended with KPW to avoid over-writes.

    # Set-up and Pre-define variables
    totAtom = len(arrAtom)

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writePDB = open(fileNameNew, 'w')

    # Error checks
    # End check
    
    # Write all Atom information in PDB format 
    line = "REMARK Created by KPrepWrite\n"
    writePDB.write(line)

    # Loop over each atom, print info in its row
    for i in range(totAtom) :
        if reset == True : line = "%6s%5i %4s%1s%3s %1s%4i    %8.3f%8.3f%8.3f%6s%6s           %2s%2s\n" % \
        (arrAtom[i].type, \
         arrAtom[i].atomNum, \
         arrAtom[i].atomName, \
         arrAtom[i].atomNameAlt, \
         arrAtom[i].resType, \
         arrAtom[i].chainType, \
         arrAtom[i].resNum, \
         arrAtom[i].coordinates[0], \
         arrAtom[i].coordinates[1], \
         arrAtom[i].coordinates[2], \
         str(arrAtom[i].o), \
         str(arrAtom[i].tf), \
         arrAtom[i].element, \
         arrAtom[i].q)

        elif reset == False : line = "%6s%5i %4s%1s%3s %1s%4i    %8.3f%8.3f%8.3f%6s%6s           %2s%2s\n" % \
        (arrAtom[i].type, \
         arrAtom[i].atomSeq, \
         arrAtom[i].atomName, \
         arrAtom[i].atomNameAlt, \
         arrAtom[i].resType, \
         arrAtom[i].chainType, \
         arrAtom[i].resSeq, \
         arrAtom[i].coordinates[0], \
         arrAtom[i].coordinates[1], \
         arrAtom[i].coordinates[2], \
         str(arrAtom[i].o), \
         str(arrAtom[i].tf), \
         arrAtom[i].element, \
         arrAtom[i].q)

        writePDB.write(line)

    line = "END                                                                             "
    writePDB.write(line)

    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################### writeGO_pdb #########################
def writeGO_pdb(fileNameNew, arrAtom) :#{{{
    # Inputs:
    # fileNameNew           This is the name of the file that will be created
    # arrAtom               This is an arrays of class variables. Make sure that you have all needed
    #                         data for the write

    # Actons:
    # Will produce a PDB file in the Go Model format
    #    DO NOT mix with all-atom data, will NOT produce an appropriate all-atom pdb

    # Outputs:
    # Nothing is returned, but the fileNameNew file is created and prepended with KPW to avoid over-writes

    # Set-up and Pre-define variables
    totAtom = len(arrAtom)

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writePDB = open(fileNameNew, 'w')

    # Go-model format check
    for i in range(totAtom) :
        if arrAtom[i].atomNum != arrAtom[i].resNum :
            print("\nWARNING: atomNum and resNum do not match, Go-model PDB might be compromised...\n")
            break
    # End check
    
    # Write all Atom information in PEG format; Assume Go-model format: 
    #    atomNum = resNum; and base-1 numbering
    line = "REMARK Created by KPrepWrite\n"
    writePDB.write(line)

    # Loop over each atom, print info in its row
    for i in range(totAtom) :
        line = "%6s%5i %4s%1s%3s %1s%4i    %8.3f%8.3f%8.3f  1.00  1.00      PROT    \n" % \
        (arrAtom[i].type, \
         arrAtom[i].atomNum, \
         arrAtom[i].atomName, \
         arrAtom[i].atomNameAlt, \
         arrAtom[i].resType, \
         arrAtom[i].chainType, \
         arrAtom[i].resNum, \
         arrAtom[i].coordinates[0], \
         arrAtom[i].coordinates[1], \
         arrAtom[i].coordinates[2])
        writePDB.write(line)

    line = "END                                                                             "
    writePDB.write(line)

    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################# write***_data ######################
# Functions #{{{
def header_data(arrHeader, title="Created by KPrep") : #{{{
    headerStr = """\
%s

\t%d\t\tatoms
\t%d\t\tbonds
\t%d\t\tangles
\t%d\t\tdihedrals
\t%d\t\timpropers

\t%d\t\tatom types
\t%d\t\tbond types
\t%d\t\tangle types
\t%d\t\tdihedral types
\t%d\t\timproper types

\t%.3f\t\t%.3f\t\txlo xhi
\t%.3f\t\t%.3f\t\tylo yhi
\t%.3f\t\t%.3f\t\tzlo zhi
        """ % (title, arrHeader[0],  arrHeader[1],  arrHeader[2],  arrHeader[3],  arrHeader[4],
               arrHeader[5],  arrHeader[6],  arrHeader[7],  arrHeader[8],  arrHeader[9],
               arrHeader[10], arrHeader[11], arrHeader[12], arrHeader[13], arrHeader[14], arrHeader[15])
    return headerStr #}}}

def masses_data(arrAtom, totAtomNum) : #{{{
    counter = 0
    massStr = "\nMasses\n\n"
    for i in range(totAtomNum) :
        if counter < arrAtom[i].atomID :
            a = str(arrAtom[i].atomID).ljust(5)
            b = str("%.2f" % (arrAtom[i].mw)).rjust(8)
            c = arrAtom[i].atomType
            massStr += "    %s  %s  # %s\n" % (a, b, c)
            counter += 1
    return massStr #}}}

def massesGO_data(arrAtom, totAtomNum) : #{{{
    massStr = "\nMasses\n\n"
    for i in range(totAtomNum) :
        a = str(arrAtom[i].atomNum).ljust(5)
        b = str("%.2f" % (arrAtom[i].mw)).rjust(8)
        c = arrAtom[i].resType
        massStr += "    %s  %s  # %s\n" % (a, b, c)
    return massStr #}}}

def paircoeff_data(arrNonbonded, totNonbondedNum) : #{{{
    counter = 0
    paircoeffStr = "\nPair Coeffs # lj/charmm/coul/long\n\n"
    for i in range(totNonbondedNum) :
        if arrNonbonded[i].epsil < 0 : arrNonbonded[i].epsil = -arrNonbonded[i].epsil
        if counter < arrNonbonded[i].atomID :
            a = str(arrNonbonded[i].atomID).ljust(5)
            b = str("%.6f" % (arrNonbonded[i].epsil)).rjust(10)
            c = str("%.6f" % (arrNonbonded[i].rmin2 * KPI.nbRminToSigma)).rjust(10) # mult by constant for proper form
            d = arrNonbonded[i].atomType
            line = "    %s  %s  %s  # %s\n" % (a, b, c, d)
            paircoeffStr += line
            counter += 1
    return paircoeffStr #}}}

def paircoeffGO_data(arrNonbonded, arrAtom, totNonbondedNum) : #{{{
    paircoeffStr = "\nPair Coeffs # lj/charmm/eten\n\n"
    for i in range(totNonbondedNum) :
        if arrNonbonded[i].epsil < 0 : arrNonbonded[i].epsil = -arrNonbonded[i].epsil
        a = str(arrNonbonded[i].atomNum).ljust(5)
        b = str("%.6f" % (arrNonbonded[i].epsil)).rjust(10)
        c = str("%.6f" % (arrNonbonded[i].rmin2 * KPI.nbRminToEten)).rjust(10) # mult by constant for proper form
        d = arrAtom[arrNonbonded[i].atomNum-1].resType
        line = "    %s  %s  %s # %s\n" % (a, b, c, d)
        paircoeffStr += line
    return paircoeffStr #}}}

def atommolecular_data(arrAtom, totAtomNum) : #{{{
    atomStr = "\nAtoms # molecular\n\n"
    for i in range(totAtomNum) :
        a = str(arrAtom[i].atomNum).ljust(5)
        b = str(arrAtom[i].resNum).ljust(5)
        c = str(arrAtom[i].atomID).ljust(5)
        x = str("%.3f" % (arrAtom[i].coordinates[0])).rjust(12)
        y = str("%.3f" % (arrAtom[i].coordinates[1])).rjust(12)
        z = str("%.3f" % (arrAtom[i].coordinates[2])).rjust(12)
        line = "    %s  %s  %s  %s  %s  %s\n" % (a, b, c, x, y, z)
        atomStr += line
    return atomStr #}}}

def atommolecularGO_data(arrAtom, totAtomNum) : #{{{
    atomStr = "\nAtoms # molecular\n\n"
    for i in range(totAtomNum) :
        a = str(arrAtom[i].atomNum).ljust(5)
        b = str(1).ljust(5) # Go-model proteins are all the same moleucle
        c = str(arrAtom[i].atomNum).ljust(5)
        x = str("%.3f" % (arrAtom[i].coordinates[0])).rjust(12)
        y = str("%.3f" % (arrAtom[i].coordinates[1])).rjust(12)
        z = str("%.3f" % (arrAtom[i].coordinates[2])).rjust(12)
        line = "    %s  %s  %s  %s  %s  %s\n" % (a, b, c, x, y, z)
        atomStr += line
    return atomStr #}}}

def atomfull_data(arrAtom, totAtomNum) : #{{{
    atomStr = "\nAtoms # full\n\n"
    for i in range(totAtomNum) :
        a = str(arrAtom[i].atomNum).ljust(5)
        b = str(arrAtom[i].resNum).ljust(5)
        c = str(arrAtom[i].atomID).ljust(5)
        d = str(arrAtom[i].q).rjust(5)
        x = str("%.3f" % (arrAtom[i].coordinates[0])).rjust(12)
        y = str("%.3f" % (arrAtom[i].coordinates[1])).rjust(12)
        z = str("%.3f" % (arrAtom[i].coordinates[2])).rjust(12)
        line = "    %s  %s  %s  %s  %s  %s  %s\n" % (a, b, c, d, x, y, z)
        atomStr += line
    return atomStr #}}}

def bondcoeff_data(arrBond, totBondNum) : #{{{
    counter = 0
    bondcoeffStr = "\nBond Coeffs\n\n"
    for i in range(totBondNum) :
        if counter < arrBond[i].bondID :
            a = str(arrBond[i].bondID).ljust(5)
            b = str("%.2f" % (arrBond[i].epsil)).rjust(7)
            c = str("%.6f" % (arrBond[i].sigma)).rjust(7)
            d = "%s %s" % (arrBond[i].bondType[0], arrBond[i].bondType[1])
            line = "    %s    %s    %s  # %s\n" % (a, b, c, d)
            bondcoeffStr += line
            counter += 1
    return bondcoeffStr #}}}

def bondcoeffGO_data(arrBond, totBondNum) : #{{{
    bondcoeffStr = "\nBond Coeffs\n\n"
    for i in range(totBondNum) :
        a = str(arrBond[i].bondNum).ljust(5)
        b = str("%.2f" % (arrBond[i].epsil)).rjust(7)
        c = str("%.6f" % (arrBond[i].sigma)).rjust(7)
        d = "%s %s" % (str(arrBond[i].bond[0]), str(arrBond[i].bond[1]))
        line = "    %s    %s    %s  # %s\n" % (a, b, c, d)
        bondcoeffStr += line
    return bondcoeffStr #}}}

def bond_data(arrBond, totBondNum) : #{{{
    bondStr = "\nBonds\n\n"
    for i in range(totBondNum) :
        a = str(arrBond[i].bondNum).ljust(5)
        b = str(arrBond[i].bondID).ljust(5)
        x = str(arrBond[i].bond[0]).ljust(5)
        y = str(arrBond[i].bond[1]).ljust(5)
        line = "    %s    %s    %s    %s\n" % (a, b, x, y)
        bondStr += line
    return bondStr #}}}

def anglecoeff_data(arrAngle, totAngleNum) : #{{{
    counter = 0
    anglecoeffStr = "\nAngle Coeffs\n\n"
    for i in range(totAngleNum) :
        if counter < arrAngle[i].angleID :
            a = str(arrAngle[i].angleID).ljust(5)
            b = str("%.2f" % (arrAngle[i].epsil)).rjust(7)
            c = str("%.6f" % (arrAngle[i].sigma)).rjust(12) 
            d = str(0).ljust(3)  # !!!!! I do not know what this value means
            e = str(0).ljust(3)  # !!!!! I do not know what this value means
            f = "%s %s %s" % (arrAngle[i].angleType[0], arrAngle[i].angleType[1], arrAngle[i].angleType[2])
            line = "    %s    %s %s    %s  %s  # %s\n" % (a, b, c, d, e, f)
            anglecoeffStr += line
            counter += 1
    return anglecoeffStr #}}}

def anglecoeffGO_data(arrAngle, totAngleNum) : #{{{
    anglecoeffStr = "\nAngle Coeffs\n\n"
    for i in range(totAngleNum) :
        a = str(arrAngle[i].angleNum).ljust(5)
        b = str("%.2f" % (arrAngle[i].epsil)).rjust(7)
        c = str("%.6f" % (arrAngle[i].sigma)).rjust(12) 
        d = str(0).ljust(3)  # !!!!! I do not know what this value means
        e = str(0).ljust(3)  # !!!!! I do not know what this value means
        f = "%s %s %s" % (str(arrAngle[i].angle[0]), str(arrAngle[i].angle[1]), str(arrAngle[i].angle[2]))
        line = "    %s    %s %s    %s  %s  # %s\n" % (a, b, c, d, e, f)
        anglecoeffStr += line
    return anglecoeffStr #}}}

def angle_data(arrAngle, totAngleNum) : #{{{
    angleStr = "\nAngles\n\n"
    for i in range(totAngleNum) :
        a = str(arrAngle[i].angleNum).ljust(5)
        b = str(arrAngle[i].angleID).ljust(5)
        x = str(arrAngle[i].angle[0]).ljust(5)
        y = str(arrAngle[i].angle[1]).ljust(5)
        z = str(arrAngle[i].angle[2]).ljust(5)
        line = "    %s    %s    %s    %s    %s\n" % (a, b, x, y, z)
        angleStr += line
    return angleStr #}}}

def dihedralcoeff_data(arrDihedral, totDihedralNum) : #{{{
    counter = 0
    if totDihedralNum == 0 : return ""
    dihedralcoeffStr = "\nDihedral Coeffs\n\n"
    for i in range(totDihedralNum) :
        if counter < arrDihedral[i].dihedralID :
            a = str(arrDihedral[i].dihedralID).ljust(5)
            b = str("%.6f" % (arrDihedral[i].K)).rjust(10)
            c = str(arrDihedral[i].n).ljust(3)
            d = str("%i" % (arrDihedral[i].d)).rjust(10)
            e = str(arrDihedral[i].wt).ljust(3)
            f = "%s %s %s %s" % (arrDihedral[i].dihedralType[0], arrDihedral[i].dihedralType[1], \
                                 arrDihedral[i].dihedralType[2], arrDihedral[i].dihedralType[3])
            line = "    %s  %s    %s  %s    %s  # %s\n" % (a, b, c, d, e, f)
            dihedralcoeffStr += line
            counter += 1
    return dihedralcoeffStr #}}}

def dihedralcoeffGO_data(arrDihedral, totDihedralNum) : #{{{
    if totDihedralNum == 0 : return ""
    dihedralcoeffStr = "\nDihedral Coeffs\n\n"
    for i in range(totDihedralNum) :
        a = str(arrDihedral[i].dihedralNum).ljust(5)
        b = str("%.6f" % (arrDihedral[i].K)).rjust(10)
        c = str(arrDihedral[i].n).ljust(3)
        d = str("%i" % (arrDihedral[i].d)).rjust(10)
        e = str(arrDihedral[i].wt).ljust(3)
        f = "%s %s %s %s" % (str(arrDihedral[i].dihedral[0]), str(arrDihedral[i].dihedral[1]), \
                             str(arrDihedral[i].dihedral[2]), str(arrDihedral[i].dihedral[3]))
        line = "    %s  %s    %s  %s    %s  # %s\n" % (a, b, c, d, e, f)
        dihedralcoeffStr += line
    return dihedralcoeffStr #}}}

def dihedral_data(arrDihedral, totDihedralNum) : #{{{
    if totDihedralNum == 0 : return ""
    dihedralStr = "\nDihedrals\n\n"
    for i in range(totDihedralNum) :
        a = str(arrDihedral[i].dihedralNum).ljust(5)
        b = str(arrDihedral[i].dihedralID).ljust(5)
        x = str(arrDihedral[i].dihedral[0]).ljust(5)
        y = str(arrDihedral[i].dihedral[1]).ljust(5)
        z = str(arrDihedral[i].dihedral[2]).ljust(5)
        w = str(arrDihedral[i].dihedral[3]).ljust(5)
        line = "    %s    %s    %s    %s    %s    %s\n" % (a, b, x, y, z, w)
        dihedralStr += line
    return dihedralStr #}}}

def impropercoeff_data(arrImproper, totImproperNum) : #{{{
    counter = 0
    if totImproperNum == 0 : return ""
    impropercoeffStr = "\nImproper Coeffs\n\n"
    for i in range(totImproperNum) :
        if counter < arrImproper[i].improperID :
            a = str(arrImproper[i].improperID).ljust(5)
            b = str("%.6f" % (arrImproper[i].K)).rjust(10)
            c = str("%.6f" % (arrImproper[i].chi)).rjust(10)
            d = "%s %s %s %s" % (arrImproper[i].improperType[0], arrImproper[i].improperType[1], \
                                 arrImproper[i].improperType[2], arrImproper[i].improperType[3])
            line = "    %s    %s    %s  # %s\n" % (a, b, c, d)
            impropercoeffStr += line
            counter += 1
    return impropercoeffStr #}}}

def improper_data(arrImproper, totImproperNum) : #{{{
    if totImproperNum == 0 : return ""
    improperStr = "\nImpropers\n\n"
    for i in range(totImproperNum) :
        a = str(arrImproper[i].improperNum).ljust(5)
        b = str(arrImproper[i].improperID).ljust(5)
        x = str(arrImproper[i].improper[0]).ljust(5)
        y = str(arrImproper[i].improper[1]).ljust(5)
        z = str(arrImproper[i].improper[2]).ljust(5)
        w = str(arrImproper[i].improper[3]).ljust(5)
        line = "    %s    %s    %s    %s    %s    %s\n" % (a, b, x, y, z, w)
        improperStr += line
    return improperStr #}}}
#}}}

#{{{ writeLAMMPS_data
def writeLAMMPS_data(fileNameNew, arrAtom=[], arrBond=[], arrAngle=[], arrDihedral=[], \
                     arrImproper=[], arrNonbonded=[], boxMin=-791, boxMax=791) :
    # Inputs:
    # fileNameNew           This is the name of the file that will be created
    # arrAtom-arrImproper   These are arrays of class variables. Make sure that you have all needed
    #                         data for the write
    # boxMin, boxMax        These variables define the dimensions of a cubic simulation box

    # Actions:
    # Will take all the input data and produce a LAMMPS-compatable data file
    # If data is neglected it can produce a partial data file you can copy/paste into another file should you want.
    # For reference see read_data command in literature for specifics on *_style requirements
    # NOTE: boxMin/Max are important; when run, if it hangs-up give it a more reasonable box size
    #       -> Add functionality where it just does max/min ranges of the coordinates as box ranges

    # Outputs:
    # Nothing is returned, but fileNameNew is created and populated with relavent information.

    # Set-up and Pre-define variables
    #{{{
    chainDict = dict(A=1, B=2, C=3, D=4, E=5, F=6, G=7, H=8, \
                     I=9, J=10, K=11, L=12, M=13, N=14, O=15)
    totAtomNum = len(arrAtom)
    totAtomID  = 0 
    for j in range(totAtomNum) :
        if totAtomID < arrAtom[j].atomID : totAtomID = arrAtom[j].atomID

    totBondNum = len(arrBond)
    totBondID  = 0
    for j in range(totBondNum) :
        if totBondID < arrBond[j].bondID : totBondID = arrBond[j].bondID

    totAngleNum = len(arrAngle)
    totAngleID  = 0
    for j in range(totAngleNum) :
        if totAngleID < arrAngle[j].angleID : totAngleID = arrAngle[j].angleID

    totDihedralNum = len(arrDihedral)
    totDihedralID  = 0
    for j in range(totDihedralNum) :
        if totDihedralID < arrDihedral[j].dihedralID : totDihedralID = arrDihedral[j].dihedralID

    totNonbondedNum = len(arrNonbonded)

    totImproperNum = len(arrImproper)
    totImproperID  = 0
    for j in range(totImproperNum) :
        if totImproperID < arrImproper[j].improperID : totImproperID = arrImproper[j].improperID#}}}
    
    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeData = open(fileNameNew, 'w')

    # Write header for data file 
    print("Writing header...")
    arrHeader = [totAtomNum, totBondNum, totAngleNum, totDihedralNum, totImproperNum,
                 totAtomID,  totBondID,  totAngleID,  totDihedralID,  totImproperID,
                 boxMin, boxMax, boxMin, boxMax, boxMin, boxMax]
    headerStr = header_data(arrHeader)
    writeData.write(headerStr)

    # Write "Masses" to data file
    print("Writing 'Masses' section...")
    massStr = masses_data(arrAtom, totAtomNum)
    writeData.write(massStr)

    # Write "Pair Coeff" to data file
    print("Writing 'Pair Coeff' section...")
    paircoeffStr = paircoeff_data(arrNonbonded, totNonbondedNum)
    writeData.write(paircoeffStr)

    # Write "Atoms" (atom coordinates) to data file
    # Future, maybe include "atom_style full" that specifies charge before the coordinates?
    print("Writing 'Atoms' section...")
    atomStr = atommolecular_data(arrAtom, totAtomNum)
    writeData.write(atomStr)

    # Write "Bond Coeffs" to data file
    print("Writing 'Bond Coeffs' section...")
    bondcoeffStr = bondcoeff_data(arrBond, totBondNum)
    writeData.write(bondcoeffStr)

    # Write "Bonds" to data file
    print("Writing 'Bonds' section...")
    bondStr = bond_data(arrBond, totBondNum)
    writeData.write(bondStr)
         
    # Write "Angle Coeffs" to data file
    print("Writing 'Angle Coeffs' section...")
    anglecoeffStr = anglecoeff_data(arrAngle, totAngleNum)
    writeData.write(anglecoeffStr)

    # Write "Angles" to data file
    print("Writing 'Angles' section...")
    angleStr = angle_data(arrAngle, totAngleNum)
    writeData.write(angleStr)

    # Write "Dihedral Coeffs" to data file
    print("Writing 'Dihedral Coeffs' section...")
    dihedralcoeffStr = dihedralcoeff_data(arrDihedral, totDihedralNum)
    writeData.write(dihedralcoeffStr)

    # Write "Dihedrals" to data file
    print("Writing 'Dihedrals' section...")
    dihedralStr = dihedral_data(arrDihedral, totDihedralNum)
    writeData.write(dihedralStr)

    # Write "Improper Coeffs" to data file
    print("Writing 'Improper Coeffs' section...")
    impropercoeffStr = impropercoeff_data(arrImproper, totImproperNum)
    writeData.write(impropercoeffStr)

    # Write "Impropers" to data file
    print("Writing 'Impropers' section...")
    improperStr = improper_data(arrImproper, totImproperNum)
    writeData.write(improperStr)

    writeData.close()
    print("Finished writing %s" % (fileNameNew))#}}}

#{{{ writeLAMMPS_GO_data
def writeLAMMPS_GO_data(fileNameNew, arrAtom=[], arrBond=[], arrAngle=[], arrDihedral=[], \
                     arrImproper=[], arrNonbonded=[], boxMin=-791, boxMax=791) :
    # Inputs:
    # fileNameNew           This is the name of the file that will be created
    # arrAtom-arrImproper   These are arrays of class variables. Make sure that you have all needed
    #                         data for the write
    # boxMin, boxMax        These variables define the dimensions of a cubic simulation box

    # Actions:
    # Will take all the input data and produce a LAMMPS-compatable data file
    # If data is neglected it can produce a partial data file you can copy/paste into another file should you want.
    # For reference see read_data command in literature for specifics on *_style requirements
    # NOTE: boxMin/Max are important; when run, if it hangs-up give it a more reasonable box size
    #       -> Add functionality where it just does max/min ranges of the coordinates as box ranges

    # Outputs:
    # Nothing is returned, but fileNameNew is created and populated with relavent information.

    # Set-up and Pre-define variables
    #{{{
    chainDict = dict(A=1, B=2, C=3, D=4, E=5, F=6, G=7, H=8, \
                     I=9, J=10, K=11, L=12, M=13, N=14, O=15)
    totAtomNum = len(arrAtom)
    totAtomID  = totAtomNum

    totBondNum = len(arrBond)
    totBondID  = totBondNum

    totAngleNum = len(arrAngle)
    totAngleID  = totAngleNum

    totDihedralNum = len(arrDihedral)
    totDihedralID  = totDihedralNum

    totImproperNum = 0
    totImproperID  = 0

    totNonbondedNum = len(arrNonbonded)
    #}}}

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeData = open(fileNameNew, 'w')

    # Write header for data file 
    print("Writing header...")
    arrHeader = [totAtomNum, totBondNum, totAngleNum, totDihedralNum, totImproperNum,
                 totAtomID,  totBondID,  totAngleID,  totDihedralID,  totImproperID,
                 boxMin, boxMax, boxMin, boxMax, boxMin, boxMax]
    headerStr = header_data(arrHeader, title="Created by KPrep for a Go-model system")
    writeData.write(headerStr)

    # Write "Masses" to data file
    print("Writing 'Masses' section...")
    massStr = massesGO_data(arrAtom, totAtomNum)
    writeData.write(massStr)

    # Write "Pair Coeff" to data file
    print("Writing 'Pair Coeff' section...")
    paircoeffStr = paircoeffGO_data(arrNonbonded, arrAtom, totNonbondedNum)
    writeData.write(paircoeffStr)

    # Write "Atoms" (atom coordinates) to data file
    # Future, maybe include "atom_style full" that specifies charge before the coordinates?
    print("Writing 'Atoms' section...")
    atomStr = atommolecularGO_data(arrAtom, totAtomNum)
    writeData.write(atomStr)

    # Write "Bond Coeffs" to data file
    print("Writing 'Bond Coeffs' section...")
    bondcoeffStr = bondcoeffGO_data(arrBond, totBondNum)
    writeData.write(bondcoeffStr)

    # Write "Bonds" to data file
    print("Writing 'Bonds' section...")
    bondStr = bond_data(arrBond, totBondNum)
    writeData.write(bondStr)
         
    # Write "Angle Coeffs" to data file
    print("Writing 'Angle Coeffs' section...")
    anglecoeffStr = anglecoeffGO_data(arrAngle, totAngleNum)
    writeData.write(anglecoeffStr)

    # Write "Angles" to data file
    print("Writing 'Angles' section...")
    angleStr = angle_data(arrAngle, totAngleNum)
    writeData.write(angleStr)

    # Write "Dihedral Coeffs" to data file
    print("Writing 'Dihedral Coeffs' section...")
    dihedralcoeffStr = dihedralcoeffGO_data(arrDihedral, totDihedralNum)
    writeData.write(dihedralcoeffStr)

    # Write "Dihedrals" to data file
    print("Writing 'Dihedrals' section...")
    dihedralStr = dihedral_data(arrDihedral, totDihedralNum)
    writeData.write(dihedralStr)

    writeData.close()
    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################# writeLAMMPS_GO_nc #####################
def writeLAMMPS_GO_nc(fileNameNew, arrNBfix) :#{{{
    # Inputs:
    # fileNameNew           This must be the file name that will be created that contains the native contact data
    # arrNbfix              This is an arrays of class variables. Make sure that you have all needed
    #                         data for the write

    # Actons:
    # Will produce a list of LAMMPS pair_coeff commands that enable native contact interactions

    # Outputs:
    # Nothing is returned, but the fileNameNc file is generated.

    # Set-up and Pre-define variables
    totNBfixNum = len(arrNBfix)

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeData = open(fileNameNew, 'w')

    # Wrtie pair_coeff NBfix modifications
    print("Writing 'NBfix mod file' section...")
    for i in range(totNBfixNum) :
        x = str(arrNBfix[i].nbfix[0]).ljust(5)
        y = str(arrNBfix[i].nbfix[1]).ljust(5)
        a = str("%.6f" % (arrNBfix[i].epsil)).rjust(12)
        b = str("%.6f" % (arrNBfix[i].sigma)).rjust(12)
        line = ("pair_coeff      %s    %s  %s  %s\n" % (x, y, a, b))
        writeData.write(line)

    writeData.close()
    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################# writeLAMMPS_GO_aa #####################
def writeLAMMPS_GO_aa(fileNameNew, arrAtom) :#{{{
    # Inputs:
    # fileNameNew           This must be the file name that will be created that contains the native contact data
    # arrAtom               This is an arrays of class variables. Make sure that you have all needed
    #                         data for the write

    # Actons:
    # Will produce a list of LAMMPS set commands that encode each amino acid for hydrophobic interaction

    # Outputs:
    # Nothing is returned, but the fileNameNew file is generated.

    # Set-up and Pre-define variables
    totAA = len(arrAtom)

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeData = open(fileNameNew, 'w')

    # Wrtie pair_coeff NBfix modifications
    print("Writing 'HISAA amino acid mod file' section...")
    for i in range(totAA) :
        a = str("%i" % (arrAtom[i].atomNum)).ljust(4)
        b = str("%i" % (KPI.dictAA[arrAtom[i].resType])).ljust(4)
        line = ("set             atom  %s  i_AAcid  %s\n" % (a, b))
        writeData.write(line)

    writeData.close()
    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################# writeLAMMPS_GO_wc #####################
def writeLAMMPS_GO_wc(fileNameNew, arrAtom, arrNonbonded) :#{{{
    # Inputs:
    # fileNameNew           This must be the file name that will be created that contains the native contact data
    # arrAtom               This is an arrays of class variables. Make sure that you have all needed
    #                         data for the write

    # Actons:
    # Will produce a list of LAMMPS wall_coeff commands that enable wall hydrophobicity interaction

    # Outputs:
    # Nothing is returned, but the fileNameNew file is generated.

    # Set-up and Pre-define variables
    totAA = len(arrAtom)

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeData = open(fileNameNew, 'w')

    # Wrtie pair_coeff NBfix modifications
    print("Writing 'HISAA wall coefficient mod file' section...")
    for i in range(totAA) :
        a = str("%i" % (arrAtom[i].atomNum)).ljust(4)
        b = str("%.6f" % (arrNonbonded[i].epsil)).rjust(9)
        line = ("set             atom  %s  d_wc_epsilon  %s\n" % (a, b))
        writeData.write(line)
        b = str("%.6f" % (arrNonbonded[i].rmin2 * KPI.nbRminToEten)).rjust(9)
        line = ("set             atom  %s  d_wc_sigma    %s\n" % (a, b))
        writeData.write(line)

    writeData.close()
    print("Finished writing %s" % (fileNameNew))#}}}
#########################################################


################### writeKPREP_kp #######################
#writeKPREP_kp #{{{
def writeKPREP_kp(fileNameNew, arrAtom=[], arrBond=[], arrAngle=[], arrDihedral=[], \
                  arrImproper=[], arrNonbonded=[]) :
    # Inputs

    # Actions

    # Outputs

    # Open file for writing
    print("Writing file: %s" % fileNameNew)
    writeFile = open(fileNameNew, 'w')

    # Write atom data
    line = "KPI.Atom\n%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % \
            ("type".ljust(8), "atomNum".ljust(8), "atomSeq".ljust(8), "atomName".ljust(9), \
             "atomNameAlt".ljust(12), "resID".ljust(6), "resType".ljust(8), "chainType".ljust(10), \
             "resNum".ljust(8), "resSeq".ljust(8), "coordinates".ljust(29), "o".ljust(6), \
             "tf".ljust(6), "element".ljust(8), "q".ljust(6), "mw".ljust(8), \
             "atomID".ljust(8), "atomType".ljust(8))
    writeFile.write(line)
    for atom in arrAtom :
        a = atom.type.ljust(8)
        if a.replace(' ','') == '' : a = "N/A".ljust(8)
        b = str(atom.atomNum).ljust(8)
        if b.replace(' ','') == '' : b = "N/A".ljust(8)
        c = str(atom.atomSeq).ljust(8)
        if c.replace(' ','') == '' : c = "N/A".ljust(8)
        d = atom.atomName.ljust(9)
        if d.replace(' ','') == '' : d = "N/A".ljust(9)
        e = atom.atomNameAlt.ljust(12)
        if e.replace(' ','') == '' : e = "N/A".ljust(12)
        f = str(atom.resID).ljust(6)
        if f.replace(' ','') == '' : f = "N/A".ljust(6)
        g = atom.resType.ljust(8)
        if g.replace(' ','') == '' : g = "N/A".ljust(8)
        h = atom.chainType.ljust(10)
        if h.replace(' ','') == '' : h = "N/A".ljust(10)
        i = str(atom.resNum).ljust(8)
        if i.replace(' ','') == '' : i = "N/A".ljust(8)
        j = str(atom.resSeq).ljust(8)
        if j.replace(' ','') == '' : j = "N/A".ljust(8)
        k = str(atom.coordinates[0]).ljust(9)
        if k.replace(' ','') == '' : k = "N/A".ljust(9)
        l = str(atom.coordinates[1]).ljust(9)
        if l.replace(' ','') == '' : l = "N/A".ljust(9)
        m = str(atom.coordinates[2]).ljust(9)
        if m.replace(' ','') == '' : m = "N/A".ljust(9)
        n = str(atom.o).ljust(6)
        if n.replace(' ','') == '' : n = "N/A".ljust(6)
        o = str(atom.tf).ljust(6)
        if o.replace(' ','') == '' : o = "N/A".ljust(6)
        p = atom.element.ljust(8)
        if p.replace(' ','') == '' : p = "N/A".ljust(8)
        q = str(atom.q).ljust(6)
        if q.replace(' ','') == '' : q = "N/A".ljust(6)
        r = str(atom.mw).ljust(8)
        if r.replace(' ','') == '' : r = "N/A".ljust(8)
        s = str(atom.atomID).ljust(8)
        if s.replace(' ','') == '' : s = "N/A".ljust(8)
        t = atom.atomType.ljust(8)
        if t.replace(' ','') == '' : t = "N/A".ljust(8)
        line = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % \
                (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write bond data
    line = "\nKPI.Bond\n%s %s %s %s %s %s\n" % \
            ("bondID".ljust(8), "bondNum".ljust(8), "bond".ljust(19), \
             "bondType".ljust(19), "sigma".ljust(8), "epsil".ljust(8))
    writeFile.write(line)
    for bond in arrBond :
        a = str(bond.bondID).ljust(8)
        if a.replace(' ','') == '' : a = "N/A".ljust(8)
        b = str(bond.bondNum).ljust(8)
        if b.replace(' ','') == '' : b = "N/A".ljust(8)
        c = str(bond.bond[0]).ljust(9)
        if c.replace(' ','') == '' : c = "N/A".ljust(9)
        d = str(bond.bond[1]).ljust(9)
        if d.replace(' ','') == '' : d = "N/A".ljust(9)
        e = bond.bondType[0].ljust(9)
        if e.replace(' ','') == '' : e = "N/A".ljust(9)
        f = bond.bondType[1].ljust(9)
        if f.replace(' ','') == '' : f = "N/A".ljust(9)
        g = str(bond.sigma).ljust(8)
        if g.replace(' ','') == '' : g = "N/A".ljust(8)
        h = str(bond.epsil).ljust(8)
        if h.replace(' ','') == '' : g = "N/A".ljust(8)
        line = "%s %s %s %s %s %s %s %s\n" % \
                (a,b,c,d,e,f,g,h)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write angle data
    line = "\nKPI.Angle\n%s %s %s %s %s %s\n" % \
            ("angleID".ljust(8), "angleNum".ljust(8), "angle".ljust(29), \
             "angleType".ljust(29), "sigma".ljust(8), "epsil".ljust(8))
    writeFile.write(line)
    for angle in arrAngle :
        a = str(angle.angleID).ljust(8)
        if a.replace(' ','') == '' : a = "N/A".ljust(8)
        b = str(angle.angleNum).ljust(8)
        if b.replace(' ','') == '' : b = "N/A".ljust(8)
        c = str(angle.angle[0]).ljust(9)
        if c.replace(' ','') == '' : c = "N/A".ljust(9)
        d = str(angle.angle[1]).ljust(9)
        if d.replace(' ','') == '' : d = "N/A".ljust(9)
        e = str(angle.angle[2]).ljust(9)
        if e.replace(' ','') == '' : e = "N/A".ljust(9)
        f = angle.angleType[0].ljust(9)
        if f.replace(' ','') == '' : f = "N/A".ljust(9)
        g = angle.angleType[1].ljust(9)
        if g.replace(' ','') == '' : g = "N/A".ljust(9)
        h = angle.angleType[2].ljust(9)
        if h.replace(' ','') == '' : h = "N/A".ljust(9)
        i = str(angle.sigma).ljust(8)
        if i.replace(' ','') == '' : i = "N/A".ljust(8)
        j = str(angle.epsil).ljust(8)
        if j.replace(' ','') == '' : j = "N/A".ljust(8)
        line = "%s %s %s %s %s %s %s %s %s %s\n" % \
                (a,b,c,d,e,f,g,h,i,j)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write dihedral data
    line = "\nKPI.Dihedral\n%s %s %s %s %s %s %s %s\n" % \
            ("dihedralID".ljust(10), "dihedralNum".ljust(11), "dihedral".ljust(39), \
             "dihedralType".ljust(39), "K".ljust(8), "n".ljust(8), "d".ljust(8), "wt".ljust(8))
    writeFile.write(line)
    for dihedral in arrDihedral :
        a = str(dihedral.dihedralID).ljust(10)
        if a.replace(' ','') == '' : a = "N/A".ljust(10)
        b = str(dihedral.dihedralNum).ljust(11)
        if b.replace(' ','') == '' : b = "N/A".ljust(11)
        c = str(dihedral.dihedral[0]).ljust(9)
        if c.replace(' ','') == '' : c = "N/A".ljust(9)
        d = str(dihedral.dihedral[1]).ljust(9)
        if d.replace(' ','') == '' : d = "N/A".ljust(9)
        e = str(dihedral.dihedral[2]).ljust(9)
        if e.replace(' ','') == '' : e = "N/A".ljust(9)
        f = str(dihedral.dihedral[3]).ljust(9)
        if f.replace(' ','') == '' : f = "N/A".ljust(9)
        g = dihedral.dihedralType[0].ljust(9)
        if g.replace(' ','') == '' : g = "N/A".ljust(9)
        h = dihedral.dihedralType[1].ljust(9)
        if h.replace(' ','') == '' : h = "N/A".ljust(9)
        i = dihedral.dihedralType[2].ljust(9)
        if i.replace(' ','') == '' : i = "N/A".ljust(9)
        j = dihedral.dihedralType[3].ljust(9)
        if j.replace(' ','') == '' : j = "N/A".ljust(9)
        k = str(dihedral.K).ljust(8)
        if k.replace(' ','') == '' : k = "N/A".ljust(8)
        l = str(dihedral.n).ljust(8)
        if l.replace(' ','') == '' : l = "N/A".ljust(8)
        m = str(dihedral.d).ljust(8)
        if m.replace(' ','') == '' : m = "N/A".ljust(8)
        n = str(dihedral.wt).ljust(8)
        if n.replace(' ','') == '' : n = "N/A".ljust(8)
        line = "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % \
                (a,b,c,d,e,f,g,h,i,j,k,l,m,n)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write nonbonded data
    line = "\nKPI.Nonbonded\n%s %s %s %s %s\n" % \
            ("atomID".ljust(8), "atomNum".ljust(8), \
             "atomType".ljust(8), "sigma".ljust(8), "epsil".ljust(8))
    writeFile.write(line)
    for nonbonded in arrNonbonded :
        a = str(nonbonded.atomID).ljust(8)
        if a.replace(' ','') == '' : a = "N/A".ljust(8)
        b = str(nonbonded.atomNum).ljust(8)
        if b.replace(' ','') == '' : b = "N/A".ljust(8)
        c = nonbonded.atomType.ljust(8)
        if c.replace(' ','') == '' : c = "N/A".ljust(8)
        d = str(nonbonded.sigma).ljust(8)
        if d.replace(' ','') == '' : d = "N/A".ljust(8)
        e = str(nonbonded.epsil).ljust(8)
        if e.replace(' ','') == '' : e = "N/A".ljust(8)
        line = "%s %s %s %s %s\n" % \
                (a,b,c,d,e)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write improper data
    line = "\nKPI.Improper\n%s %s %s %s %s %s\n" % \
            ("improperID".ljust(10), "improperNum".ljust(11), "improper".ljust(39), \
             "improperType".ljust(39), "K".ljust(8), "chi".ljust(8))
    writeFile.write(line)
    for improper in arrImproper :
        a = str(improper.improperID).ljust(10)
        if a.replace(' ','') == '' : a = "N/A".ljust(10)
        b = str(improper.improperNum).ljust(11)
        if b.replace(' ','') == '' : b = "N/A".ljust(11)
        c = str(improper.improper[0]).ljust(9)
        if c.replace(' ','') == '' : c = "N/A".ljust(9)
        d = str(improper.improper[1]).ljust(9)
        if d.replace(' ','') == '' : d = "N/A".ljust(9)
        e = str(improper.improper[2]).ljust(9)
        if e.replace(' ','') == '' : e = "N/A".ljust(9)
        f = str(improper.improper[3]).ljust(9)
        if f.replace(' ','') == '' : f = "N/A".ljust(9)
        g = improper.improperType[0].ljust(9)
        if g.replace(' ','') == '' : g = "N/A".ljust(9)
        h = improper.improperType[1].ljust(9)
        if h.replace(' ','') == '' : h = "N/A".ljust(9)
        i = improper.improperType[2].ljust(9)
        if i.replace(' ','') == '' : i = "N/A".ljust(9)
        j = improper.improperType[3].ljust(9)
        if j.replace(' ','') == '' : j = "N/A".ljust(9)
        k = str(improper.K).ljust(8)
        if k.replace(' ','') == '' : k = "N/A".ljust(8)
        l = str(improper.chi).ljust(8)
        if l.replace(' ','') == '' : l = "N/A".ljust(8)
        line = "%s %s %s %s %s %s %s %s %s %s %s %s\n" % \
                (a,b,c,d,e,f,g,h,i,j,k,l)
        writeFile.write(line)
    writeFile.write("END\n")

    #Write UreyBradley

    #Write NBfix
    
    writeFile.close()
#}}}
#########################################################


