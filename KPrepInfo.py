#
# KPrepInfo.py
# --------------------
# Authors: Addison Smith and Derek Bush
# Contributors: Anthony Gillepsie
# Written to be run on python 2 or 3

#########################################################
#--------------- Index of Classes -----------------------
# Atom          Contains information about each model atom; Smallest reacting site
# Bond          Contains bond information about the molecule
# Angle         Contains angle information about the molecule
# Dihedral      Contains dihedral information about the modlecule
# Nonbonded     Contains CHARMM-NONBONDED parameters; compair to LAMMPS "Pair Coeff" parameters 
# Improper      Contains Improper parameters
# UreyBradley   Contains Urey-Bradley parameters
# NBfix         Contains CHARMM-NBFIX PARAMETERS; Used in Go-model native contact Pair-like parameters
#--------------------------------------------------------
# Note: to access class data: ex. access Atom resID value: UserAtomVarName.resID
#--------------------------------------------------------
# Note: once initialized any of the class values can be changed/redefined (see above)
#-------------------------------------------------------
# Note: these structures use "Num" and "ID" identifiers:
# -ID typically means this is a possible repeating type identifier, for example all CHARMM-type "PHE" 
#   atoms would have the same resID
# -Num means the sequencial and unique class entry in order of appearance
# You may have multiple "ID"s but you should only code one unique "Num" per atom per system/molecule
#-------------------------------------------------------
# Note: "Num" identifiers are a 1-base identification
#--------------------------------------------------------
# Note: for a detailed description on how the classes were set-up, see bottom of page
#--------------------------------------------------------
# Note: to toggle VIM folds: put cursor at your location -> command mode -> "za"
#       to create a fold: command mode -> ex. ":12,34fold" (the numbers are line numbers)
#--------------------------------------------------------
#########################################################
#--------------- Index of Functions ---------------------
# showAtom()
# showAtomInit()
# showBond()
# showBondInit()
# showAngle()
# showAngleInit()
# showDihedral()
# showDihedralInit()
# showNonbonded()
# showNonbondedInit()
# showImproper()
# showImproperInit()
# showUreyBradley()
# showUreyBradleyInit()
# showNBfix()
# showNBfixInit()
## showHISAAInit()           # (TBD)
#--------------------------------------------------------
# Note: All previous "show***()" functions print to screen the Atom variable's class values
# Note: All previous "show***Init()" functions print to screen thet class' variables and default values
#########################################################


####################### Constants #######################
# Constants    #{{{
PI = 3.1415926535898
nbRminToSigma = 2**(5/6) # Conversion between CHARMM Rmin/2 to LAMMPS LJ 12-6 sigma for pair_coeff calculations
                         # sigma = (Rmin/2) * nbRmintoSigma; where Rmin/2 = CHARMM value reported in NONBONDED section
nbRminToEten = 2         # Conversion between CHARMM Rmin/2 to LAMMPS eten 12-10-6 sigma for pair_coeff calculations
                         # sigma = (Rmin/2) * nbRmintoEten; where Rmin/2 = CHARMM value reported in NONBONDED section
                         # KPrepInfo stores Rmin/2 values in nonbonded.rmin2     #}}}
tsConv = 0.04888821      # Conversion between CHARMM time units in param file to ps. Used in timestep conversions
#########################################################


############## Dictionaries and Arrays ##################
#{{{
dictAA = {"ALA" : 1,  #Alanine          (A)
          "ARG" : 18, #Arginine         (R)
          "ASN" : 14, #Asparagine       (N)
          "ASP" : 4,  #Aspartic Acid    (D)
          "CYS" : 3,  #Cysteine         (C)
          "GLU" : 5,  #Glutamic Acid    (E)
          "GLN" : 17, #Glutamine        (Q)
          "GLY" : 7,  #Glycine          (G)
          "HIS" : 8,  #Histidine        (H)
          "HSD" : 8,  #Histidine        (H)
          "ILE" : 9,  #Isoleucine       (I)
          "LEU" : 12, #Leucine          (L)
          "LYS" : 11, #Lysine           (K)
          "MET" : 13, #Methionine       (M)
          "PHE" : 6,  #Phenylalanine    (F)
          "PRO" : 16, #Proline          (P)
          "SER" : 19, #Serine           (S)
          "THR" : 20, #Threonine        (T)
          "TRP" : 23, #Typtophan        (W)
          "TYR" : 25, #Tyrosine         (Y)
          "VAL" : 22  #Valine           (V)
          }
dictAAname = {"ALA" : "A",  #Alanine          
              "ARG" : "R",  #Arginine         
              "ASN" : "N",  #Asparagine       
              "ASP" : "D",  #Aspartic Acid    
              "CYS" : "C",  #Cysteine         
              "GLU" : "E",  #Glutamic Acid    
              "GLN" : "Q",  #Glutamine        
              "GLY" : "G",  #Glycine          
              "HIS" : "H",  #Histidine        
              "ILE" : "I",  #Isoleucine       
              "LEU" : "L",  #Leucine          
              "LYS" : "K",  #Lysine           
              "MET" : "M",  #Methionine       
              "PHE" : "F",  #Phenylalanine    
              "PRO" : "P",  #Proline          
              "SER" : "S",  #Serine           
              "THR" : "T",  #Threonine        
              "TRP" : "W",  #Typtophan        
              "TYR" : "Y",  #Tyrosine         
              "VAL" : "V"   #Valine           
              }
alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', \
            'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', \
            'U', 'V', 'W', 'X', 'Y', 'Z']

#}}}
#########################################################


####################### Atom ############################
class Atom:  #{{{        
                     # PDB col info
    type = ''        # cols  1 - 6  # CHARMM style atom type (ATOM, HETATOM, etc.) 
    atomNum = 0                     # Atom number as indicated by order in PDB (always starts at 1)
    atomSeq = 0      # cols  7 - 11 # Atom serial number as sequenced by the PDB structure file (otherwise equals atomNum)
    atomName = ''    # cols 13 - 16 # Atom name as indicated by the PDB (typcially a CHARMM atom like CA or HB1)
    atomNameAlt = '' # cols 17      # Alternate location identifier
    resID = 0                       # The type of residue, all resIDs have the same resType
    resType = ''     # cols 18 - 20 # Atom residue name (typically a CHARMM residue name like TYR)
    chainType = ''   # cols 22      # The chain identifier associated with this atom
    resNum = 0                      # Residue number as indicated by order in PDB (always starts at 1)
    resSeq = 0       # cols 23-26   # Residue serial numer as sequenced by the PDB structure file
    coordinates = [0.0, 0.0, 0.0]   # x, y, z
    o = 0.0          # cols 55 - 60 # Occupancy
    tf = 0.0         # cols 61 - 66 # Temperature factor
    element = ''     # cols 77 - 78 # Element symbol, right justified
    q = 0.0          # cols 79 - 80 # Atom charge
    mw = 0.0                        # Molecular weight of the atom. This is not defined in the PDB
    atomID = 0                      # The type of atom, all atomIDs have the same atomType
    atomType = ''                   # The string code that identifies that atom type in .rtf and .prm files

    # Potential growth areas that could be added:
    #angles                         # saves all the angleIDs it is a part of
    #bonds, dihedrals, nbfixes      # same as angle, saves IDs this atom is a part of

    def __init__(self, type='', atomNum=0, atomSeq=0, atomName='', atomNameAlt='', resID=0, resType='', chainType='', \
                       resNum=0, resSeq=0, x=0.0, y=0.0, z=0.0, o=0.0, tf=0.0, element='', q=0.0, mw=0.0, \
                       atomID=0, atomType='') :
        self.type = type
        self.atomNum = atomNum
        self.atomSeq = atomSeq
        self.atomName = atomName
        self.atomNameAlt = atomNameAlt
        self.resID = resID
        self.resType = resType
        self.chainType = chainType
        self.resNum = resNum
        self.resSeq = resSeq
        self.coordinates = [x,y,z]
        self.o = o
        self.tf = tf
        self.element = element
        self.q = q
        self.mw = mw
        self.atomID = atomID
        self.atomType = atomType #}}}

def showAtom(atomName, varName="Atom") : #{{{
    # This function prints to screen all the class values in the passed-in Atom variable.
    # If defined, varName is supposed to be the name of the Atom class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'type = \'%s\'        \n' % atomName.type
    str += 'atomNum = %i         \n' % atomName.atomNum
    str += 'atomSeq = %i         \n' % atomName.atomSeq
    str += 'atomName = \'%s\'    \n' % atomName.atomName
    str += 'atomNameAlt = \'%s\' \n' % atomName.atomNameAlt
    str += 'resID = %i           \n' % atomName.resID
    str += 'resType = \'%s\'     \n' % atomName.resType
    str += 'chainType = \'%s\'   \n' % atomName.chainType
    str += 'resNum = %i          \n' % atomName.resNum
    str += 'resSeq = %i          \n' % atomName.resSeq
    str += 'coordinates = [%f, %f, %f]  \n' % (atomName.coordinates[0], atomName.coordinates[1], atomName.coordinates[2])
    str += 'o = %f               \n' % atomName.o
    str += 'tf = %f              \n' % atomName.tf
    str += 'element = \'%s\'     \n' % atomName.element
    str += 'q = %f               \n' % atomName.q
    str += 'mw = %f              \n' % atomName.mw
    str += 'atomID = %i          \n' % atomName.atomID
    str += 'atomType = \'%s\'    \n' % atomName.atomType
    print(str)#}}}

def showAtomInit() : #{{{
    # This function prints to screen all the variables in this class and their descriptions
    str  = '\n'
    str += 'Atom class variables incude:\n'
    str += 'type = \'\'                       - CHARMM style atom type (ATOM, HETATOM, etc.) \n'
    str += 'atomNum = 0                     - Atom number as indicated by order in PDB (always starts at 1) \n'
    str += 'atomSeq = 0                     - Atom serial number as sequenced by the PDB structure file \n'
    str += 'atomName = \'\'                   - Atom name as indicated by the PDB (typcially a CHARMM atom like CA or HB1) \n'
    str += 'atomNameAlt = \'\'                - Alternate location identifier \n'
    str += 'resID = 0                       - The type of residue, all resIDs have the same resType \n'
    str += 'resType = \'\'                    - Atom residue name (typically a CHARMM residue name like TYR) \n'
    str += 'chainType = \'\'                  - The chain identifier associated with this atom \n'
    str += 'resNum = 0                      - Residue number as indicated by order in PDB (always starts at 1) \n'
    str += 'resSeq = 0                      - Residue serial numer as sequenced by the PDB structure file \n'
    str += 'coordinates = [0.0, 0.0, 0.0]   - Atom x, y, z coordinates \n'
    str += 'o = 0.0                         - Occupancy \n'
    str += 'tf = 0.0                        - Temperature factor \n'
    str += 'element = \'\'                    - Element symbol, right justified \n'
    str += 'q = 0.0                         - Atom charge \n'
    str += 'mw = 0.0                        - Molecular weight of the atom. This is not defined in the PDB \n'
    str += 'atomID = 0                      - The type of atom, all atomIDs have the same atomType \n'
    str += 'atomType = \'\'                   - The string code that identifies that atom type in .rtf and .prm files \n'
    str += '\nTo create a new Atom variable: \n'
    str += 'UserAtomName = KPrepInfo.Atom(type="ATOM", x=1.23, atomID=45, ...)\n'
    str += '\nTo edit an Atom: \n'
    str += 'UserAtomName.atomNum = 1234 \n'
    print(str)#}}}
#########################################################


######################## Bond ###########################
class Bond: #{{{
    bondID = 0                      # The type of bond, all bondIDs have the same bondType
    bondNum = 0                     # Sequencial identifier of the bond
    bond = [0, 0]                   # Containts the atomIDs of the two atoms involved
    bondType = ['', '']             # Containts the atomTypes of the two atoms involved
    sigma = 0.0                     # The LJ bond distance
    epsil = 0.0                     # The LJ bond energy
    
    def __init__(self, bondID=0, bondNum=0, atom1=0, atom2=0, atomType1='', atomType2='', sigma=0.0, epsil=0.0) :
        self.bondID = bondID
        self.bondNum = bondNum
        self.bond = [atom1, atom2]
        self.bondType = [atomType1, atomType2]
        self.sigma = sigma
        self.epsil = epsil #}}}

def showBond(bondName, varName="Bond") : #{{{
    # This function prints to screen all the class values in the passed-in Bond variable.
    # If defined, varName is supposed to be the name of the Bond class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'bondID = %i       \n' % bondName.bondID
    str += 'bondNum = %i      \n' % bondName.bondNum
    str += 'bond = [%f, %f]   \n' % (bondName.bond[0], bondName.bond[1])
    str += 'bondType = [\'%s\', \'%s\']  \n' % (bondName.bondType[0], bondName.bondType[1])
    str += 'sigma = %f        \n' % bondName.sigma
    str += 'epsil = %f        \n' % bondName.epsil
    print(str) #}}}

def showBondInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'Bond class variables include:\n'
    str += 'bondID = 0                      - The type of bond, all bondIDs have the same bondType\n'
    str += 'bondNum = 0                     - Sequencial identifier of the bond\n'
    str += 'bond = [0, 0]                   - Containts the atomIDs of the two atoms involved\n'
    str += 'bondType = [\'\', \'\']             - Containts the atomTypes of the two atoms involved\n'
    str += 'sigma = 0.0                     - The LJ bond distance\n'
    str += 'epsil = 0.0                     - The LJ bond energy\n'
    str += '\nTo create a new bond variable:\n'
    str += 'UserBondName = KPrepInfo.Bond(bondNum=71, atom1=5, atom2=6, ...)\n'
    str += '\nTo edit a Bond: \n'
    str += 'UserBondName.sigma = 3.33 \n'
    print(str) #}}}
#########################################################


######################## Angle ##########################
class Angle: #{{{
    angleID = 0                     # The type of angle, all angleIDs have the same angleType
    angleNum = 0                    # Sequencial identifier of the angle
    angle = [0, 0, 0]               # Contains the atomIDs of the three atoms involved
    angleType = ['', '', '']        # Contains the atomTypes of the three atoms involved
    sigma = 0.0                     # The LJ angle distance
    epsil = 0.0                     # The LJ angle energy

    def __init__(self, angleID=0, angleNum=0, atom1=0, atom2=0, atom3=0, atomType1='', \
                 atomType2='', atomType3='', sigma=0.0, epsil=0.0) :
        self.angleID = angleID
        self.angleNum = angleNum
        self.angle = [atom1, atom2, atom3]
        self.angleType = [atomType1, atomType2, atomType3]
        self.sigma = sigma
        self.epsil = epsil #}}}

def showAngle(angleName, varName="Angle") : #{{{
    # This function prints to screen all the class values in the passed-in Angle variable.
    # If defined, varName is supposed to be the name of the Angle class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'angleID = %i         \n' % angleName.angleID
    str += 'angleNum = %i        \n' % angleName.angleNum
    str += 'angle = [%i, %i, %i]   \n' % (angleName.angle[0], angleName.angle[1], angleName.angle[2])
    str += 'angleType = [\'%s\', \'%s\', \'%s\'] \n' % \
            (angleName.angleType[0], angleName.angleType[1], angleName.angleType[2])
    str += 'sigma = %f           \n' % angleName.sigma
    str += 'epsil = %f           \n' % angleName.epsil
    print(str) #}}}

def showAngleInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'Angle class variables include:\n'
    str += 'angleID = 0                     - The type of angle, all angleIDs have the same angleType\n'
    str += 'angleNum = 0                    - Sequencial identifier of the angle\n'
    str += 'angle = [0, 0, 0]               - Contains the atomIDs of the three atoms involved\n'
    str += 'angleType = [\'\', \'\', \'\']        - Contains the atomTypes of the three atoms involved\n'
    str += 'sigma = 0.0                     - The LJ angle distance\n'
    str += 'epsil = 0.0                     - The LJ angle energy\n'
    str += '\nTo create a new angle variable:\n'
    str += 'UserAngleName = KPrepInfo.Angle(angleNum=71, atom1=5, atom2=6, ...)\n'
    str += '\nTo edit an Angle: \n'
    str += 'UserAngleName.angleType[2] = \'HA1\' \n'
    print(str) #}}}
#########################################################


###################### Dihedral #########################
class Dihedral: #{{{
    dihedralID = 0                  # The type of dihedralral, all dihedralIDs have the same dihedralType
    dihedralNum = 0                 # Sequencial identifier of the dihedral
    dihedral = [0, 0, 0, 0]         # Contains the atomIDs of the four atoms involved
    dihedralType = ['', '', '','']  # Contains the atomTypes of the four atoms involved
    K = 0.0                         # Energy of the interaction
    n = 0                           # Phase
    d = 0                           # Dihedral phase angle (0 or 180 degrees)
    wt = 0                          # Weighting factor to avoid double counting (0<wt<1)
    
    def __init__(self, dihedralNum=0, dihedralID=0, atom1=0, atom2=0, atom3=0, atom4=0, \
                 atomType1='', atomType2='', atomType3='', atomType4='', K=0.0, n=0, d=0, wt=0) :
        self.dihedralNum = dihedralNum
        self.dihedralID = dihedralID
        self.dihedral = [atom1, atom2, atom3, atom4]
        self.dihedralType = [atomType1, atomType2, atomType3, atomType4]
        self.K = K
        self.n = n
        self.d = d
        self.wt = wt #}}}

def showDihedral(dihedName, varName="Dihedral") : #{{{
    # This function prints to screen all the class values in the passed-in Dihedral variable.
    # If defined, varName is supposed to be the name of the Dihedral class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'dihedralID = %i   \n' % dihedName.dihedralID
    str += 'dihedralNum = %i  \n' % dihedName.dihedralNum
    str += 'dihedral = [%i, %i, %i, %i] \n' % \
            (dihedName.dihedral[0], dihedName.dihedral[1], dihedName.dihedral[2], dihedName.dihedral[3])
    str += 'dihedralType = [\'%s\', \'%s\', \'%s\',\'%s\'] \n' % \
            (dihedName.dihedralType[0], dihedName.dihedralType[1], dihedName.dihedralType[2], dihedName.dihedralType[3])
    str += 'K = %f            \n' % dihedName.K
    str += 'n = %i            \n' % dihedName.n
    str += 'd = %i            \n' % dihedName.d
    str += 'wt = %i           \n' % dihedName.wt
    print(str) #}}}

def showDihedralInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'Dihedral class variables include:\n'
    str += 'dihedralID = 0                  - The type of dihedralral, all dihedralIDs have the same dihedralType\n'
    str += 'dihedralNum = 0                 - Sequencial identifier of the dihedral\n'
    str += 'dihedral = [0, 0, 0, 0]         - Contains the atomIDs of the four atoms involved\n'
    str += 'dihedralType = [\'\', \'\', \'\',\'\']  - Contains the atomTypes of the four atoms involved\n'
    str += 'K = 0.0                         - Energy of the interaction\n'
    str += 'n = 0                           - Phase\n'
    str += 'd = 0                           - Dihedral phase angle (0 or 180 degrees)\n'
    str += 'wt = 0                          - Weighting factor to avoid double counting (0<wt<1)\n'
    str += '\nTo create a new dihedral variable:\n'
    str += 'UserDihedName = KPrepInfo.Dihed(dihedralNum=45, dihedralID=789, atom1=67, ...)\n'
    str += '\nTo edit a Dihedral: \n'
    str += 'UserDihedName.K = 9.87 \n'
    print(str) #}}}
#########################################################


##################### Nonbonded #########################
class Nonbonded: #{{{
    atomID = 0                      # The type of atom, all atomIDs have the same atomType
    atomNum = 0                     # Atom number as indicated by order in PDB (always starts at 1)
    atomType = ''                   # The string code that identifies that atom type in .rtf and .prm files
    rmin2 = 0.0                     # The nonbonded sigma distance
    epsil = 0.0                     # The nonbonded epsilon energy
    
    def __init__(self, atomID=0, atomNum=0, atomType='', rmin2=0.0, epsil=0.0) :
        self.atomID = atomID
        self.atomNum = atomNum
        self.atomType = atomType
        self.rmin2 = rmin2
        self.epsil = epsil #}}}

def showNonbonded(nonbondedName, varName="Nonbonded") : #{{{
    # This function prints to screen all the class values in the passed-in Nonbonded variable.
    # If defined, varName is supposed to be the name of the Dihedral class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'atomID = %i        \n' % nonbondedName.atomID
    str += 'atomNum = %i       \n' % nonbondedName.atomNum
    str += 'atomType = \'%s\'  \n' % nonbondedName.atomType
    str += 'sigma = %f         \n' % nonbondedName.sigma
    str += 'epsil = %f         \n' % nonbondedName.epsil
    print(str)#}}}

def showNonbondedInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'Nonbonded class variables include:\n'
    str += 'atomID = 0                      - The type of atom, all atomIDs have the same atomType\n'
    str += 'atomNum = 0                     - Atom number as indicated by order in PDB (always starts at 1)\n'
    str += 'atomType = \'\'                   - The string code that identifies that atom type in .rtf and .prm files\n'
    str += 'sigma = 0.0                     - The nonbonded sigma distance\n'
    str += 'epsil = 0.0                     - The nonbonded epsilon energy\n'
    str += '\nTo create a new nonbonded:\n'
    str += 'UserNonbondedName = KPrepInfo.Nonbonded(atomID=65, nonbondedSigma=2.33, nonbondedEpsil=0.1)\n'
    str += '\nTo edit a Nonbonded: \n'
    str += 'UserNonbondedName.ID = 42 \n'
    print(str)#}}}
#########################################################


###################### Improper #########################
class Improper: #{{{
    improperID = 0                  # The type of improper
    improperNum = 0                 # Sequencial identifier of the improper
    improper = [0, 0, 0, 0]         # Contains the atomIDs of the four atoms involved
    improperType = ['', '', '', ''] # Contains the atomTypes of the four atoms involved
    K = 0.0                         # Energy of interaction
    chi = 0.0                       # Improper angle in degrees

    def __init__(self, improperID=0, improperNum=0, atom1=0, atom2=0, atom3=0, atom4=0, atomType1='', \
                 atomType2='', atomType3='', atomType4='', K=0.0, chi=0.0) :
        self.improperID = improperID
        self.improperNum = improperNum
        self.improper = [atom1, atom2, atom3, atom4]
        self.improperType = [atomType1, atomType2, atomType3, atomType4]
        self.K = K
        self.chi = chi #}}}

def showImproper(impName, varName="Improper") : #{{{
    # This function prints to screen all the class values in the passed-in Improper variable.
    # If defined, varName is supposed to be the name of the Improper class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'improperID = %i     \n' % impName.improperID
    str += 'improperNum = %i    \n' % impName.improperNum
    str += 'improper = [%i, %i, %i, %i]  \n' % \
            (impName.improper[0], impName.improper[1], impName.improper[2], impName.improper[3])
    str += 'improperType = [\'%s\', \'%s\', \'%s\', \'%s\'] \n' % \
            (impName.improperType[0], impName.improperType[1], impName.improperType[2], impName.improperType[3])
    str += 'K = %f              \n' % impName.K
    str += 'chi = %f            \n' % impName.chi
    print(str)#}}}

def showImproperInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'Improper class variables include:\n'
    str += 'improperID = 0                  - The type of improper\n'
    str += 'improperNum = 0                 - Sequencial identifier of the improper\n'
    str += 'improper = [0, 0, 0, 0]         - Contains the atomIDs of the four atoms involved\n'
    str += 'improperType = [\'\', \'\', \'\', \'\'] - Contains the atomTypes of the four atoms involved\n'
    str += 'K = 0.0                         - Energy of interaction\n'
    str += 'chi = 0.0                       - Improper angle in degrees\n'
    str += '\nTo create a new improper:\n'
    str += 'UserImpName = KPrepInfo.Improper(improperNum=4, improperID=89, ...)\n'
    str += '\nTo edit an Improper: \n'
    str += 'UserInpName.improper[0] = 99 \n'
    print(str)#}}}
#########################################################


##################### UreyBradley #######################
class UreyBradley: #{{{
    #Note: I don't know where/when we actually use this force field
    #      this means that the constants and descriptions could be wrong
    ubID = 0                        # The type of Urey-Bradley force
    ubNum = 0                       # Sequencial identifier of the Urty-Bradley force
    ub = [0, 0, 0]                  # Contains the atomIDs of the three atoms involved
    ubType = ['', '', '']           # Contains the atomTypes of the three atoms involved
    c1 = 0.0                        # From what I've seen in past uses there are
    c2 = 0.0                        # two constants
    
    def __init__(self, ubID=0, ubNum=0, atom1=0, atom2=0, atom3=0, atomType1='', \
                 atomType2='', atomType3='', c1=0.0, c2=0.0) :
        self.ubID = ubID
        self.ubNum = ubNum
        self.ub = [atom1, atom2, atom3]
        self.ubType = [atomType1, atomType2, atomType3]
        self.c1 = c1
        self.c2 = c2 #}}}

def showUreyBradley(ubName, varName="UreyBradley") : #{{{
    # This function prints to screen all the class values in the passed-in UreyBradley variable.
    # If defined, varName is supposed to be the name of the UreyBradley class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'ubID = %i         \n' % ubName.ubID
    str += 'ubNum = %i        \n' % ubName.ubNum
    str += 'ub = [%i, %i, %i] \n' % (ubName.ub[0], ubName.ub[1], ubName.ub[2])
    str += 'ubType = [\'%s\', \'%s\', \'%s\']  \n' % (ubName.ubType[0], ubName.ubType[1], ubName.ubType[2])
    str += 'c1 = %f           \n' % ubName.c1
    str += 'c2 = %f           \n' % ubName.c2
    print(str)#}}}

def showUreyBradleyInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'UreyBradley class variables include:\n'
    str += 'ubID = 0                        - The type of Urey-Bradley force\n'
    str += 'ubNum = 0                       - Sequencial identifier of the Urey-Bradley force\n'
    str += 'ub = [0, 0, 0]                  - Contains the atomIDs of the three atoms involved\n'
    str += 'ubType = [\'\', \'\', \'\']           - Contains the atomTypes of the three atoms involved\n'
    str += 'c1 = 0.0                        - From what I\'ve seen in past, there are\n'
    str += 'c2 = 0.0                        - two constants\n'
    str += '\nTo create a new UreyBradley force:\n'
    str += 'UserUbName = KPrepInfo.UreyBradley(ubNum=4, ubID=89, ...)\n'
    str += '\nTo edit a UreyBradley: \n'
    str += 'UserUbName.ubID = 777 \n'
    print(str)#}}}
#########################################################


######################## NBfix ##########################
class NBfix: #{{{
    nbfixID = 0                     # The nbfixID associated with this pair interaction
    nbfix = [0, 0]                  # Contains the atomIDs of the two atoms involved
    nbfixType = ['', '']            # Contains the atomTypes of the two atoms involved
    sigma = 0.0                     # The NBfix sigma distance
    epsil = 0.0                     # The NBfix epsilon energy

    def __init__(self, nbfixID=0, atom1=0, atom2=0, atomType1='', atomType2='', sigma=0.0, epsil=0.0) :
        self.nbfixID = nbfixID
        self.nbfix = [atom1, atom2]
        self.nbfixType = [atomType1, atomType2]
        self.sigma = sigma
        self.epsil = epsil #}}}

def showNBfix(nbfixName, varName="NBfix") : #{{{
    # This function prints to screen all the class values in the passed-in NBfix variable.
    # If defined, varName is supposed to be the name of the NBfix class variable as a string; useful for troubleshooting.
    str  = '\n'
    str += 'Assigned class values to %s variable:\n' % varName
    str += 'nbfixID = %i      \n' % nbfixName.nbfixID
    str += 'nbfix = [%i, %i]  \n' % (nbfixName.nbfix[0], nbfixName.nbfix[1])
    str += 'nbfixType = [\'%s\', \'%s\']  \n' % (nbfixName.nbfixType[0], nbfixName.nbfixType[1])
    str += 'sigma = %f        \n' % nbfixName.sigma
    str += 'epsil = %f        \n' % nbfixName.epsil
    print(str)#}}}

def showNBfixInit() : #{{{
    # This function prints to screen all the class variables in this class and the default values
    str  = '\n'
    str += 'NBfix class variables include:\n'
    str += 'nbfixID = 0                     - The nbfixID associated with this pair interaction\n'
    str += 'nbfix = [0, 0]                  - Contains the atomIDs of the two atoms involved\n'
    str += 'nbfixType = [\'\', \'\']            - Contains the atomTypes of the two atoms involved\n'
    str += 'sigma = 0.0                     - The NBfix sigma distance\n'
    str += 'epsil = 0.0                     - The NBfix epsilon energy\n'
    str += '\nTo create a new NBfix:\n'
    str += 'UserNBfixName = KPrepInfo.NBfix(nbfixID=65, atom1=8, atom2=56 sigma=2.33, ...)\n'
    str += '\nTo edit an NBfix: \n'
    str += 'UserNBfixName.epsil = 4.567 \n'
    print(str)#}}}
#########################################################


#--------- Description of Class set-up ------------------
############### Example: Bond class #####################
'''#{{{
class Bond: #This defines the class name 
            #To define/initialize a new "Bond" variable:   VarName = KPrepInfo.Bond()
            #You can also have arrays of "Bond" variables: VarName[i] = KPrepInfo.Bond()
            #.append works well if you are building your "Bond" array: VarName.append(KPrepInto.Bond())
    
    #This next field defines the class variables (and functions)
    #To reference such a variable: ex. VarName.bondID
    bondNum = 0  
    bondID = 0   
    .
    .
    .

    #All classes must first be initialized
    #We set-up __init__ so that it doesn't require any inputs
    def __init__(self, bondNum=0, bondID=0, ...) :
        self.bondNum = bondNum          #self.bondNum references the class var, 
                                        #bondNum is the function input
        self.bondID = bondID
        .
        .
        .

def showBondInit() :
    # This function prints to screen all the class values defined in this class variable.
    # If defined, varName is supposed to be the name of the NBfix class variable as a string; useful for troubleshooting.

def showBondInit() :
    # This function prints to screen all the class variables in this class and the default values
'''#}}}
#########################################################

