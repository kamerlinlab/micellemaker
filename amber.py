#!/usr/bin/python
# -*- coding: utf-8 -*- 

import os, sys, string, commands, shutil

AMBERHOME = "/home/helix/amber16/"
AMBEREXEC = "bin/pmemd.cuda"
os.environ['AMBERHOME'] = "/home/helix/amber16/"
os.environ['CUDA_HOME'] = "/usr/local/cuda-8.0/"
os.environ['PATH'] = "/usr/local/cuda-8.0/bin:$PATH"
os.environ['LD_LIBRARY_PATH'] = "/usr/local/cuda-8.0/lib64/:$LD_LIBRARY_PATH"
os.environ['LD_LIBRARY_PATH'] = "/home/helix/amber16/lib/:$LD_LIBRARY_PATH"

dict_lipid = {}

dict_lipid["HGA"] = "AMT C07"
dict_lipid["HGB"] = "BMT C07"
dict_lipid["OGA"] = "AMT C08"
dict_lipid["OGB"] = "BMT C08"
dict_lipid["NGA"] = "AMT C09"
dict_lipid["NGB"] = "BMT C09"
dict_lipid["DGA"] = "AMT C10"
dict_lipid["DGB"] = "BMT C10"
dict_lipid["UGA"] = "AMT C11"
dict_lipid["UGB"] = "BMT C11"
dict_lipid["EGA"] = "AMT C12"
dict_lipid["EGB"] = "BMT C12"

dict_lipid["HMA"] = "AMT AMG C07"
dict_lipid["HMB"] = "AMT BMG C07"
dict_lipid["OMA"] = "AMT AMG C08"
dict_lipid["OMB"] = "AMT BMG C08" 
dict_lipid["NMA"] = "AMT AMG C09"
dict_lipid["NMB"] = "AMT BMG C09"
dict_lipid["DMA"] = "AMT AMG C10"
dict_lipid["DMB"] = "AMT BMG C10"
dict_lipid["UMA"] = "AMT AMG C11"
dict_lipid["UMB"] = "AMT BMG C11"
dict_lipid["EMA"] = "AMT AMG C12"
dict_lipid["EMB"] = "AMT BMG C12"

dict_lipid["SDS"] = "SDS"

dict_lipid["UGN"] = "UGN"
dict_lipid["UGV"] = "UGV"
dict_lipid["UGL"] = "UGL"

dict_lipid["UAG"] = "UAG"
dict_lipid["UAA"] = "UAA"
dict_lipid["UAV"] = "UAV"
dict_lipid["UAL"] = "UAL"

dict_lipid["UVG"] = "UVG"
dict_lipid["UVA"] = "UVA"
dict_lipid["UVV"] = "UVV"
dict_lipid["UVL"] = "UVL"

dict_lipid["ULG"] = "ULG"
dict_lipid["ULA"] = "ULA"
dict_lipid["ULV"] = "ULV"
dict_lipid["ULL"] = "ULL"

dict_lipid["SUA"] = "SUA"
dict_lipid["SUV"] = "SUV"
dict_lipid["SUL"] = "SUL"

aAAS = ['UGN','UGV','UGL','UAG','UAA','UAV','UAL','UVG','UVA','UVV','UVL','ULG','ULA','ULV','ULL','SUA','SUV','SUL']

class do_amber():

  def __init__(self,sType,sMName,sIonP,sIonN,sConc,sBoxDim,iNo,sGPU_ID):
    self.sType = sType
    self.sMName = sMName
    self.sIonP = sIonP
    self.sIonN = sIonN
    self.sConc = sConc
    self.sBoxDim = sBoxDim
    self.iNo = iNo
    self.sGPU_ID = sGPU_ID
      
  def calc_ions(self,sBoxSize,sConc):
    # (0.5 Mol) * (6.02*10^23) / 1dm^3 / ((10^9)^3 * (Ang^3/dm^3)) -> 3*E-4 molecules/Ang^3 -> box = 77760 Ang^3 -> 23.328 molecules
    fMolAng = float(sConc) * (6.02*10**23) / ((10**9)**3) # molecules/Ang^3
    sNoIons = str(int(round(fMolAng * float(sBoxSize),0)))
    return sNoIons,sNoIons

  def check_quality(self):
    iC = int(commands.getstatusoutput("/usr/bin/grep NaN equil*.out | /usr/bin/wc -l")[1])
    iD = int(commands.getstatusoutput("/usr/bin/grep ERROR *.log | /usr/bin/wc -l")[1])
    iS = iC + iD
    if(iS > 0):
      print("  Error during equilibration. Restart.\n")
      return False
    else:
      return True

  def runLeap(self):
    commands.getstatusoutput(AMBERHOME+"bin/tleap -f leap_"+self.sMName+".ff14SB")

  def genLeapFile(self):  
  
    fXB = float(self.sBoxDim.split()[0])
    fYB = float(self.sBoxDim.split()[1])
    fZB = float(self.sBoxDim.split()[2])
    sBoxSize = str(fXB*fYB*fZB)
    sNoP,sNoN = self.calc_ions(sBoxSize,self.sConc)
    sPNam = self.sIonP
    
    if(self.sIonP == 'Mg2+' or self.sIonP == 'Ca2+'):
      sNoN = str(2*int(sNoN))
      self.sIonP = self.sIonP[:2].upper()
      
    print "  Add "+sNoP+" "+sPNam+" salt ions"
    print "  Add "+sNoN+" "+self.sIonN+" salt ions\n"
    
    if(self.sType == "SDS" or self.sType in aAAS):
      if(self.sIonP == 'MG' or self.sIonP == 'CA'):
        sNoP = str( int(sNoP) + int(round(float(self.iNo)/2.0,0)) )
        print "  Add "+str( int(round(float(self.iNo)/2.0,0)) )+" "+sPNam+" ions to neutralize the system\n"
      else:
        sNoP = str( int(sNoP)+int(self.iNo) )
        print "  Add "+str(self.iNo)+" "+sPNam+" ions to neutralize the system\n"
      
    sContent = """
logFile leap.log
#
# ----- leaprc for loading the ff14SB force field
# ----- NOTE: this is designed for PDB format 3!
#    Uses frcmod.ff14SB for proteins; ff99bsc0 for DNA; ff99bsc0_chiOL3 for RNA
#
#	load atom type hybridizations
#
addAtomTypes {
	{ "H"   "H" "sp3" }
	{ "HO"  "H" "sp3" }
	{ "HS"  "H" "sp3" }
	{ "H1"  "H" "sp3" }
	{ "H2"  "H" "sp3" }
	{ "H3"  "H" "sp3" }
	{ "H4"  "H" "sp3" }
	{ "H5"  "H" "sp3" }
	{ "HW"  "H" "sp3" }
	{ "HC"  "H" "sp3" }
	{ "HA"  "H" "sp3" }
	{ "HP"  "H" "sp3" }
	{ "HZ"  "H" "sp3" }
	{ "OH"  "O" "sp3" }
	{ "OS"  "O" "sp3" }
	{ "O"   "O" "sp2" }
	{ "O2"  "O" "sp2" }
	{ "OP"  "O" "sp2" }
	{ "OW"  "O" "sp3" }
	{ "CT"  "C" "sp3" }
	{ "CX"  "C" "sp3" }
	{ "C8"  "C" "sp3" }
	{ "2C"  "C" "sp3" }
	{ "3C"  "C" "sp3" }
	{ "CH"  "C" "sp3" }
	{ "CS"  "C" "sp2" }
	{ "C"   "C" "sp2" }
	{ "CO"   "C" "sp2" }
	{ "C*"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CB"  "C" "sp2" }
	{ "CC"  "C" "sp2" }
	{ "CN"  "C" "sp2" }
	{ "CM"  "C" "sp2" }
	{ "CK"  "C" "sp2" }
	{ "CQ"  "C" "sp2" }
	{ "CD"  "C" "sp2" }
	{ "C5"  "C" "sp2" }
	{ "C4"  "C" "sp2" }
	{ "CP"  "C" "sp2" }
	{ "CI"  "C" "sp3" }
	{ "CJ"  "C" "sp2" }
	{ "CW"  "C" "sp2" }
	{ "CV"  "C" "sp2" }
	{ "CR"  "C" "sp2" }
	{ "CA"  "C" "sp2" }
	{ "CY"  "C" "sp2" }
	{ "C0"  "Ca" "sp3" }
	{ "N"   "N" "sp2" }
	{ "NA"  "N" "sp2" }
	{ "N2"  "N" "sp2" }
	{ "N*"  "N" "sp2" }
	{ "NP"  "N" "sp2" }
	{ "NQ"  "N" "sp2" }
	{ "NB"  "N" "sp2" }
	{ "NC"  "N" "sp2" }
	{ "NT"  "N" "sp3" }
	{ "NY"  "N" "sp2" }
	{ "N3"  "N" "sp3" }
	{ "S"   "S" "sp3" }
	{ "SH"  "S" "sp3" }
	{ "P"   "P" "sp3" }
	{ "LP"  ""  "sp3" }
	{ "F"   "F" "sp3" }
	{ "Cl"  "Cl" "sp3" }
	{ "Br"  "Br" "sp3" }
	{ "I"   "I"  "sp3" }
	{ "F-"  "F" "sp3" }
	{ "Cl-"  "Cl" "sp3" }
	{ "Br-"  "Br" "sp3" }
	{ "I-"   "I"  "sp3" }
	{ "Li+"  "Li"  "sp3" }
	{ "Na+"  "Na"  "sp3" }
	{ "K+"  "K"  "sp3" }
	{ "Mg2+"  "Mg"  "sp3" }
        { "Zn2+"  "Zn"  "sp3" }
        { "Ca2+"  "Ca"  "sp3" }

# glycam
#	{ "OG"  "O" "sp3" }
#	{ "OL"  "O" "sp3" }
#	{ "AC"  "C" "sp3" }
#	{ "EC"  "C" "sp3" }
# 	{ "DA"  "Ca" "sp3" }
#        { "DB"  "Ca" "sp3" } 
#        { "DC"  "Ca" "sp3" } 
#        { "DD"  "Ca" "sp3" }
#        { "DE"  "Ca" "sp3" }
#        { "DF"  "Ca" "sp3" }
#        { "Ca" "Ca" "sp3" } 
}
#
#	Load the main parameter set.
#
parm10 = loadamberparams parm10.dat
frcmod14SB = loadamberparams frcmod.ff14SB
#
#	Load main chain and terminating amino acid libraries, nucleic acids
#
loadOff amino12.lib
loadOff aminoct12.lib
loadOff aminont12.lib
#
#       Load water and ions
# 
loadOff atomic_ions.lib
loadOff solvents.lib
HOH = TP3
WAT = TP3
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams frcmod.ions234lm_126_tip3p
#
#	Define the PDB name map for the amino acids and nucleic acids
#
addPdbResMap {
  { 0 "HYP" "NHYP" } { 1 "HYP" "CHYP" }
  { 0 "ALA" "NALA" } { 1 "ALA" "CALA" }
  { 0 "ARG" "NARG" } { 1 "ARG" "CARG" }
  { 0 "ASN" "NASN" } { 1 "ASN" "CASN" }
  { 0 "ASP" "NASP" } { 1 "ASP" "CASP" }
  { 0 "CYS" "NCYS" } { 1 "CYS" "CCYS" }
  { 0 "CYX" "NCYX" } { 1 "CYX" "CCYX" }
  { 0 "GLN" "NGLN" } { 1 "GLN" "CGLN" }
  { 0 "GLU" "NGLU" } { 1 "GLU" "CGLU" }
  { 0 "GLY" "NGLY" } { 1 "GLY" "CGLY" }
  { 0 "HID" "NHID" } { 1 "HID" "CHID" }
  { 0 "HIE" "NHIE" } { 1 "HIE" "CHIE" }
  { 0 "HIP" "NHIP" } { 1 "HIP" "CHIP" }
  { 0 "ILE" "NILE" } { 1 "ILE" "CILE" }
  { 0 "LEU" "NLEU" } { 1 "LEU" "CLEU" }
  { 0 "LYS" "NLYS" } { 1 "LYS" "CLYS" }
  { 0 "MET" "NMET" } { 1 "MET" "CMET" }
  { 0 "PHE" "NPHE" } { 1 "PHE" "CPHE" }
  { 0 "PRO" "NPRO" } { 1 "PRO" "CPRO" }
  { 0 "SER" "NSER" } { 1 "SER" "CSER" }
  { 0 "THR" "NTHR" } { 1 "THR" "CTHR" }
  { 0 "TRP" "NTRP" } { 1 "TRP" "CTRP" }
  { 0 "TYR" "NTYR" } { 1 "TYR" "CTYR" }
  { 0 "VAL" "NVAL" } { 1 "VAL" "CVAL" }
  { 0 "HIS" "NHIS" } { 1 "HIS" "CHIS" }
}

#  try to be good about reading in really old atom names as well:
addPdbAtomMap {
  { "O5*" "O5'" }
  { "C5*" "C5'" }
  { "C4*" "C4'" }
  { "O4*" "O4'" }
  { "C3*" "C3'" }
  { "O3*" "O3'" }
  { "C2*" "C2'" }
  { "O2*" "O2'" }
  { "C1*" "C1'" }
  { "C5M" "C7"  }
  { "H1*" "H1'" }
  { "H2*1" "H2'" }
  { "H2*2" "H2''" }
  { "H2'1" "H2'" }
  { "H2'2" "H2''" }
  { "H3*" "H3'" }
  { "H4*" "H4'" }
  { "H5*1" "H5'" }
  { "H5*2" "H5''" }
  { "H5'1" "H5'" }
  { "H5'2" "H5''" }
  { "HO'2" "HO2'" }
  { "H5T"  "HO5'" }
  { "H3T"  "HO3'" }
  { "O1'" "O4'" }
  { "OA"  "OP1" }
  { "OB"  "OP2" }
  { "O1P" "OP1" }
  { "O2P" "OP2" }
}

#
# assume that most often proteins use HIE
#
#NHIS = NHIE
#HIS = HIE
#CHIS = CHIE
#
# My stuff
#
#
# Bondi radii for igb = 5
#
set default PBradii mbondi2

loadamberparams GLYCAM_06j.dat
"""
    sContent2 = ""
    aOff = dict_lipid[self.sType].split()
    for i in range(len(aOff)):
       sContent2 += "loadOff ../amber/glycam_lib/"+aOff[i]+".off\n"
    sTag = aOff[len(aOff)-1]
    
    if(aOff[0] in aAAS):
      sContent2 += "source leaprc.gaff2\n"
      sContent2 += "loadamberparams gaff2.dat\n"
      sContent2 += "loadamberparams ../amber/glycam_lib/"+aOff[0]+".frcmod\n"
      
    sContent2 += "PDB = loadpdb "+self.sMName+".pdb\n"
    
    if(aOff[0] in aAAS):
      iZC = 1
      for k in range(self.iNo):
        sContent2 += "bond PDB."+str(iZC)+".C11 PDB."+str(iZC+1)+".N\n"
        iZC +=3
      commands.getstatusoutput("/bin/sed -i \'s/CLI/"+aOff[0]+"/g\' "+self.sMName+".pdb")
      sNoP = str(0)
      sNoN = str(0)
    
    if(aOff[0] in ["UVV","ULL","UAA"]):
      ff2 = open("tmpf.pdb","w")
      ff = open(self.sMName+".pdb","r")
      aRes = []
      iZ = 1
      iCo = 0
      for line in ff.readlines():
        if(line[:4] == "ATOM"):
          sAA = line[16:20]
          aRes.append(sAA)
          if(len(aRes) > 1):
            if(aRes[iCo] != aRes[iCo-1]):
              iZ +=1
          iCo+=1
          ff2.write(line[:16]+" "+line[17:25]+str(iZ)+line[27:])
        else:
          ff2.write(line)
      ff.close()
      ff2.close()
      shutil.copy("tmpf.pdb",self.sMName+".pdb")
          
    commands.getstatusoutput("/bin/sed -i 's/CLI/"+sTag+"/g' "+self.sMName+".pdb")
    
    sContent2 += "set PDB box {"+self.sBoxDim+"}\n"
    sContent2 += "solvatebox PDB TIP3PBOX 0.0\n"

#    sContent2 += "solvateoct PDB TIP3PBOX 11.0\n"
    sContent2 += "addions PDB "+self.sIonP+" "+sNoP+"\n"
    sContent2 += "addions PDB "+self.sIonN+" "+sNoN+"\n"

    sContent2 += "savepdb PDB "+self.sMName+"_solv_ions.pdb\n"
    sContent2 += "saveamberparm PDB "+self.sMName+"_solv_ions.prmtop "+self.sMName+"_solv_ions.inpcrd\n"
    sContent2 += "quit\n"
    
    FLeap = open("leap_"+self.sMName+".ff14SB","w")
    FLeap.write(sContent)
    FLeap.write(sContent2)
    FLeap.close()
    
    return self.sType

  def minimize(self,iMaxcyc):

    sCon = ""
    aOff = dict_lipid[self.sType].split()
    iOff = len(aOff)
    if(iOff < 3):
      sCon = "* * * "+aOff[0]+"\n* * * C12"
    else:
      sCon = "* * * "+aOff[0]+"\n* * * "+aOff[1]+"\n* * * C12"
    
    if(iOff == 1):
      if(aOff[0] == "UGN"):
        sCon = "* * * GLY\n* * * ALA\n* * * C10"
      if(aOff[0] == "UGV"):
        sCon = "* * * GLY\n* * * VAL\n* * * C10"
      if(aOff[0] == "UGL"):
        sCon = "* * * GLY\n* * * LEU\n* * * C10"
      
      if(aOff[0] == "UAG"):
        sCon = "* * * ALA\n* * * GLY\n* * * C10"
      if(aOff[0] == "UAA"):
        sCon = "* * * ALA\n* * * ALA\n* * * C10"
      if(aOff[0] == "UAV"):
        sCon = "* * * ALA\n* * * VAL\n* * * C10"
      if(aOff[0] == "UAL"):
        sCon = "* * * ALA\n* * * LEU\n* * * C10"
      
      if(aOff[0] == "UVG"):
        sCon = "* * * VAL\n* * * GLY\n* * * C10"
      if(aOff[0] == "UVA"):
        sCon = "* * * VAL\n* * * ALA\n* * * C10"
      if(aOff[0] == "UVV"):
        sCon = "* * * VAL\n* * * VAL\n* * * C10"
      if(aOff[0] == "UVL"):
        sCon = "* * * VAL\n* * * LEU\n* * * C10"
      
      if(aOff[0] == "ULG"):
        sCon = "* * * LEU\n* * * GLY\n* * * C10"
      if(aOff[0] == "ULA"):
        sCon = "* * * LEU\n* * * ALA\n* * * C10"
      if(aOff[0] == "ULV"):
        sCon = "* * * LEU\n* * * VAL\n* * * C10"
      if(aOff[0] == "ULL"):
        sCon = "* * * LEU\n* * * LEU\n* * * C10"
      
      if(aOff[0] == "SUA"):
        sCon = "* * * ALA\n* * * C10"
      if(aOff[0] == "SUV"):
        sCon = "* * * VAL\n* * * C10"
      if(aOff[0] == "SUL"):
        sCon = "* * * LEU\n* * * C10"
              
    sInput = """Initial repair minimization without position restraints
 &cntrl
  nmropt = 0,

  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 100,     ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 1,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 25,

  ipol   = 0,

  ibelly = 0,       ntr    = 0,

  imin   = 1,
  maxcyc = """+str(iMaxcyc)+""",
  ncyc   = 100000,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 10000,
  nscm   = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,
  ig     = 71277,
  ntt    = 0,
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 1,       tol    = 0.000001

 &end
END

""" 
    sFile = "min1_"+self.sMName
    FMin = open(sFile+".in","w")
    FMin.write(sInput)
    FMin.close()
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
#    commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i "+sFile+".in -o "+sFile+".out -p "+self.sMName+"_solv_ions.prmtop -c "+self.sMName+"_solv_ions.inpcrd -r "+sFile+".restrt -ref "+self.sMName+"_solv_ions.inpcrd")
    os.system(AMBERHOME+AMBEREXEC+" -O -i "+sFile+".in -o "+sFile+".out -p "+self.sMName+"_solv_ions.prmtop -c "+self.sMName+"_solv_ions.inpcrd -r "+sFile+".restrt -ref "+self.sMName+"_solv_ions.inpcrd")
    commands.getstatusoutput(AMBERHOME+"/bin/ambpdb -ctr -p "+self.sMName+"_solv_ions.prmtop < "+sFile+".restrt > "+sFile+".pdb")
    

    sInput = """Water minimization with position restraints
 &cntrl
  nmropt = 0,

  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 100,     ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 1,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 25,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 1,
  maxcyc = 5000,
  ncyc   = 2000,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 10000,
  nscm   = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,
  ig     = 71277,
  ntt    = 0,
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 1,       tol    = 0.000001
 &end
posres
25.0
FIND 
"""+sCon+"""
SEARCH
RES 1 10000
END  
END

""" 
    sFile = "min2_"+self.sMName
    FMin = open(sFile+".in","w")
    FMin.write(sInput)
    FMin.close()
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i "+sFile+".in -o "+sFile+".out -p "+self.sMName+"_solv_ions.prmtop -c min1_"+self.sMName+".restrt -r "+sFile+".restrt -ref min1_"+self.sMName+".restrt")
    commands.getstatusoutput(AMBERHOME+"/bin/ambpdb -ctr -p "+self.sMName+"_solv_ions.prmtop < "+sFile+".restrt > "+sFile+".pdb")
    

    sInput = """Water minimization with position restraints
 &cntrl
  nmropt = 0,

  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 100,     ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 1,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 25,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 1,
  maxcyc = 5000,
  ncyc   = 2000,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 10000,
  nscm   = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,
  ig     = 71277,
  ntt    = 0,
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 1,       tol    = 0.000001

 &end
posres
5.0
FIND 
"""+sCon+"""
SEARCH
RES 1 10000
END  
END

""" 
    sFile = "min3_"+self.sMName
    FMin = open(sFile+".in","w")
    FMin.write(sInput)
    FMin.close()
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i "+sFile+".in -o "+sFile+".out -p "+self.sMName+"_solv_ions.prmtop -c min2_"+self.sMName+".restrt -r "+sFile+".restrt -ref min2_"+self.sMName+".restrt")
    commands.getstatusoutput(AMBERHOME+"/bin/ambpdb -ctr -p "+self.sMName+"_solv_ions.prmtop < "+sFile+".restrt > "+sFile+".pdb")
    
    return self


  def equilibration(self):

    sCon = ""
    aOff = dict_lipid[self.sType].split()
    iOff = len(aOff)
    if(iOff < 3 and iOff > 1):
      sCon = "C1 * * "+aOff[0]+"\nC12 * * C12"
    if(iOff == 3):
      sCon = "C1 * * "+aOff[0]+"\nC1 * * "+aOff[1]+"\nC12 * * C12"
    if(iOff == 1):
      sCon = "C1 * * "+aOff[0]+"\nC12 * * SDS"  
    
    if(iOff == 1):
      if(aOff[0] in aAAS):
        sCon = "C1 * * "+aOff[0]+"\nC10 * * "+aOff[0]
    
    
    sInput ="""Minimization with position restraints
 &cntrl
  nmropt = 1, 

  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 100,     ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 1,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 25,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 1,       
  maxcyc = 5000,   
  ncyc   = 200,   
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 10000,
  nscm   = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,   
  ig     = 71277,
  ntt    = 0,
  tautp  = 0.5,

  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 1,       tol    = 0.000001,

  watnam='WAT', ! Water residues are named WAT
  owtnm='O',   ! Water oxygens are named O
 &end
 &wt
    type='END'
 &end
Membrane posres
5.0
FIND 
"""+sCon+"""
SEARCH
RES 1 10000
END  
END  
DISANG=restraints01.rest
LISTIN=POUT
LISTOUT=POUT
"""
    sFile = "equil0_"+self.sMName
    FEQ = open(sFile+".in","w")
    FEQ.write(sInput)
    FEQ.close()

    sInput ="""NVT MD with position restraints
 &cntrl
  nmropt = 1,

  ntx    = 1,       irest  = 0,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 5000,    ntwx   = 5000,    ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 2,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 25,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 0,
  maxcyc = 5000,
  ncyc   = 200,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 25000,
  nscm   = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 100.0,
  ig     = 71277,
  ntt    = 1,
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 2,       tol    = 0.000001,
  
  watnam='WAT', ! Water residues are named WAT
  owtnm='O',   ! Water oxygens are named O

 &end

 &wt
  type   = 'TEMP0', istep1 = 0,       istep2 = 20000,
                    value1 = 100.0,   value2 = 300.0,
 &end

 &wt
  type   = 'TEMP0', istep1 = 20001,   istep2 = 25000,
                    value1 = 300.0,   value2 = 300.0,
 &end
                                                      
 &wt
  type   = 'END'
 &end 
Membrane posres
25.0
FIND 
"""+sCon+"""
SEARCH
RES 1 10000
END  
END  
DISANG=restraints01.rest
LISTIN=POUT   
LISTOUT=POUT"""

    sFile = "equil1_"+self.sMName 
    FEQ = open(sFile+".in","w")
    FEQ.write(sInput)
    FEQ.close()  

    sInput="""NPT MD with position restraints 
 &cntrl
  nmropt = 1, 

  ntx    = 7,       irest  = 1,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 5000,    ntwx   = 5000,    ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 2,       ntb    = 2,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 10,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 0,       
  maxcyc = 5000,   
  ncyc   = 200,   
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 25000,
  nscm   = 0,
  t      = 50.0,    dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,   
  ig     = 71277,
  ntt    = 1,
  tautp  = 0.5,
  
  barostat = 1,
  ntp    = 1,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 2,       tol    = 0.000001,

  watnam='WAT', ! Water residues are named WAT
  owtnm='O',   ! Water oxygens are named O

 &end
 &wt
    type='END'
 &end  
Membrane posres  
25.0
FIND 
"""+sCon+"""
SEARCH
RES 1 10000
END  
END  
DISANG=restraints01.rest
LISTIN=POUT
LISTOUT=POUT"""
     
    for i in range(3):
      sFile = "equil"+str(i+2)+"_"+self.sMName
      FEQ = open(sFile+".in","w")
      FEQ.write(sInput)
      FEQ.close()
    
    iN = 24.0
    for i in range(12):
      sFile = "equil"+str(i+5)+"_"+self.sMName
      FEQ = open(sFile+".in","w")
      iN -= 2.0
      if(i <= 2):
        sRest = "restraints02.rest"
      if(i > 2):
        sRest = "restraints03.rest"
      if(i > 6):
        sRest = "restraints04.rest"
      if(i > 9):
        sRest = "restraints05.rest"
      
      sInput = """NVT MD with position restraints
 &cntrl
  nmropt = 1,

  ntx    = 7,       irest  = 1,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 5000,    ntwx   = 5000,    ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 2,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 10,

  ipol   = 0,

  ibelly = 0,       ntr    = 1,

  imin   = 0,
  maxcyc = 5000,
  ncyc   = 100,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 50000,
  nscm   = 0,
  t      = 100.0,   dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,
  ig     = 71277,
  ntt    = 1, 
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 2,       tol    = 0.000001,

  watnam='WAT', ! Water residues are named WAT
  owtnm='O',   ! Water oxygens are named O

  jfastw = 0,

  ivcap  = 0,       fcap   = 1.5,
 
 &end
 &wt
    type='END'
 &end
Membrane posres
"""+str(iN)+"""
FIND
"""+sCon+"""
SEARCH
RES 1 10000
END  
END      
DISANG="""+sRest+"""
LISTIN=POUT
LISTOUT=POUT"""
      FEQ.write(sInput)
      FEQ.close()
                            
    sInput = """NVT MD without position restraints
 &cntrl
  nmropt = 0,

  ntx    = 7,       irest  = 1,       ntrx   = 1,      ntxo   = 1,
  ntpr   = 5000,    ntwx   = 5000,    ntwv   = 0,      ntwe   = 0,
  ioutfm = 0,       ntwprt = 0,

  ntf    = 2,       ntb    = 1,       dielc  = 1,      igb    = 0,
  cut    = 8.0,     nsnb   = 10,

  ipol   = 0,

  ibelly = 0,       ntr    = 0,

  imin   = 0,
  maxcyc = 5000,
  ncyc   = 200,
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 50000,
  nscm   = 0,
  t      = 200.0,   dt     = 0.002,

  temp0  = 300.0,   tempi  = 300.0,
  ig     = 71277,
  ntt    = 1,
  tautp  = 0.5,
  
  ntp    = 0,       pres0  = 1.0,     comp   = 44.6,
  taup   = 5.0,

  ntc    = 2,       tol    = 0.000001,

  watnam='WAT', ! Water residues are named WAT
  owtnm='O',   ! Water oxygens are named O

 &end
END"""

    sFile = "equil17_"+self.sMName
    FEQ = open(sFile+".in","w")
    FEQ.write(sInput)
    FEQ.close()
    
    FEptraj = open(self.sMName+"_ptraj.in","w")
    FElog = open(self.sMName+"_equil.log","a")
    
    print("  Equilibration in progress... 10%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil0_"+self.sMName+".in -o equil0_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c min3_"+self.sMName+".restrt -r equil0_"+self.sMName+".restrt -ref "+self.sMName+"_solv_ions.inpcrd")[1]+"\n")
    commands.getstatusoutput(AMBERHOME+"/bin/ambpdb -ctr -p "+self.sMName+"_solv_ions.prmtop < equil0_"+self.sMName+".restrt > equil0_"+self.sMName+".pdb")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 15%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil1_"+self.sMName+".in -o equil1_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil0_"+self.sMName+".restrt -r equil1_"+self.sMName+".restrt -ref equil0_"+self.sMName+".restrt -x equil1_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil1_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
        
    print("  Equilibration in progress... 20%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil2_"+self.sMName+".in -o equil2_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil1_"+self.sMName+".restrt -r equil2_"+self.sMName+".restrt -ref equil1_"+self.sMName+".restrt -x equil2_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil2_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
              
    print("  Equilibration in progress... 25%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil3_"+self.sMName+".in -o equil3_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil2_"+self.sMName+".restrt -r equil3_"+self.sMName+".restrt -ref equil2_"+self.sMName+".restrt -x equil3_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil3_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 30%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil4_"+self.sMName+".in -o equil4_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil3_"+self.sMName+".restrt -r equil4_"+self.sMName+".restrt -ref equil3_"+self.sMName+".restrt -x equil4_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil4_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 35%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil5_"+self.sMName+".in -o equil5_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil4_"+self.sMName+".restrt -r equil5_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil5_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil5_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 40%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil6_"+self.sMName+".in -o equil6_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil5_"+self.sMName+".restrt -r equil6_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil6_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil6_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 45%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil7_"+self.sMName+".in -o equil7_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil6_"+self.sMName+".restrt -r equil7_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil7_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil7_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 50%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil8_"+self.sMName+".in -o equil8_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil7_"+self.sMName+".restrt -r equil8_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil8_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil8_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 55%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil9_"+self.sMName+".in -o equil9_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil8_"+self.sMName+".restrt -r equil9_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil9_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil9_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 60%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil10_"+self.sMName+".in -o equil10_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil9_"+self.sMName+".restrt -r equil10_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil10_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil10_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 65%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil11_"+self.sMName+".in -o equil11_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil10_"+self.sMName+".restrt -r equil11_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil11_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil11_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 70%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil12_"+self.sMName+".in -o equil12_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil11_"+self.sMName+".restrt -r equil12_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil12_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil12_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 75%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil13_"+self.sMName+".in -o equil13_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil12_"+self.sMName+".restrt -r equil13_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil13_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil13_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 80%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil14_"+self.sMName+".in -o equil14_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil13_"+self.sMName+".restrt -r equil14_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil14_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil14_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 85%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil15_"+self.sMName+".in -o equil15_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil14_"+self.sMName+".restrt -r equil15_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil15_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil15_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 90%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil16_"+self.sMName+".in -o equil16_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil15_"+self.sMName+".restrt -r equil16_"+self.sMName+".restrt -ref equil4_"+self.sMName+".restrt -x equil16_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil16_"+self.sMName+".mdcrd 1 100 1\n")
    if(self.check_quality() == False):
      return False
    
    print("  Equilibration in progress... 95%")
    os.environ['CUDA_VISIBLE_DEVICES'] = self.sGPU_ID
    FElog.write(commands.getstatusoutput(AMBERHOME+AMBEREXEC+" -O -i equil17_"+self.sMName+".in -o equil17_"+self.sMName+".out -p "+self.sMName+"_solv_ions.prmtop -c equil16_"+self.sMName+".restrt -r equil17_"+self.sMName+".restrt -x equil17_"+self.sMName+".mdcrd")[1]+"\n")
    FEptraj.write("trajin equil17_"+self.sMName+".mdcrd 1 1000 1\n")
    if(self.check_quality() == False):
      return False
      
    print("  Equilibration in progress... 100%")
    print("  Preparing results...")
    
    FEptraj.write("trajout equil_traj_"+self.sMName+".pdb pdb append\n")
    FEptraj.write("center :"+aOff[0]+" mass origin\n")
    FEptraj.write("image origin center\n")
    FEptraj.write("go\n")
    
    FElog.write(commands.getstatusoutput(AMBERHOME+"/bin/cpptraj "+self.sMName+"_solv_ions.prmtop "+self.sMName+"_ptraj.in")[1]+"\n")
    
    FFptraj = open(self.sMName+"_final_ptraj.in","w")
    FFptraj.write("trajin equil17_"+self.sMName+".mdcrd 1 1000 1000\n")
    FFptraj.write("trajout "+self.sMName+"_final.pdb pdb\n")
    FFptraj.write("center :"+aOff[0]+" mass origin\n")
    FFptraj.write("image origin center\n")
    FFptraj.write("go\n")       
    FFptraj.close()
    
    FElog.write(commands.getstatusoutput(AMBERHOME+"/bin/cpptraj "+self.sMName+"_solv_ions.prmtop "+self.sMName+"_final_ptraj.in")[1]+"\n")
    
    FElog.close()
    
    FEptraj.close()
    
    commands.getstatusoutput("cat "+self.sMName+"_final.pdb | /usr/bin/grep -v Na | /usr/bin/grep -v Cl | /usr/bin/grep -v K| /usr/bin/grep -v WAT | /usr/bin/grep -v ' H ' > "+self.sMName+".pdb")    
    
    return self
    
  def restraints(self):

    aForce = ["100","50","25","5","1"]
#    aForce = ["250","100","50","10","2"]

    aAMT = ["O5 C1 C2 C3","C1 C2 C3 C4","C2 C3 C4 C5","C3 C4 C5 O5","C4 C5 O5 C1","C5 O5 C1 C2"]
    aBMT = ["O5 C1 C2 C3","C1 C2 C3 C4","C2 C3 C4 C5","C3 C4 C5 O5","C4 C5 O5 C1","C5 O5 C1 C2"]
    aAMG = ["C1 C2 C3 C4","C2 C3 C4 C5","C3 C4 C5 O5","C4 C5 O5 C1","C5 O5 C1 C2"]
    aBMG = ["C1 C2 C3 C4","C2 C3 C4 C5","C3 C4 C5 O5","C4 C5 O5 C1","C5 O5 C1 C2"]
    aSDS = [""]

    for s in range(len(aForce)):
      iRest = aForce[s]

      str1 = """
    r1=47.5, r2=57.5,
    r3=62.5, r4=72.5,
    rk2="""+iRest+""",
    rk3="""+iRest+""",
/
"""

      str2 = """
    r1=-72.5, r2=-62.5,
    r3=-57.5, r4=-47.5,
    rk2="""+iRest+""",
    rk3="""+iRest+""",
/
"""   
      aOff = dict_lipid[self.sType].split()
      aDihe = []
                                          
      aNr = []
      aID = []

      fR = open("restraints0"+str(s+1)+".rest","w")
      ff = open("min3_"+self.sMName+".pdb","r")

      for line in ff.readlines():
       
       for v in range(len(aOff)-1):
        
        if(aOff[v] == "AMT"):
          aDihe = aAMT
        if(aOff[v] == "BMT"):
          aDihe = aBMT
        if(aOff[v] == "AMG"):
          aDihe = aAMG
        if(aOff[v] == "BMG"):
          aDihe = aBMG
         
        if(line[:4] == "ATOM"):
          ww = string.split(line)
          
          if(ww[3] == aOff[v]): 
            resid = ww[4]
            aID.append(ww[4])
    
            if(ww[4] != aID[len(aID)-2]):
              for i in range(len(aDihe)):
                str3 = ""
                list = aDihe[i].split()
                for j in range(len(list)):
                  for k in range(len(aNr)):
                    if(list[j] == aNr[k].split()[0]):
                      str3 +=aNr[k].split()[1]+","
 
                if(i%2 == 0):
                  fR.write("&rst\n    iat="+str3+str1)
                else:
                  fR.write("&rst\n    iat="+str3+str2)
        
              aNr = []
    
            aNr.append(ww[2]+' '+ww[1])
    
      ff.close()
      fR.close()
        
    return self
