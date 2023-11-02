#!/usr/bin/python
# -*- coding: utf-8 -*- 

import os, sys, math, random, argparse, time, shutil
from math import sin, asin, cos, acos, pi, sqrt, ceil, degrees, radians
from amber import *
from super import *

commands.getstatusoutput("ulimit")

start_time = time.time()

__author__ = 'Dennis M. Krueger'
__copyright__ = 'Copyright 2016, Dennis M. Krueger'
__version__ = '1.0'

class MyParser(argparse.ArgumentParser):
  def error(self, message): 
    sys.stderr.write('error: %s\n' % message)
    self.print_help()
    sys.exit(2)

def float_number(x):
  x = float(x)
  if x > 4.0:
    raise argparse.ArgumentTypeError("Maximum distance between lipids is 4 A")
  if x < 3.0:
      raise argparse.ArgumentTypeError("Minimum distance between lipids is 3 A")
  return x

def float_conc(x):
  x = float(x)
  if x > 1.0:    
    raise argparse.ArgumentTypeError("Maximum salt concentration is 1 mol/L")                        
  if x < 0.0:
    raise argparse.ArgumentTypeError("Minimum salt concentration is 0 mol/L")
  return x

def int_lipno(x):
  x = int(x)              
  if x > 200:               
    raise argparse.ArgumentTypeError("Maximum no. of lipids is 200")                             
  if x < 20:               
      raise argparse.ArgumentTypeError("Minimum no. of lipids is 20") 
  return x

parser=MyParser(description='Syntax description')
parser.add_argument('-l', help='Lipid type', required=True, choices=['HG','OG','NG','DG','UG','EG','HM','OM','NM','DM','UM','EM','SDS','UGN','UGV','UGL','UAG','UAA','UAV','UAL','UVG','UVA','UVV','UVL','ULG','ULA','ULV','ULL','SUA','SUV','SUL'])
parser.add_argument('-a', help='Sugar stereochemistry for glycolipids (A)lpha or (B)eta', choices=['A','B'])
parser.add_argument('-n', type=int_lipno, help='No. of lipids (min. 20 ; max. 200)', required=True)
parser.add_argument('-d', type=float_number, help='Minimum distance between lipids in A (min. 3 ; max. 4)', required=True)
parser.add_argument('-s', help='Salt type', required=True, choices=['NaCl','KCl','MgCl','CaCl'])
parser.add_argument('-c', type=float_conc, help='Salt concentration in mol/L (min 0.0 ; max 1.0)', required=True)
parser.add_argument('-m', help='Perform minimization', action='store_true', required=False)
parser.add_argument('-q', help='Perform equilibration', action='store_true', required=False)
parser.add_argument('-f', help=argparse.SUPPRESS, required=False)
parser.add_argument('-g', help=argparse.SUPPRESS, required=False) 

print("""
  __  __ _          _ _        __  __       _             
 |  \/  (_)        | | |      |  \/  |     | |            
 | \  / |_  ___ ___| | | ___  | \  / | __ _| | _____ _ __ 
 | |\/| | |/ __/ _ \ | |/ _ \ | |\/| |/ _` | |/ / _ \ '__|
 | |  | | | (_|  __/ | |  __/ | |  | | (_| |   <  __/ |   
 |_|  |_|_|\___\___|_|_|\___| |_|  |_|\__,_|_|\_\___|_|   
 \n ____________________________________________\n
  written by Dennis M. Krueger, November 2016\n
  ICM, Uppsala University, Uppsala, Sweden\n

  Lipid library:

  HG : Heptyl-a/b-D-glucopyranoside     HM : Heptyl-a/b-D-maltopyranoside
  OG : Octyl-a/b-D-glucopyranoside      OM : Octyl-a/b-D-maltopyranoside
  NG : Nonyl-a/b-D-glucopyranoside      NM : Nonyl-a/b-D-maltopyranoside
  DG : Decyl-a/b-D-glucopyranoside      DM : Decyl-a/b-D-maltopyranoside
  UG : Undecyl-a/b-D-glucopyranoside    UM : Undecyl-a/b-D-maltopyranoside
  EG : Dodecyl-a/b-D-glucopyranoside    EM : Dodecyl-a/b-D-maltopyranoside
  
  SDS : Sodium-dodecyl-sulfate
  
  *************************************************************************
  """)

try:
  args = parser.parse_args()
except:
  print("\n")
  sys.exit()

aAAS = ['UGN','UGV','UGL','UAG','UAA','UAV','UAL','UVG','UVA','UVV','UVL','ULG','ULA','ULV','ULL','SUA','SUV','SUL']

if(args.l == "SDS" or args.l in aAAS):
  sType = args.l
  if(args.a != ""):
    print("  No stereochemistry required, argument -a will be ignored")
else:
  if(args.a != "A" and args.a != "B"):
    print("\n")
    parser.error("Stereochemistry required")
  else:
    sType = args.l+args.a #"EMA"

if(str(args.g) == "None"):
  sGPU_ID = "0"
else:
  sGPU_ID = args.g

sSaltNo = ""
iNo = int(args.n)   #100
fDist = args.d #2.0
if(args.s[:-2] == "Mg" or args.s[:-2] == "Ca"):
  sIonP = args.s[:-2]+"2+"
  sSaltNo = "₂"
else:
  sIonP = args.s[:-2]+"+"
sIonN = args.s[-2:]+"-"
sConc = str(args.c)

sPath = os.getcwd()+"/molecules/"+sType+"/"
fArea = (fDist**2)*iNo
fRad = sqrt(fArea/(4*pi)) # sphere surface area = 4*pi*r^2
fRad = ceil(fRad)
#fAngle = asin(fDist/fRad)
#iRot = int(ceil(360.0/degrees(fAngle)))

sMName = "micelle_"+sType+"_"+str(iNo)+"_"+str(fDist)+"_"+args.s+"_"+sConc
sDate = time.strftime('%x').split('/')
if(str(args.f) == "None"):
  sFolder = sDate[1]+sDate[0]+sDate[2]+"_"+time.strftime('%X').replace(':','')+"_"+sMName
else:
  sFolder = args.f

print("  Lipid: "+sType)
print("  No.lipids: "+str(iNo))
print("  Lipid distance: "+str(fDist)+" A")
#print("  Core surface area: "+str(fArea)+" A^3")
#print("  Core radius: "+str(fRad)+" A")
print("  Salt: "+str(args.s))+sSaltNo
print("  Salt concentration: "+sConc+" mol/L")
print("\n")
#print "Rotation angle: "+str(round(degrees(fAngle),2))
#print "No.rotations: "+str(iRot)

aLipid = []
aLipidC = []

def fReadLipid(sLipPath,sLipType):
  
  sContent = os.listdir(sLipPath)
  iN = 0
  for i in range(len(sContent)):
    if(sContent[i][-4:] == ".pdb"):
      iN+=1
  iNoC = 0 
  file = open(sLipPath+sLipType+".pdb",'r')
  for line in file.readlines():
    fX = line[31:38].strip()
    fY = line[39:46].strip()
    fZ = line[47:54].strip()    
    aLipidC.append(fX+" "+fY+" "+fZ)
    aLipid.append(line[:-1])
    if(line[17:20] == "CLI"):
      iNoC+=1
  file.close()
  return iNoC
    
def GetSphereCoords(iNoPoints):
    #each point will be of form 'x, y, z'; in cartesian coordinates, the distance from the origion [0., 0., 0.] for each point will be 1.0 
    #converted from:  http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere ) 
    dGold = pi*(3.0-sqrt(5.0))  # ~2.39996323 
    dZ   =  2.0/iNoPoints
    dL = random.random() * iNoPoints
    dZC    =  1.0 - dZ/2.0
    ptsOnSphere =[]
    for k in range(iNoPoints): 
      dR  = sqrt(1.0-dZC*dZC)
      fX = cos(dL)*dR
      fY = sin(dL)*dR
      fZ = dZC
      fS = sqrt(fX**2+fY**2+fZ**2)    
      fX = fX/fS*fRad
      fY = fY/fS*fRad
      fZ = fZ/fS*fRad
      fXs = "%8.3f" % fX
      fYs = "%8.3f" % fY
      fZs = "%8.3f" % fZ
      ptsOnSphere.append(fXs+" "+fYs+" "+fZs)
      dZC    = dZC - dZ
      dL = dGold + dL
    
    fA = random.uniform(-180.0,180.0)
    fB = random.uniform(-180.0,180.0)
    fG = random.uniform(-180.0,180.0)
    
    ptsOnSphereR =[]
    for i in range(len(ptsOnSphere)):
      fX = float(ptsOnSphere[i].split()[0])
      fY = float(ptsOnSphere[i].split()[1])
      fZ = float(ptsOnSphere[i].split()[2])

      fXsA = fX
      fYsA = fY*cos(fA)-fZ*sin(fA)
      fZsA = fY*sin(fA)+fZ*cos(fA)
    
      fXsB = fXsA*cos(fB)+fZsA*sin(fB)
      fYsB = fYsA
      fZsB = -fXsA*sin(fB)+fZsA*cos(fB) 
    
      fXsG = fXsB*cos(fG)-fYsB*sin(fG)
      fYsG = fXsB*sin(fG)+fYsB*cos(fG)
      fZsG = fZsB

      fXs = "%8.3f" % fXsG
      fYs = "%8.3f" % fYsG
      fZs = "%8.3f" % fZsG
      ptsOnSphereR.append(fXs+" "+fYs+" "+fZs)      
    
    return ptsOnSphereR
     
def RandomRotationLipid(aLip):

  fA = random.uniform(-180.0,180.0)
  fB = random.uniform(-180.0,180.0)
  fG = random.uniform(-180.0,180.0)
  
  FInLip = open(aLip,"r")
  FOutLip = open(sType+"_random.pdb",'w')
  for line in FInLip.readlines():
    fX = float(line[31:38].strip())
    fY = float(line[39:46].strip())
    fZ = float(line[47:54].strip())

    fXsA = fX
    fYsA = fY*cos(fA)-fZ*sin(fA)
    fZsA = fY*sin(fA)+fZ*cos(fA)
    
    fXsB = fXsA*cos(fB)+fZsA*sin(fB)
    fYsB = fYsA
    fZsB = -fXsA*sin(fB)+fZsA*cos(fB) 
    
    fXsG = fXsB*cos(fG)-fYsB*sin(fG)
    fYsG = fXsB*sin(fG)+fYsB*cos(fG)
    fZsG = fZsB

    fXs = "%8.3f" % fXsG
    fYs = "%8.3f" % fYsG
    fZs = "%8.3f" % fZsG
        
    if(len(fXs) == 5):
      fXs = "   "+fXs
    if(len(fXs) == 6):
      fXs = "  "+fXs
    if(len(fXs) == 7):
      fXs = " "+fXs
    if(len(fYs) == 5):
      fYs = "   "+fYs
    if(len(fYs) == 6):
      fYs = "  "+fYs
    if(len(fYs) == 7):
      fYs = " "+fYs
    if(len(fZs) == 5):
      fZs = "   "+fZs
    if(len(fZs) == 6):
      fZs = "  "+fZs
    if(len(fZs) == 7):
      fZs = " "+fZs
    
    FOutLip.write(line[:30]+fXs+fYs+fZs+line[54:])
  FOutLip.close()
  FInLip.close()

def get_diameter(sName):
  aPDB = []
  FPDB = open(sName+".pdb",'r')
  for line in FPDB.readlines():
    if(line[:4] == "ATOM"):
      fX = line[31:38].strip()
      fY = line[39:46].strip()
      fZ = line[47:54].strip() 
      aPDB.append(fX+" "+fY+" "+fZ)
  FPDB.close()
  fDiam = 0.0
  for i in range(len(aPDB)):
    for j in range(len(aPDB)):
      fDX = float(aPDB[i].split()[0]) - float(aPDB[j].split()[0])
      fDY = float(aPDB[i].split()[1]) - float(aPDB[j].split()[1])
      fDZ = float(aPDB[i].split()[2]) - float(aPDB[j].split()[2])
      fDiff = sqrt(fDX**2+fDY**2+fDZ**2)
      if(fDiff > fDiam):
        fDiam = fDiff    
  return fDiam

def check_clashes(lipid,neighbours,fRmsd):
  
  bFail = False
  aL = []
  aN = []
  
  F1 = open(lipid,"r")
  for line in F1.readlines():
    if(line[:4] == "ATOM"):# and line[17:20] != "CLI"):
      fX = line[31:38].strip()
      fY = line[39:46].strip()
      fZ = line[47:54].strip()
      aL.append(fX+" "+fY+" "+fZ)
  F1.close
  
  F2 = open(neighbours,"r")
  for line in F2.readlines():
    if(line[:4] == "ATOM"):# and line[17:20] != "CLI"):
      fX = line[31:38].strip()       
      fY = line[39:46].strip()       
      fZ = line[47:54].strip()       
      aN.append(fX+" "+fY+" "+fZ)  
  F2.close()

  for i in range(len(aL)):
    for j in range(len(aN)):
      fXL = float(aL[i].split()[0])
      fYL = float(aL[i].split()[1])
      fZL = float(aL[i].split()[2])
  
      fXN = float(aN[j].split()[0])
      fYN = float(aN[j].split()[1])
      fZN = float(aN[j].split()[2])
      
      rms = sqrt( (fXL-fXN)**2 + (fYL-fYN)**2 + (fZL-fZN)**2 )
      
      if(rms < fRmsd): #larger value means more sensitive to clashes!
        bFail = True
        break
  
  return bFail

def build_micelle():
  
  print("  Setting up micelle...\n") 
  
  iCcount = fReadLipid(sPath,sType)
  
  aSphereCoords = GetSphereCoords(iNo)

  sCD1 = ""
  sCD10 = ""
  iL1 = 0
  iL2 = 0
  sCcount_l = "2" 
  sCcount_h = str(iCcount)
    
  for i in range(len(aLipid)):
    
    if(aLipid[i][13:16].strip() == "C"+sCcount_h and aLipid[i][17:20].strip() == "CLI"):
      sCD1 = aLipid[i][31:38].strip()+" "+aLipid[i][39:46].strip()+" "+aLipid[i][47:54].strip()  
      iL1 = i
    if(aLipid[i][13:16].strip() == "C"+sCcount_l and aLipid[i][17:20].strip() == "CLI"):
      sCD10 = aLipid[i][31:38].strip()+" "+aLipid[i][39:46].strip()+" "+aLipid[i][47:54].strip()
      iL2 = i

  fLen = sqrt( (float(sCD10.split()[0])-float(sCD1.split()[0]) )**2 + ( float(sCD10.split()[1])-float(sCD1.split()[1]) )**2 + (float(sCD10.split()[2])-float(sCD1.split()[2]) )**2 )
  fMove = fLen+fRad

  FMicelle = open(sMName+".pdb","w")
  FMicelle.close()
  FSPhere = open("sphere.pdb","w")
  
  iK = 0
  
  for j in range(len(aSphereCoords)):
    
    fX1 = float(aSphereCoords[j].split()[0])
    fY1 = float(aSphereCoords[j].split()[1])
    fZ1 = float(aSphereCoords[j].split()[2])
    
    sfX1 = "%8.3f" % (fX1)
    sfY1 = "%8.3f" % (fY1)
    sfZ1 = "%8.3f" % (fZ1)
            
    fVB = sqrt(fX1**2+fY1**2+fZ1**2)
    fX1n = fX1/fVB*fMove
    fY1n = fY1/fVB*fMove
    fZ1n = fZ1/fVB*fMove  
    
    sfX2 = "%8.3f" % (fX1n)
    sfY2 = "%8.3f" % (fY1n)                       
    sfZ2 = "%8.3f" % (fZ1n)
    
    if(len(sfX1) == 5):
      sfX1 = "   "+sfX1
    if(len(sfX1) == 6):
      sfX1 = "  "+sfX1
    if(len(sfX1) == 7):
      sfX1 = " "+sfX1
    if(len(sfY1) == 5):
      sfY1 = "   "+sfY1
    if(len(sfY1) == 6):
      sfY1 = "  "+sfY1
    if(len(sfY1) == 7):
      sfY1 = " "+sfY1
    if(len(sfZ1) == 5):
      sfZ1 = "   "+sfZ1
    if(len(sfZ1) == 6):
      sfZ1 = "  "+sfZ1
    if(len(sfZ1) == 7):
      sfZ1 = " "+sfZ1
    
    if(len(sfX2) == 5):
      sfX2 = "   "+sfX2
    if(len(sfX2) == 6):
      sfX2 = "  "+sfX2
    if(len(sfX2) == 7):
      sfX2 = " "+sfX2
    if(len(sfY2) == 5):
      sfY2 = "   "+sfY2
    if(len(sfY2) == 6):
      sfY2 = "  "+sfY2
    if(len(sfY2) == 7):
      sfY2 = " "+sfY2
    if(len(sfZ2) == 5):
      sfZ2 = "   "+sfZ2
    if(len(sfZ2) == 6):
      sfZ2 = "  "+sfZ2
    if(len(sfZ2) == 7):
      sfZ2 = " "+sfZ2
    
    FTmp = open("tmp_vec.pdb","w")
    if(int(sCcount_h) > 9):
      FTmp.write(aLipid[iL1][:13]+"C"+sCcount_h+aLipid[iL1][16:30]+sfX1+sfY1+sfZ1+aLipid[iL1][54:]+"\n")
    else:
      FTmp.write(aLipid[iL1][:13]+"C"+sCcount_h+" "+aLipid[iL1][16:30]+sfX1+sfY1+sfZ1+aLipid[iL1][54:]+"\n")      
    FTmp.write(aLipid[iL2][:13]+"C"+sCcount_l+" "+aLipid[iL2][16:30]+sfX2+sfY2+sfZ2+aLipid[iL2][54:]+"\n")
    FTmp.close()
    
    FSPhere.write(aLipid[iL1][:13]+"C"+sCcount_h+aLipid[iL1][17:30]+sfX1+sfY1+sfZ1+aLipid[iL1][54:]+"\n")
    FSPhere.write(aLipid[iL2][:13]+"C"+sCcount_l+" "+aLipid[iL2][17:30]+sfX2+sfY2+sfZ2+aLipid[iL2][54:]+"\n")
    
    RandomRotationLipid(sPath+sType+".pdb")
    iTrans1 = random.randrange(int(sCcount_h)-1,int(sCcount_h)+1)
    iTrans2 = iTrans1-int(sCcount_h)+int(sCcount_l)

#    iTrans1 = sCcount_h 
#    iTrans2 = sCcount_l
    
    sele1 = ["C"+str(iTrans2),"C"+str(iTrans1)]
    sele2 = ["C"+sCcount_l,"C"+sCcount_h]
    
    superposition(sType,sele1,sele2).super()
    
    iFail = 0
    fRmsd = 2.5
    while True:
      bClash = check_clashes(sType+"_tmp.pdb",sMName+".pdb",fRmsd)
      
      if(iFail > 20):
        fRmsd = fRmsd-0.1
        iFail = 0
#        print("  Found clash for lipid "+str(j+1)+", change threshold to "+str(fRmsd) ) 
#        print("  Please check clashes for lipid "+str(j+1)+".")
#        break
        
      if(bClash == True):
        iFail +=1
#        print("  Found clash for lipid "+str(j+1)+", replacing..")
        RandomRotationLipid(sPath+sType+".pdb")
    
#        iTrans1 = sCcount_h 
#        iTrans2 = sCcount_l
        
        iTrans1 = random.randrange(int(sCcount_h)-1,int(sCcount_h)+1)
        iTrans2 = iTrans1-int(sCcount_h)+int(sCcount_l)
        
        sele1 = ["C"+str(iTrans2),"C"+str(iTrans1)]
        sele2 = ["C"+sCcount_l,"C"+sCcount_h]
        
        superposition(sType,sele1,sele2).super()    
      
      else:
        break
    
    print "  Placed lipid "+str(j+1)
      
    iZ = ""
    iK += 1
    if(iK >= 10):
      iZ = " "+str(iK)
    if(iK < 10):
      iZ = "  "+str(iK)
    if(iK >= 100):
      iZ = str(iK)
    
    FMicelle = open(sMName+".pdb","a")  
    FDDMtmp = open(sType+"_tmp.pdb","r")
    for line in FDDMtmp.readlines():
      if(line[:4] == "ATOM"):
        FMicelle.write(line[:23]+iZ+line[26:])
    FMicelle.write("\nTER\n")
    FDDMtmp.close()
    FMicelle.close()
    shutil.copyfile(sMName+".pdb", sMName+"_init.pdb")
    
  FSPhere.close()
    
  os.remove("tmp_vec.pdb")
  os.remove(sType+"_tmp.pdb")  
  os.remove(sType+"_random.pdb")
  os.remove("sphere.pdb")
  
  print("  Micelle generated!\n")
  
  return sMName,fMove

def check_quality():
  try:
    os.rename('leap.log','leap.logfile')
  except:
    ""
  print("  Check quality of the structure...")
  iA = int(commands.getstatusoutput("/usr/bin/grep FAILURE min2*.out | /usr/bin/wc -l")[1])
  print("  "+str(iA)+" minimization errors.")
  iB = int(commands.getstatusoutput("/usr/bin/grep '\*\*\*\*\*\*\*\*\*\*\*\*\*' min2*.out | /usr/bin/wc -l")[1])
  print("  "+str(iB)+" energy errors.")
  iC = int(commands.getstatusoutput("/usr/bin/grep NaN min2*.out | /usr/bin/wc -l")[1])
  print("  "+str(iC)+" coordination errors.")
  try:
    iD = int(commands.getstatusoutput("/usr/bin/grep ERROR *.log | /usr/bin/wc -l")[1])
  except:
    iD = 0
  print("  "+str(iD)+" density errors.")
  iS = iA + iB + iC + iD
  if(iS > 0):
    print("  Input structure has bad quality. Restart.\n")
    return False
  else:
    print("  Structure has good quality. Continue...\n")
    return True

def build_system():
  
  try:
    os.system("/usr/bin/rm *.log")
  except:
    ""
  try:
    os.system("/usr/bin/rm *.out")
  except:
    ""
  
  while(True):      
    sMName,fMove = build_micelle()

    #fDiam = get_diameter(sMName)

    sDiam = str(int(fMove*1.5*2))
    sBoxLen = str(int(fMove*1.5*2+2*11))
    sBoxDim = sBoxLen+" "+sBoxLen+" "+sBoxLen
    #print("  Micelle diameter: ~"+str(sDiam)+" A")
    print("  Box size: "+str(sBoxDim)+" A\n")

    MD = do_amber(sType,sMName,sIonP,sIonN,sConc,sBoxDim,iNo,sGPU_ID)
    MD.genLeapFile()
    
    #print("  Run Leap...")
    print("  Setting up box...")
    MD.runLeap()
    print("  Done!\n")

    if(args.m):
      if(args.q == False):
        print("  Run minimization...")
        MD.minimize(10000)
      print("  Done!\n")

    if(args.q):
      
      print("  Run minimization...")
      MD.minimize(10000)
      print("  Done!\n") 
      
      if(check_quality() == False):
        return False
        break
      
      #MD.restraints() 
      
      print("  Run equilibration...")# part 1 of 2...")
      if(MD.equilibration() == False):
        return False
        break
#      MD.runLeap()
#      MD.minimize(1000)
      print("  Done!\n")
#      print("  Run equilibration part 2 of 2...")
#      if(MD.equilibration() == False):
#        return False
#        break
#      print("  Done!\n")
    
    print("  Compressing output files...")
    commands.getstatusoutput("/usr/bin/gzip -9 *.inpcrd *.restrt *.mdcrd")
    print("  Done!\n")

    return True
    break
  
  return True

def main():
  os.mkdir(sFolder)
  os.chdir(sFolder)
  iExit = 0
  global fDist
  while(True):
    print("  Trial "+str(iExit+1)+"\n")
    if(build_system() == True):
      print("  Trial "+str(iExit+1)+" successful!\n")
      break
    else:
      iExit += 1
    if(iExit > 2):
      if(fDist+1.0 <= 4): 
        fDist += 1.0
      else:
        fDist -= 1.0
        print("  Distance between lipids changed to "+str(fDist)+"!\n")
    if(iExit > 4):
      print("  Micelle generation failed too often. Please check input parameters!")
      break

main()
print("  Runtime: %s minutes \n" % round(((time.time() - start_time)/60),2) )
