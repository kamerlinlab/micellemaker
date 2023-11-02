#!/usr/bin/python

import os, sys, math, numpy
from math import sqrt

class superposition():

  def __init__(self,sName,aSele1,aSele2):
     self.name = sName
     self.sele1 = aSele1
     self.sele2 = aSele2

  def super(self):

    aMol = []
    aContent = []
    aSele1 = []
    aSele2 = []
    lenSele1 = 0.0
    lenSele2 = 0.0

    lip_rand = open(self.name+"_random.pdb","r")
    for line in lip_rand.readlines():
      fX = float(line[31:38].strip())
      fY = float(line[39:46].strip())
      fZ = float(line[47:54].strip())
      aMol.append([fX,fY,fZ])
      aContent.append(line)
      if(line[13:16].strip() == self.sele1[0] and line[17:20] == "CLI"):
        aSele1.append([fX,fY,fZ])
        lenSele1 = sqrt(fX**2+fY**2+fZ**2)
      if(line[13:16].strip() == self.sele1[1] and line[17:20] == "CLI"):
        aSele1.append([fX,fY,fZ])
    lip_rand.close()

    tmp_vec = open("tmp_vec.pdb","r")
    for line in tmp_vec.readlines():
      fX = float(line[31:38].strip()) 
      fY = float(line[39:46].strip())     
      fZ = float(line[47:54].strip())
      if(line[13:16].strip() == self.sele2[0] and line[17:20] == "CLI"):
        aSele2.append([fX,fY,fZ])
        lenSele2 = sqrt(fX**2+fY**2+fZ**2)
      if(line[13:16].strip() == self.sele2[1] and line[17:20] == "CLI"): 
        aSele2.append([fX,fY,fZ])
    tmp_vec.close()
    
    # make numpy arrays
    aMol = numpy.array(aMol)
    aSele1 = numpy.array(aSele1)
    aSele2 = numpy.array(aSele2)

    # set sele1 and molecule to origin, sele2 is already in origin
    aMol = numpy.subtract(aMol,aSele1[1])
    aSele1 = numpy.subtract(aSele1,aSele1[1])
    
    # define vectors
    v1 = aSele1[0]
    v2 = aSele2[0]

    # calculate angle between vectors
    cosang = numpy.dot(v1, v2) 
    sinang = numpy.linalg.norm(numpy.cross(v1, v2))
    angle = numpy.arctan2(sinang, cosang)

    def rotation_matrix(axis, theta):
      axis = numpy.asarray(axis)
      axis = axis/math.sqrt(numpy.dot(axis, axis))
      a = math.cos(theta/2.0)
      b, c, d = -axis*math.sin(theta/2.0)
      aa, bb, cc, dd = a*a, b*b, c*c, d*d
      bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
      return numpy.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                          [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                          [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

    # rotation axis, orthonal to v1,v2 plane
    axis = numpy.cross(v1,v2)

    # calculate rotation matrix to v1
    rMat = rotation_matrix(axis,angle)

    #apply rotation to molecule
    aMol = numpy.transpose(aMol)
    aRot = numpy.dot(rMat, aMol)

    #apply translation to molecule
    aRotX = numpy.add(aRot[0],v2[0])
    aRotY = numpy.add(aRot[1],v2[1])
    aRotZ = numpy.add(aRot[2],v2[2])

    fTrans = open(self.name+"_tmp.pdb","w") 
    for i in range(len(aRotX)):
      fX = aRotX[i]
      fY = aRotY[i]
      fZ = aRotZ[i]
      fXs = "%8.3f" % fX
      fYs = "%8.3f" % fY
      fZs = "%8.3f" % fZ
      fTrans.write(aContent[i][:30]+fXs+fYs+fZs+aContent[i][54:])
    fTrans.close()

