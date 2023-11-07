MicelleMaker

https://pubs.acs.org/doi/10.1021/acsomega.7b00820

written by Dennis M. Krueger, November 2016

Kamerlin Lab, ICM, Uppsala University, Uppsala, Sweden


Lipid library:

HG : Heptyl-a/b-D-glucopyranoside
HM : Heptyl-a/b-D-maltopyranoside
OG : Octyl-a/b-D-glucopyranoside
OM : Octyl-a/b-D-maltopyranoside
NG : Nonyl-a/b-D-glucopyranoside
NM : Nonyl-a/b-D-maltopyranoside
DG : Decyl-a/b-D-glucopyranoside
DM : Decyl-a/b-D-maltopyranoside
UG : Undecyl-a/b-D-glucopyranoside
UM : Undecyl-a/b-D-maltopyranoside
EG : Dodecyl-a/b-D-glucopyranoside
EM : Dodecyl-a/b-D-maltopyranoside
  
SDS : Sodium-dodecyl-sulfate
  
*************************************************************************
  
usage: python micelle_maker.py [-h] 

Syntax description

-h, --help	show this help message and exit

-l	{HG,OG,NG,DG,UG,EG,HM,OM,NM,DM,UM,EM,SDS,UGN,UGV,UGL,UAG,UAA,UAV,UAL,UVG,UVA,UVV,UVL,ULG,ULA,ULV,ULL,SUA,SUV,SUL}	Lipid type

-a	{A,B}	Sugar stereochemistry for glycolipids (A)lpha or (B)eta

-n N	No. of lipids (min. 20 ; max. 200)

-d D	Minimum distance between lipids in A (min. 3 ; max. 4)

-s	{NaCl,KCl,MgCl,CaCl}	Salt type

-c C	Salt concentration in mol/L (min 0.0 ; max 1.0)

-m	Perform minimization

-q	Perform equilibration

License and Disclaimer:

This software is published under a GNU General Public License v3.

As this software was made in part by people employed by Georgia Tech University we must also clarify: “The software is provided “as is.” Neither the Georgia Institute of Technology nor any of its units or its employees, nor the software developers of Micelle Maker or any other person affiliated with the creation, implementation, and upkeep of the software’s code base, knowledge base, and servers (collectively, the “Entities”) shall be held liable for your use of the platform or any data that you enter. The Entities do not warrant or make any representations of any kind or nature with respect to the System, and the Entities do not assume or have any responsibility or liability for any claims, damages, or losses resulting from your use of the platform. None of the Entities shall have any liability to you for use charges related to any device that you use to access the platform or use and receive the platform, including, without limitation, charges for Internet data packages and Personal Computers. THE ENTITIES DISCLAIM ALL WARRANTIES WITH REGARD TO THE SERVICE,INCLUDING WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE TO THE FULLEST EXTENT ALLOWED BY LAW.”
