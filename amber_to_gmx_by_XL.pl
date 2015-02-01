#!/usr/bin/perl 

#*************************************************************************
#                ==== AMBER >>> to >>> GMX ====
#  Convert AMBER top and crd files to GROMACS format
#  Compatible with GLYCAM force field 
#
#  Copyright (C) 2012-2015 Xin Li <lixin.reco@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#*************************************************************************

   use strict;
   use warnings;

   printf "\n";
   printf "***********************************************\n";
   printf "*       ==== AMBER >>> to >>> GMX ====        *\n";
   printf "*            Copyright (C) Xin Li             *\n";
   printf "*                                             *\n";
   printf "*      This program is free software and      *\n";
   printf "*      comes with ABSOLUTELY NO WARRANTY.     *\n";
   printf "*                                             *\n";
   printf "*               Licensed under                *\n";
   printf "*         GNU General Public License          *\n";
   printf "*                (Version 2)                  *\n";
   printf "***********************************************\n";
   printf "\n";
   printf "    Usage: perl amber_to_gmx_by_XL.pl \\\n";
   printf "           -top [name_of_amber_top_file] \\\n";
   printf "           -crd [name_of_amber_crd_file] \\\n";
   printf "           -name [name_of_your_molecule] \\\n";
   printf "           -glycam [yes_or_no]\n";
   printf "\n";

   if (@ARGV==0) { exit; }

#  options

   my %options = @ARGV;
   die "Error: undefined -top option!\n" unless (defined $options{"-top"});
   my $amb_top = $options{"-top"};

   die "Error: undefined -crd option!\n" unless (defined $options{"-crd"});
   my $amb_crd = $options{"-crd"};

   die "Error: undefined -name option!\n" unless (defined $options{"-name"});
   my $gmx_top = "gmx_".$options{"-name"}.".top";
   my $gmx_gro = "gmx_".$options{"-name"}.".gro";

   my $use_glycam = 0;
   if ( defined $options{"-glycam"} )
   {
      $use_glycam = $options{"-glycam"};
      if ( $use_glycam eq "yes" or $use_glycam eq "no" )
      {
         printf "Use glycam 1-4 scaling factors? $use_glycam\n";
      }
      else
      {
         die "Error: invalid -glycam option (should be \"yes\" or \"no\")!\n";
      }
   }
   else
   {
      printf "Will NOT use glycam 1-4 scaling factors. (Default)\n";
   }
   printf "\n";

#  constant: pi
   my $pi = 3.14159265358979323846;

#*********************************************************
#  Read amber top file
#*********************************************************

   my $string = "";
   open ATOP,"$amb_top" or die 
      "Error: cannot open amber top file $amb_top !\n";
   while ( <ATOP> ) 
   {
#     read flags
      if ( /^\%/ ) 
      {
         if ( /^\%FLAG/ ) 
         {
            chomp;
            s/^\%FLAG/<FLAG>/;
            s/\s+//g;
            $string .= $_."<FLAG> ";
         }
      }
#     read data 
      else
      {
         chomp;
         s/\s+/ /g;
         $string .= $_." ";
      }

   }
   close ATOP;

#  remove the first <FLAG> 
   $string =~ s/^<FLAG>//;
#  split the string 
   my @array = split "<FLAG>",$string;
#  change the array into hash
   my %hash = @array;

#
#  For amber top format, 
#  see http://ambermd.org/formats.html#topology
#

#  POINTERS
   $_ = $hash{"POINTERS"};
   @_ = (split);

   my $NATOM    = $_[0] ;  # total number of atoms 
   my $NTYPES   = $_[1] ;  # total number of distinct atom types (distinct LJ params - X.L.)
   my $NBONH    = $_[2] ;  # number of bonds containing hydrogen
   my $MBONA    = $_[3] ;  # number of bonds not containing hydrogen
   my $NTHETH   = $_[4] ;  # number of angles containing hydrogen
   my $MTHETA   = $_[5] ;  # number of angles not containing hydrogen
   my $NPHIH    = $_[6] ;  # number of dihedrals containing hydrogen
   my $MPHIA    = $_[7] ;  # number of dihedrals not containing hydrogen
   my $NHPARM   = $_[8] ;  # currently not used
   my $NPARM    = $_[9] ;  # used to determine if addles created prmtop
   my $NNB      = $_[10];  # number of excluded atoms
   my $NRES     = $_[11];  # number of residues
   my $NBONA    = $_[12];  # MBONA + number of constraint bonds
   my $NTHETA   = $_[13];  # MTHETA + number of constraint angles
   my $NPHIA    = $_[14];  # MPHIA + number of constraint dihedrals
   my $NUMBND   = $_[15];  # number of unique bond types
   my $NUMANG   = $_[16];  # number of unique angle types
   my $NPTRA    = $_[17];  # number of unique dihedral types
   my $NATYP    = $_[18];  # number of atom types in parameter file, see SOLTY below
   my $NPHB     = $_[19];  # number of distinct 10-12 hydrogen bond pair types
   my $IFPERT   = $_[20];  # set to 1 if perturbation info is to be read in
   my $NBPER    = $_[21];  # number of bonds to be perturbed
   my $NGPER    = $_[22];  # number of angles to be perturbed
   my $NDPER    = $_[23];  # number of dihedrals to be perturbed
   my $MBPER    = $_[24];  # number of bonds with atoms completely in perturbed group
   my $MGPER    = $_[25];  # number of angles with atoms completely in perturbed group
   my $MDPER    = $_[26];  # number of dihedrals with atoms completely in perturbed groups
   my $IFBOX    = $_[27];  # set to 1 if standard periodic box, 2 when truncated octahedral
   my $NMXRS    = $_[28];  # number of atoms in the largest residue
   my $IFCAP    = $_[29];  # set to 1 if the CAP option from edit was specified
   my $NUMEXTRA = $_[30];  # number of extra points found in topology
   my $NCOPY    = $_[31];  # number of PIMD slices / number of beads

#  print information on screen

   printf "%-10s%-10s%-10s%-10s%-10s\n", 
      "Residues", "Atoms", "Bonds", "Angles", "Dihedrals";
   printf "%-10d%-10d%-10d%-10d%-10d\n", 
      $NRES, $NATOM, $NBONH+$NBONA, $NTHETH+$NTHETA, $NPHIH+$NPHIA; 
   printf "\n";

   printf "%-12s%-12s%-12s%-12s\n", 
      "LJ-Types", "BondTypes", "AngleTypes", "DihedralTypes";
   printf "%-12d%-12d%-12d%-12d\n", 
      $NTYPES, $NUMBND, $NUMANG, $NPTRA;
   printf "\n";

#*********************************************************
#  Read ATOM_NAME, ATOM_TYPE_INDEX, 
#       NONBONDED_PARM_INDEX,
#       LENNARD_JONES_ACOEF, 
#       LENNARD_JONES_BCOEF,
#       AMBER_ATOM_TYPE
#*********************************************************

#  ATOM_NAME (Note: add an arbitrary element to $_ 
#             so that the array starts from index 1)
   $_ = "zero  ".$hash{"ATOM_NAME"};
#  save data in corresponding array
   my @IGRAPH = (split);
#  check the size of the array
   die "Error: incorrect number of ATOM_NAME!\n"
      if ( $NATOM != @IGRAPH-1 );

#  ATOM_TYPE_INDEX
   $_ = "zero  ".$hash{"ATOM_TYPE_INDEX"};
   my @IAC = (split);
   die "Error: incorrect number of ATOM_TYPE_INDEX!\n"
      if ( $NATOM != @IAC-1 );

#  NONBONDED_PARM_INDEX
   $_ = "zero  ".$hash{"NONBONDED_PARM_INDEX"};
   my @ICO = (split);
   die "Error: incorrect number of NONBONDED_PARM_INDEX!\n"
      if ( $NTYPES*$NTYPES != @ICO-1 );

#  LENNARD_JONES_ACOEF
   $_ = "zero  ".$hash{"LENNARD_JONES_ACOEF"};
   my @CN1 = (split);
   die "Error: incorrect number of LENNARD_JONES_ACOEF!\n"
      if ( $NTYPES*($NTYPES+1)/2 != @CN1-1 );

#  LENNARD_JONES_BCOEF
   $_ = "zero  ".$hash{"LENNARD_JONES_BCOEF"};
   my @CN2 = (split);
   die "Error: incorrect number of LENNARD_JONES_BCOEF!\n"
      if ( $NTYPES*($NTYPES+1)/2 != @CN2-1 );

#  AMBER_ATOM_TYPE
   $_ = "zero  ".$hash{"AMBER_ATOM_TYPE"};
   my @ISYMBL = (split);
   die "Error: incorrect number of AMBER_ATOM_TYPE!\n"
      if ( $NATOM != @ISYMBL-1 );

#*********************************************************
#  Read CHARGE, MASS, RESIDUE_LABEL, RESIDUE_POINTER
#*********************************************************

#  CHARGE
   $_ = "zero  ".$hash{"CHARGE"};
   my @CHARGE = (split);
   die "Error: incorrect number of CHARGE!\n"
      if ( $NATOM != @CHARGE-1 );

#  MASS
   $_ = "zero  ".$hash{"MASS"};
   my @AMASS = (split);
   die "Error: incorrect number of MASS!\n"
      if ( $NATOM != @AMASS-1 );

#  RESIDUE_LABEL
   $_ = "zero  ".$hash{"RESIDUE_LABEL"};
   my @LBRES = (split);
   die "Error: incorrect number of RESIDUE_LABEL!\n"
      if ( $NRES != @LBRES-1 );

#  RESIDUE_POINTER
   $_ = "zero  ".$hash{"RESIDUE_POINTER"};
   my @IPRES = (split);
   push @IPRES, $NATOM+1; # add NATOM+1 to @IPRES
   die "Error: incorrect number of RESIDUE_POINTER!\n"
      if ( $NRES != @IPRES-2 );

#*********************************************************
#  Read BOND_FORCE_CONSTANT, BOND_EQUIL_VALUE
#       BONDS_INC_HYDROGEN, BONDS_WITHOUT_HYDROGEN
#*********************************************************

   printf "Reading bonds ...\n";

   my $ibond;

#  BOND_FORCE_CONSTANT
   $_ = "zero  ".$hash{"BOND_FORCE_CONSTANT"};
   my @RK = (split);
   die "Error: incorrect number of BOND_FORCE_CONSTANT!\n"
      if ( $NUMBND != @RK-1 );

#  BOND_EQUIL_VALUE
   $_ = "zero  ".$hash{"BOND_EQUIL_VALUE"};
   my @REQ = (split);
   die "Error: incorrect number of BOND_EQUIL_VALUE!\n"
      if ( $NUMBND != @REQ-1 );

#  BONDS_INC_HYDROGEN
   $_ = "zero  ".$hash{"BONDS_INC_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of BONDS_INC_HYDROGEN!\n"
      if ( $NBONH*3 != @_-1 );

   my @IBH; my @JBH; my @ICBH;

   for $ibond (1..$NBONH)
   {
      $IBH[$ibond]  = $_[($ibond-1)*3+1];
      $JBH[$ibond]  = $_[($ibond-1)*3+2];
      $ICBH[$ibond] = $_[($ibond-1)*3+3];
   }

#  BONDS_WITHOUT_HYDROGEN
   $_ = "zero  ".$hash{"BONDS_WITHOUT_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of BONDS_WITHOUT_HYDROGEN!\n"
      if ( $NBONA*3 != @_-1 );

   my @IB; my @JB; my @ICB;

   for $ibond (1..$NBONA)
   {
      $IB[$ibond]  = $_[($ibond-1)*3+1];
      $JB[$ibond]  = $_[($ibond-1)*3+2];
      $ICB[$ibond] = $_[($ibond-1)*3+3];
   }

#*********************************************************
#  Read ANGLE_FORCE_CONSTANT, ANGLE_EQUIL_VALUE,
#       ANGLES_INC_HYDROGEN, ANGLES_WITHOUT_HYDROGEN
#*********************************************************

   printf "Reading angles ...\n";

   my $iang;

#  ANGLE_FORCE_CONSTANT
   $_ = "zero  ".$hash{"ANGLE_FORCE_CONSTANT"};
   my @TK = (split);
   die "Error: incorrect number of ANGLE_FORCE_CONSTANT!\n"
      if ( $NUMANG != @TK-1 );

#  ANGLE_EQUIL_VALUE
   $_ = "zero  ".$hash{"ANGLE_EQUIL_VALUE"};
   my @TEQ = (split);
   die "Error: incorrect number of ANGLE_EQUIL_VALUE!\n"
      if ( $NUMANG != @TEQ-1 );

#  ANGLES_INC_HYDROGEN
   $_ = "zero  ".$hash{"ANGLES_INC_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of ANGLES_INC_HYDROGEN!\n"
      if ( $NTHETH*4 != @_-1 );

   my @ITH; my @JTH; my @KTH; my @ICTH;

   for $iang (1..$NTHETH)
   {
      $ITH[$iang]  = $_[($iang-1)*4+1];
      $JTH[$iang]  = $_[($iang-1)*4+2];
      $KTH[$iang]  = $_[($iang-1)*4+3];
      $ICTH[$iang] = $_[($iang-1)*4+4];
   }

#  ANGLES_WITHOUT_HYDROGEN
   $_ = "zero  ".$hash{"ANGLES_WITHOUT_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of ANGLES_WITHOUT_HYDROGEN!\n"
      if ( $NTHETA*4 != @_-1 );

   my @IT; my @JT; my @KT; my @ICT;

   for $iang (1..$NTHETA)
   {
      $IT[$iang]  = $_[($iang-1)*4+1];
      $JT[$iang]  = $_[($iang-1)*4+2];
      $KT[$iang]  = $_[($iang-1)*4+3];
      $ICT[$iang] = $_[($iang-1)*4+4];
   }

#*********************************************************
#  Read DIHEDRAL_FORCE_CONSTANT, 
#       DIHEDRAL_PERIODICITY,
#       DIHEDRALS_INC_HYDROGEN, 
#       DIHEDRALS_WITHOUT_HYDROGEN
#*********************************************************

   printf "Reading dihedrals ...\n";
   printf "\n";

   my $idih;

#  DIHEDRAL_FORCE_CONSTANT
   $_ = "zero  ".$hash{"DIHEDRAL_FORCE_CONSTANT"};
   my @PK = (split);
   die "Error: incorrect number of DIHEDRAL_FORCE_CONSTANT!\n"
      if ( $NPTRA != @PK-1 );

#  DIHEDRAL_PERIODICITY
   $_ = "zero  ".$hash{"DIHEDRAL_PERIODICITY"};
   my @PN = (split);
   die "Error: incorrect number of DIHEDRAL_PERIODICITY!\n"
      if ( $NPTRA != @PN-1 );

#  DIHEDRAL_PHASE
   $_ = "zero  ".$hash{"DIHEDRAL_PHASE"};
   my @PHASE = (split);
   die "Error: incorrect number of DIHEDRAL_PHASE!\n"
      if ( $NPTRA != @PHASE-1 );

#  DIHEDRALS_INC_HYDROGEN
   $_ = "zero  ".$hash{"DIHEDRALS_INC_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of DIHEDRALS_INC_HYDROGEN!\n"
      if ( $NPHIH*5 != @_-1 );

   my @IPH; my @JPH; my @KPH; my @LPH; my @ICPH;

   for $idih (1..$NPHIH)
   {
      $IPH[$idih]  = $_[($idih-1)*5+1];
      $JPH[$idih]  = $_[($idih-1)*5+2];
      $KPH[$idih]  = $_[($idih-1)*5+3];
      $LPH[$idih]  = $_[($idih-1)*5+4];
      $ICPH[$idih] = $_[($idih-1)*5+5];
   }

#  DIHEDRALS_WITHOUT_HYDROGEN
   $_ = "zero  ".$hash{"DIHEDRALS_WITHOUT_HYDROGEN"};

   @_ = (split);
   die "Error: incorrect number of DIHEDRALS_WITHOUT_HYDROGEN!\n"
      if ( $NPHIA*5 != @_-1 );

   my @IP; my @JP; my @KP; my @LP; my @ICP;

   for $idih (1..$NPHIA)
   {
      $IP[$idih]  = $_[($idih-1)*5+1];
      $JP[$idih]  = $_[($idih-1)*5+2];
      $KP[$idih]  = $_[($idih-1)*5+3];
      $LP[$idih]  = $_[($idih-1)*5+4];
      $ICP[$idih] = $_[($idih-1)*5+5];
   }

#*********************************************************
#  Create gmx top file
#  Write [ defaults ]
#*********************************************************

#  create gmx top file
   open GTOP,">$gmx_top" or die
      "Error: cannot creat gmx top file $gmx_top !\n";

   my $date = `date`;
   chomp $date;
   printf GTOP ";$gmx_top created at $date\n";
   printf GTOP "\n";
   printf GTOP "[ defaults ]\n";
   printf GTOP "%-15s%-15s%-15s%-15s%-15s\n", 
      ";nbfunct", "comb-rule", "gen-pairs", "fudgeLJ", "fudgeQQ";
   printf GTOP "%-15d%-15d%-15s%-15.8f%-15.8f\n", 
      1, 2, "yes", 0.5, 1.0/1.2;
   printf GTOP "\n";
   
#*********************************************************
#  Write [ atomtypes ]
#*********************************************************

   printf "Writing atomtypes ...\n";

#  change amber ACOEF and BCOEF into sigma and epsilon
#  ACOEF = 4.0 * epsilon * sigma**12 ;
#  BCOEF = 4.0 * epsilon * sigma**6 ;
#  sigma = ( ACOEF / BCOEF ) ** (1.0/6.0) ;
#  epsilon = BCOEF**2 / ( 4.0 * ACOEF ) ;
#  remember to change unit: angstrom->nm, kcal/mol->kJ/mol

   my $i; my $itype; 
   my $acoef; my $bcoef;
   my @sigma; my @epsilon;

   for $i (1..$NATOM)
   {
      die "Error: negative value of ICO for atom $i !\n"
         if ( $ICO[$NTYPES*($IAC[$i]-1)+$IAC[$i]] < 0.0 );
      $acoef = $CN1[$ICO[$NTYPES*($IAC[$i]-1)+$IAC[$i]]];
      $bcoef = $CN2[$ICO[$NTYPES*($IAC[$i]-1)+$IAC[$i]]];
      if ( $acoef==0.0 and $bcoef==0.0 )
      {
         $sigma[$i] = 0.0;
         $epsilon[$i] = 0.0;
      }
      elsif ( $acoef>0.0 and $bcoef>0.0 )
      {
         $sigma[$i] = ( $acoef / $bcoef ) ** (1.0/6.0) ;
         $epsilon[$i] = $bcoef**2 / ( 4.0 * $acoef ) ;
         $sigma[$i] *= 0.1; # from angstrom to nm
         $epsilon[$i] *= 4.184; # from kcal/mol to kJ/mol
      }
      else
      {
         die "Error: negative value of CN1 or CN2 for atom $i !\n";
      }
   }

   my %atomtype;
   for $i (1..$NATOM)
   {
      $atomtype{$ISYMBL[$i]} = 0;
   }
   printf GTOP "[ atomtypes ]\n";
   printf GTOP "%-8s%-10s%-10s%-10s%-7s%15s%15s\n",
      ";name", "bond_type", "mass", "charge", "ptype", "sigma", "epsilon";
   for $i (1..$NATOM)
   {
      if ($atomtype{$ISYMBL[$i]}==0)
      {
         for $itype (1..$NTYPES)
         {
            if ( $IAC[$i]==$itype )
            {
               printf GTOP "%-8s%-10s%-10.4f%-10.4f%-7s%15.6e%15.6e\n",
                  $ISYMBL[$i], $ISYMBL[$i], 0.0, 0.0, "A", $sigma[$i], $epsilon[$i];
               last;
            }
          }
         $atomtype{$ISYMBL[$i]} = 1;
      }
   }
   printf GTOP "\n";

#*********************************************************
#  Write [ moleculetype ] and [ atoms ]
#*********************************************************

   my $ires;

   printf GTOP "[ moleculetype ]\n";
   printf GTOP "%-20s%-10s\n", ";name", "nrexcl";
   printf GTOP "%-20s%-10d\n", $options{"-name"}, 3;
   printf GTOP "\n";

   printf "Writing atoms ...\n";
   printf "\n";

#  truncate the charge at predefined precision (10^-6)
   my $qtot = 0.0;
   my $qint = 0;
   my $maxqid = 1;
   for $i (1..$NATOM)
   {
      $CHARGE[$i] /= 18.2223; # change to electron charge unit
#     truncate the charge at 10^-6 precision
      if ( $CHARGE[$i] > 0.0 )
      {
         $CHARGE[$i] = int( $CHARGE[$i] * 10**6 + 0.5 ) / 10**6 ;
      }
      elsif ( $CHARGE[$i] < 0.0 )
      {
         $CHARGE[$i] = int( $CHARGE[$i] * 10**6 - 0.5 ) / 10**6 ;
      }
#     calculate sum of charges
      $qtot += $CHARGE[$i];
#     get the atom id with maximal charge
      $maxqid = $i if ( abs($CHARGE[$i]) > abs($CHARGE[$maxqid]) );
   }
#  remove non-integer charge due to truncated precision
   printf "The system has %f net charges.\n", $qtot;
   if ( $qtot > 0.0 )
   {
      $qint = int( $qtot + 0.5 );
   }
   elsif ( $qtot < 0.0 )
   {
      $qint = int( $qtot - 0.5 );
   }
   if ( $qtot != $qint )
   {
      $CHARGE[$maxqid] -= $qtot-$qint ;
      printf "The charge of %f is removed from atom id %d.\n",
         $qtot-$qint, $maxqid;
   }
   printf "\n";

#  print atoms
   printf GTOP "[ atoms ]\n";
   printf GTOP "%-6s%-6s%-6s%-6s%-6s%-6s%13s%13s ; qtot\n",
      ";nr", "type", "resnr", "res", "atom", "cgnr", "charge", "mass";
   $qtot = 0.0;
   for $ires (1..$NRES)
   {
      printf GTOP ";residue %d, %s\n", $ires, $LBRES[$ires];
      for $i ($IPRES[$ires]..$IPRES[$ires+1]-1)
      {
         $qtot += $CHARGE[$i];
         printf GTOP "%-6d%-6s%-6d%-6s%-6s%-6d%13.6f%13.6f ; qtot %f\n",
            $i, $ISYMBL[$i], $ires, $LBRES[$ires], $IGRAPH[$i], 
            $i, $CHARGE[$i], $AMASS[$i], $qtot;
      }
   }
   printf GTOP "\n";

#*********************************************************
#  Write [ bonds ] 
#*********************************************************

   printf "Writing bonds ...\n";

#  NOTE: the atom numbers in the following arrays that describe 
#  bonds, angles, and dihedrals are coordinate array indexes for 
#  runtime speed. The true atom number equals the absolute value 
#  of the number divided by three, plus one. In the case of the 
#  dihedrals, if the fourth atom is negative, this implies that 
#  the dihedral is an improper. If the third atom is negative, 
#  this implies that the end group interations are to be ignored. 
#  End group interactions are ignored, for example, in dihedrals 
#  of various ring systems (to prevent double counting of 1-4 
#  interactions) and in multiterm dihedrals. 

#  Remember to change the unit of REQ and RK.
#  REQ: angstrom -> nm
#  RK:  kcal/mol/A^2 -> kJ/mol/nm^2
#  Note: in gmx RK should be doubled due to the 1/2 prefactor.

   printf GTOP "[ bonds ]\n";
   printf GTOP "%-8s%-8s%-8s%15s%15s\n", ";ai", "aj", "funct", "r", "k";
   for $ibond (1..$NBONH)
   {
      printf GTOP "%-8d%-8d%-8d%15.6e%15.6e ; %s-%s\n",
         abs($IBH[$ibond])/3+1, 
         abs($JBH[$ibond])/3+1, 
         1,
         $REQ[$ICBH[$ibond]] * 0.1, 
         $RK[$ICBH[$ibond]] * 4.184 * 100.0 * 2.0,
         $IGRAPH[abs($IBH[$ibond])/3+1],
         $IGRAPH[abs($JBH[$ibond])/3+1];
   }
   for $ibond (1..$NBONA)
   {
      printf GTOP "%-8d%-8d%-8d%15.6e%15.6e ; %s-%s\n",
         abs($IB[$ibond])/3+1,
         abs($JB[$ibond])/3+1,
         1,
         $REQ[$ICB[$ibond]] * 0.1,
         $RK[$ICB[$ibond]] * 4.184 * 100.0 * 2.0,
         $IGRAPH[abs($IB[$ibond])/3+1],
         $IGRAPH[abs($JB[$ibond])/3+1];
   }
   printf GTOP "\n";

#*********************************************************
#  Write [ pairs ] 
#*********************************************************

   printf "Writing pairs ...\n";

#  NOTE: the atom numbers in the following arrays that describe 
#  bonds, angles, and dihedrals are coordinate array indexes for 
#  runtime speed. The true atom number equals the absolute value 
#  of the number divided by three, plus one. In the case of the 
#  dihedrals, if the fourth atom is negative, this implies that 
#  the dihedral is an improper. If the third atom is negative, 
#  this implies that the end group interations are to be ignored. 
#  End group interactions are ignored, for example, in dihedrals 
#  of various ring systems (to prevent double counting of 1-4 
#  interactions) and in multiterm dihedrals. 

   if ( $use_glycam eq "yes" )
   {

      printf GTOP "[ pairs ]\n";
      printf GTOP "%-6s%-6s%-6s%13s%13s%13s%15s%15s\n", 
         ";ai", "aj", "funct", "fudgeQQ", "qi", "qj", "sigma", "epsilon";
 
      for $idih (1..$NPHIH)
      {
#        skip when the third atom is negative
         if ( $KPH[$idih] > 0.0 )
         {
            printf GTOP "%-6d%-6d%-6d%13.6f%13.6f%13.6f%15.6e%15.6e ; %s-%s\n",
               abs($IPH[$idih])/3+1, abs($LPH[$idih])/3+1, 2,
               1.0, $CHARGE[abs($IPH[$idih])/3+1], $CHARGE[abs($LPH[$idih])/3+1],
               0.5 * ( $sigma[abs($IPH[$idih])/3+1] + $sigma[abs($LPH[$idih])/3+1] ),
               sqrt( $epsilon[abs($IPH[$idih])/3+1] * $epsilon[abs($LPH[$idih])/3+1] ),
               $IGRAPH[abs($IPH[$idih])/3+1],
               $IGRAPH[abs($LPH[$idih])/3+1];
         }
      }
 
      for $idih (1..$NPHIA)
      {
#        skip when the third atom is negative
         if ( $KP[$idih] > 0.0 )
         {
            printf GTOP "%-6d%-6d%-6d%13.6f%13.6f%13.6f%15.6e%15.6e ; %s-%s\n",
               abs($IP[$idih])/3+1, abs($LP[$idih])/3+1, 2,
               1.0, $CHARGE[abs($IP[$idih])/3+1], $CHARGE[abs($LP[$idih])/3+1],
               0.5 * ( $sigma[abs($IP[$idih])/3+1] + $sigma[abs($LP[$idih])/3+1] ),
               sqrt( $epsilon[abs($IP[$idih])/3+1] * $epsilon[abs($LP[$idih])/3+1] ),
               $IGRAPH[abs($IP[$idih])/3+1],
               $IGRAPH[abs($LP[$idih])/3+1];
         }
      }
 
      printf GTOP "\n";

   }
   else
   {

      printf GTOP "[ pairs ]\n";
      printf GTOP "%-8s%-8s%-8s\n", ";ai", "aj", "funct";
 
      for $idih (1..$NPHIH)
      {
#        skip when the third atom is negative
         if ( $KPH[$idih] > 0.0 )
         {
            printf GTOP "%-8d%-8d%-8d ; %s-%s\n",
               abs($IPH[$idih])/3+1, abs($LPH[$idih])/3+1, 1,
               $IGRAPH[abs($IPH[$idih])/3+1],
               $IGRAPH[abs($LPH[$idih])/3+1];
         }
      }
 
      for $idih (1..$NPHIA)
      {
#        skip when the third atom is negative
         if ( $KP[$idih] > 0.0 )
         {
            printf GTOP "%-8d%-8d%-8d ; %s-%s\n",
               abs($IP[$idih])/3+1, abs($LP[$idih])/3+1, 1,
               $IGRAPH[abs($IP[$idih])/3+1],
               $IGRAPH[abs($LP[$idih])/3+1];
         }
      }
 
      printf GTOP "\n";

   }

#*********************************************************
#  Write [ angles ] 
#*********************************************************

   printf "Writing angles ...\n";

#  NOTE: the atom numbers in the following arrays that describe 
#  bonds, angles, and dihedrals are coordinate array indexes for 
#  runtime speed. The true atom number equals the absolute value 
#  of the number divided by three, plus one. In the case of the 
#  dihedrals, if the fourth atom is negative, this implies that 
#  the dihedral is an improper. If the third atom is negative, 
#  this implies that the end group interations are to be ignored. 
#  End group interactions are ignored, for example, in dihedrals 
#  of various ring systems (to prevent double counting of 1-4 
#  interactions) and in multiterm dihedrals. 

#  Remember to change the unit of TEQ and TK.
#  TEQ: radian -> degree
#  TK:  kcal/mol/rad^2 -> kJ/mol/rad^2
#  Note: in gmx TK should be doubled due to the 1/2 prefactor.

   printf GTOP "[ angles ]\n";
   printf GTOP "%-7s%-7s%-7s%-7s%15s%15s\n", 
      ";ai", "aj", "ak", "funct", "theta", "cth";
   for $iang (1..$NTHETH)
   {
      printf GTOP "%-7d%-7d%-7d%-7d%15.6e%15.6e ; %s-%s-%s\n",
         abs($ITH[$iang])/3+1,
         abs($JTH[$iang])/3+1,
         abs($KTH[$iang])/3+1,
         1,
         $TEQ[$ICTH[$iang]] * 180.0 / $pi,
         $TK[$ICTH[$iang]] * 4.184 * 2.0,
         $IGRAPH[abs($ITH[$iang])/3+1],
         $IGRAPH[abs($JTH[$iang])/3+1],
         $IGRAPH[abs($KTH[$iang])/3+1];
   }
   for $iang (1..$NTHETA)
   {
      printf GTOP "%-7d%-7d%-7d%-7d%15.6e%15.6e ; %s-%s-%s\n",
         abs($IT[$iang])/3+1,
         abs($JT[$iang])/3+1,
         abs($KT[$iang])/3+1,
         1,
         $TEQ[$ICT[$iang]] * 180.0 / $pi,
         $TK[$ICT[$iang]] * 4.184 * 2.0,
         $IGRAPH[abs($IT[$iang])/3+1],
         $IGRAPH[abs($JT[$iang])/3+1],
         $IGRAPH[abs($KT[$iang])/3+1];
   }
   printf GTOP "\n";

#*********************************************************
#  Write [ dihedrals ] 
#*********************************************************

   printf "Writing dihedrals ...\n";
   printf "\n";

#  NOTE: the atom numbers in the following arrays that describe 
#  bonds, angles, and dihedrals are coordinate array indexes for 
#  runtime speed. The true atom number equals the absolute value 
#  of the number divided by three, plus one. In the case of the 
#  dihedrals, if the fourth atom is negative, this implies that 
#  the dihedral is an improper. If the third atom is negative, 
#  this implies that the end group interations are to be ignored. 
#  End group interactions are ignored, for example, in dihedrals 
#  of various ring systems (to prevent double counting of 1-4 
#  interactions) and in multiterm dihedrals. 

#  Torsion in AMBER force field
#  E_Torsion = V_n/2 * [ 1 + cos( n phi - phase ) ]
#  Note: The PK value is equal to V_n/2 in AMBER FF literature

#  Modified 2012-06-29 by Xin Li

   my $ip;
   my ($p_s, $k_p, $p_n);

   for $ip (1..$NPTRA)
   {
#     convert kcal/mol to kJ/mol
      $PK[$ip] *= 4.184;

#     convert rad to degree
      $PHASE[$ip] = $PHASE[$ip] / $pi * 180.0;
   }

#  write dihedrals

   printf GTOP "[ dihedrals ]\n";
   printf GTOP "%-6s%-6s%-6s%-6s%-6s%10s%13s%6s\n", 
      ";ai", "aj", "ak", "al", "funct", 
      "phi_s", "k_phi", "mult";

   my $funct;

#  dihedrals containing H atoms

   for $idih (1..$NPHIH)
   {
#     phase (phi_s), force constant (k_phi), and multiplicity
      $p_s = $PHASE[$ICPH[$idih]];
      $k_p = $PK[$ICPH[$idih]];
      $p_n = $PN[$ICPH[$idih]];

#     distinguish proper and improper
      if ( $LPH[$idih] >= 0.0 )
      {
         $funct = 9;
      }
      else
      {
         $funct = 4;
      }

#     print atoms and parameters (note funct of 9 and 4 are used)
      printf GTOP "%-6d%-6d%-6d%-6d%-6d%10.3f%13.6f%6d ; %s-%s-%s-%s\n",
         abs($IPH[$idih])/3+1, abs($JPH[$idih])/3+1,
         abs($KPH[$idih])/3+1, abs($LPH[$idih])/3+1,
         $funct, $p_s, $k_p, $p_n,
         $IGRAPH[abs($IPH[$idih])/3+1],
         $IGRAPH[abs($JPH[$idih])/3+1],
         $IGRAPH[abs($KPH[$idih])/3+1],
         $IGRAPH[abs($LPH[$idih])/3+1];
   }

#  dihedrals without H atoms

   for $idih (1..$NPHIA)
   {
#     phase (phi_s), force constant (k_phi), and multiplicity
      $p_s = $PHASE[$ICP[$idih]];
      $k_p = $PK[$ICP[$idih]];
      $p_n = $PN[$ICP[$idih]];

#     distinguish proper and improper
      if ( $LP[$idih] >= 0.0 )
      {
         $funct = 9;
      }
      else
      {
         $funct = 4;
      }

#     print atoms and parameters (note funct of 9 and 4 are used)
      printf GTOP "%-6d%-6d%-6d%-6d%-6d%10.3f%13.6f%6d ; %s-%s-%s-%s\n",
         abs($IP[$idih])/3+1, abs($JP[$idih])/3+1,
         abs($KP[$idih])/3+1, abs($LP[$idih])/3+1,
         $funct, $p_s, $k_p, $p_n,
         $IGRAPH[abs($IP[$idih])/3+1],
         $IGRAPH[abs($JP[$idih])/3+1],
         $IGRAPH[abs($KP[$idih])/3+1],
         $IGRAPH[abs($LP[$idih])/3+1];
   }

   printf GTOP "\n";

#  End of modification 2012-06-29

#*********************************************************
#  Write [ system ] and [ molecules ]
#*********************************************************

   printf GTOP "[ system ]\n";
   printf GTOP "%s\n", $options{"-name"};
   printf GTOP "\n";

   printf GTOP "[ molecules ]\n";
   printf GTOP "%-12s%-12s\n", ";Compound", "nmols";
   printf GTOP "%-12s%-12d\n", $options{"-name"}, 1;
   printf GTOP "\n";

#*********************************************************
#  Close gmx top file
#*********************************************************

#  close gmx top file
   close GTOP;
   printf "GMX top written to $gmx_top \n";
   printf "\n";

#*********************************************************
#  Read amber crd file
#  Write gmx gro file
#*********************************************************

#  read amber crd file
   $string = "";
   open ACRD,"$amb_crd" or die
      "Error: cannot open amber crd file $amb_crd !\n";
#  read title
   $_ = <ACRD>;
#  read number of atoms
   $_ = <ACRD>;
   if ( $NATOM != (split)[0] )
   {
      die "Error: inconsistent NATOM in $amb_top and $amb_crd !\n";
   }
#  read coordinates
   while ( <ACRD> )
   {
      chomp;
      $string .= " ".$_." ";
   }
#  close amber crd file
   close ACRD;
#  split coordinates
   $string = "zero  ".$string;
   my @xyz = split " ",$string;
   if ( $NATOM*3 != @xyz-1 )
   {
      die "Error: incorrect number of coordinates in $amb_crd !\n";
   }

#  create gmx gro file
   open GGRO,">$gmx_gro" or die
      "Error: cannot creat gmx gro file $gmx_gro !\n";

#  write title of gmx gro
   $date = `date`;
   chomp $date;
   printf GGRO "$gmx_gro created at $date\n";

#  write gmx gro file
#  remember to convert angstrom to nm
   printf GGRO "%5d\n", $NATOM;
   for $ires (1..$NRES)
   {
      for $i ($IPRES[$ires]..$IPRES[$ires+1]-1)
      {
         printf GGRO "%5d%5s%5s%5d%9.4f%9.4f%9.4f\n",
            $ires, $LBRES[$ires], $IGRAPH[$i], $i, 
            $xyz[($i-1)*3+1] * 0.1 + 5.0,
            $xyz[($i-1)*3+2] * 0.1 + 5.0,
            $xyz[($i-1)*3+3] * 0.1 + 5.0;
      }
   }
#  in a big box
   printf GGRO "%10.5f%10.5f%10.5f\n", 10.0, 10.0, 10.0;
#  close gmx gro file
   close GGRO;
   printf "GMX gro written to $gmx_gro \n";
   printf "\n";

#  All done!
   printf "All done!\n";

#*********************************************************
#  End of program
#*********************************************************
