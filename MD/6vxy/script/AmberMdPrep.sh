#!/bin/bash

# AmberMdPrep.sh
# Wrapper script for preparing explicitly solvated systems for MD with Amber.
# Daniel R. Roe
# NIH/NHLBI
# 2020-08-07

VERSION='0.5 (beta)'
MPIRUN=`which mpirun`
CPPTRAJ=`which cpptraj`

TEST=0
INPUT=''
RUNTYPE=''
EVALTYPE='cpptraj'
TOP=''
CRD=''
REF=''
S=0
TYPE=''
HEAVYMASK=''
BACKBONEMASK=''
ADDITIONALMASK='' # Additional mask to use during steps 1-8
PRODUCTIONMASK='' # Additional mask to use steps 9 and above
PRODUCTIONWT=''   # Restraint weight to use for production mask
PRODUCTIONREF=''  # Reference structure to use for production mask
ADDEDRES=''       # Residue names that if present will be added to heavy masks
MASTERCUT=''      # If set, override default cutoff
TEMP0=''
TEMPI=''
NPROCS=''  # Number of processes to use for CPU jobs
OVERWRITE=0
CHARMMWATERFLAG=0 # Set when charmm water present
NPROTEIN=0 # Number of protein residues
NDNA=0     # Number of DNA residues
NRNA=0     # Number of RNA residues
NLIPID=0   # Number of lipid residues
NCARBO=0   # Number of carbohydrate residues
NUNKNOWN=0 # Number of unknown residues
NCHARMMWATER=0
NWATER=0
LIPIDRESNAMES=''       # Comma-separated list of unique lipid residue names
THERMOTYPE='langevin'   # Thermostat type
BAROTYPE='montecarlo'   # Barostat type
FINALTHERMO='langevin'  # Thermostat for final density eq.
FINALBARO='montecarlo'  # Barostat for final density eq.
NTPFLAG=1 # 1 for isotropic scaling, 2 for anisotropic
FLEXIBLEWAT='' # Set to flexiblewat to ensure every MD step is OK for flexible waters
PRODUCTION_DT='0.002' # Default time step for 'production' steps (9 and final)
STATUSFILE=''

# ------------------------------------------------------------------------------
# DetectSystemType <topology file>
# Use cpptraj to print out residues from a topology file. Determine the
# type of each residue and keep track.
DetectSystemType() {
  if [ -z "$CPPTRAJ" ] ; then
    echo "Error. Topology detection relies on cpptraj, which is not present."
    exit 1
  fi
  if [ -z "$1" ] ; then
    echo "Error. Topology is blank."
    exit 1
  fi
  TMPLIPID='tmp.lipidResNames'
  if [ -f "$TMPLIPID" ] ; then
    rm $TMPLIPID
  fi
  TMPUNKNOWN='tmp.unknownResNames'
  if [ -f "$TMPUNKNOWN" ] ; then
    rm $TMPUNKNOWN
  fi
  systemNumbers=`$CPPTRAJ -p $1 --resmask \* | awk -v tmplipid="$TMPLIPID" -v tmpunknown="$TMPUNKNOWN" 'BEGIN{
    nprotein = 0;
    ndna = 0;
    nrna = 0;
    nlipid = 0;
    nunknown = 0;
    ncharmmwater = 0;
    nwater = 0;
    ncarbo = 0;
  }{
    if ($2 != "Name") {
      if ($2 == "ACE" ||
          $2 == "ALA" ||
          $2 == "ARG" ||
          $2 == "ASH" ||
          $2 == "AS4" ||
          $2 == "ASN" ||
          $2 == "ASP" ||
          $2 == "CALA" ||
          $2 == "CARG" ||
          $2 == "CASN" ||
          $2 == "CASP" ||
          $2 == "CCYS" ||
          $2 == "CCYX" ||
          $2 == "CGLN" ||
          $2 == "CGLU" ||
          $2 == "CGLY" ||
          $2 == "CHID" ||
          $2 == "CHIE" ||
          $2 == "CHIP" ||
          $2 == "CHIS" ||
          $2 == "CHYP" ||
          $2 == "CILE" ||
          $2 == "CLEU" ||
          $2 == "CLYS" ||
          $2 == "CMET" ||
          $2 == "CPHE" ||
          $2 == "CPRO" ||
          $2 == "CSER" ||
          $2 == "CTHR" ||
          $2 == "CTRP" ||
          $2 == "CTYR" ||
          $2 == "CVAL" ||
          $2 == "CYM" ||
          $2 == "CYS" ||
          $2 == "CYX" ||
          $2 == "GLH" ||
          $2 == "GL4" ||
          $2 == "GLN" ||
          $2 == "GLU" ||
          $2 == "GLY" ||
          $2 == "HID" ||
          $2 == "HIE" ||
          $2 == "HIP" ||
          $2 == "HIS" ||
          $2 == "HYP" ||
          $2 == "ILE" ||
          $2 == "LEU" ||
          $2 == "LYN" ||
          $2 == "LYS" ||
          $2 == "MET" ||
          $2 == "NALA" ||
          $2 == "NARG" ||
          $2 == "NASN" ||
          $2 == "NASP" ||
          $2 == "NCYS" ||
          $2 == "NCYX" ||
          $2 == "NGLN" ||
          $2 == "NGLU" ||
          $2 == "NGLY" ||
          $2 == "NHE" ||
          $2 == "NHID" ||
          $2 == "NHIE" ||
          $2 == "NHIP" ||
          $2 == "NHIS" ||
          $2 == "NILE" ||
          $2 == "NLEU" ||
          $2 == "NLYS" ||
          $2 == "NME" ||
          $2 == "NMET" ||
          $2 == "NPHE" ||
          $2 == "NPRO" ||
          $2 == "NSER" ||
          $2 == "NTHR" ||
          $2 == "NTRP" ||
          $2 == "NTYR" ||
          $2 == "NVAL" ||
          $2 == "PHE" ||
          $2 == "PRO" ||
          $2 == "SER" ||
          $2 == "THR" ||
          $2 == "TRP" ||
          $2 == "TYR" ||
          $2 == "VAL")
        nprotein++;
      else if ($2 == "DA" ||
               $2 == "DA3" ||
               $2 == "DA5" ||
               $2 == "DAN" ||
               $2 == "DC" ||
               $2 == "DC3" ||
               $2 == "DC5" ||
               $2 == "DCN" ||
               $2 == "DG" ||
               $2 == "DG3" ||
               $2 == "DG5" ||
               $2 == "DGN" ||
               $2 == "DT" ||
               $2 == "DT3" ||
               $2 == "DT5" ||
               $2 == "DTN")
        ndna++;
      else if ($2 == "A" ||
               $2 == "A3" ||
               $2 == "A5" ||
               $2 == "AMP" ||
               $2 == "AN" ||
               $2 == "C" ||
               $2 == "C3" ||
               $2 == "C5" ||
               $2 == "CMP" ||
               $2 == "CN" ||
               $2 == "G" ||
               $2 == "G3" ||
               $2 == "G5" ||
               $2 == "GMP" ||
               $2 == "GN" ||
               $2 == "OHE" ||
               $2 == "U" ||
               $2 == "U3" ||
               $2 == "U5" ||
               $2 == "UMP" ||
               $2 == "UN")
        nrna++;
      else if ($2 == "POPE" ||
               $2 == "DOPC" ||
               $2 == "AR" ||
               $2 == "CHL" ||
               $2 == "DHA" ||
               $2 == "LAL" ||
               $2 == "MY" ||
               $2 == "OL" ||
               $2 == "PA" ||
               $2 == "PC" ||
               $2 == "PE" ||
               $2 == "PGR" ||
               $2 == "PH-" ||
               $2 == "PS" ||
               $2 == "ST")
      {
        nlipid++;
        print $2 >> tmplipid;
      } else if ($2 == "0GB" ||
                 $2 == "4GB" ||
                 $2 == "0YA" ||
                 $2 == "4YA" ||
                 $2 == "0fA" ||
                 $2 == "0YB" ||
                 $2 == "2MA" ||
                 $2 == "4YB" ||
                 $2 == "NLN" ||
                 $2 == "UYB" ||
                 $2 == "VMB" ||
                 $2 == "0SA" ||
                 $2 == "6LB" ||
                 $2 == "ROH")
      {
        ncarbo++;
      } else if ($2 == "TIP3") {
        ncharmmwater++;
      } else if ($2 == "WAT") {
        nwater++;
      } else {
        nunknown++;
        print $2 >> tmpunknown;
      }
    }
  }END{
    printf("%i %i %i %i %i %i %i %i\n", nprotein, ndna, nrna, nlipid, nunknown, ncharmmwater, nwater, ncarbo);
  }'`
  #echo "DEBUG: $systemNumbers"
  if [ $? -ne 0 -o -z "$systemNumbers" ] ; then
    echo "System detection failed."
    exit 1
  fi
  NPROTEIN=`echo $systemNumbers | awk '{print $1;}'`
  NDNA=`echo $systemNumbers | awk '{print $2;}'`
  NRNA=`echo $systemNumbers | awk '{print $3;}'`
  NLIPID=`echo $systemNumbers | awk '{print $4;}'`
  NUNKNOWN=`echo $systemNumbers | awk '{print $5;}'`
  NCHARMMWATER=`echo $systemNumbers | awk '{print $6;}'`
  NWATER=`echo $systemNumbers | awk '{print $7;}'`
  NCARBO=`echo $systemNumbers | awk '{print $8;}'`
  if [ $NCHARMMWATER -gt 0 ] ; then
    if [ $NWATER -gt 0 ] ; then
      echo "Error: Charmm water and regular water present."
      exit 1
    fi
    hasCharmmWater=1
  else
    hasCharmmWater=0
  fi
  NWATER=`echo $systemNumbers | awk '{print $7;}'`
  ((NTOTALWATER = $NCHARMMWATER + $NWATER))
  printf "  %i protein, %i dna, %i rna, %i lipid, %i carbohydrate, %i water, %i other\n" $NPROTEIN $NDNA $NRNA $NLIPID $NCARBO $NTOTALWATER $NUNKNOWN
  if [ $NLIPID -gt 0 ] ; then
    for lres in `sort $TMPLIPID | uniq` ; do
      if [ -z "$LIPIDRESNAMES" ] ; then
        LIPIDRESNAMES=$lres
      else
        LIPIDRESNAMES="$LIPIDRESNAMES,$lres"
      fi
    done
    echo "  Lipid residue names: $LIPIDRESNAMES"
    rm $TMPLIPID
  fi
  amask=''
  if [ $NUNKNOWN -gt 0 ] ; then
    for ures in `sort $TMPUNKNOWN | uniq` ; do 
      if [ -z "$UNKRESNAMES" ] ; then
        UNKRESNAMES=$ures
      else
        UNKRESNAMES="$UNKRESNAMES,$ures"
      fi
      if [ ! -z "$ADDEDRES" ] ; then
        # Does ures match an ADDEDRES?
        for ares in $ADDEDRES ; do
          if [ "$ures" = "$ares" ] ; then
            if [ -z "$amask" ] ; then
              amask=":$ares"
            else
              amask="$amask,$ares"
            fi
          fi
        done
      fi
    done
    echo "  Unknown residues names: $UNKRESNAMES"
  fi
  if [ ! -z "$amask" ] ; then
    amask="$amask&!@H="
    echo "  Detected additional res mask: $amask"
    if [ ! -z "$ADDITIONALMASK" ] ; then
      echo "Error: Additional mask already set to $ADDITIONALMASK"
      exit 1
    fi
    ADDITIONALMASK=$amask
  fi
  if [ $hasCharmmWater -eq 1 -a $CHARMMWATERFLAG -eq 0 ] ; then
    echo "  TIP3 residue detected - assuming CHARMM water present."
    CHARMMWATERFLAG=1
  fi
}

# ------------------------------------------------------------------------------
KeyHelp() {
  echo "Recognized input file keys"
  echo "  'ntmin'"
  echo "  'maxcyc'"
  echo "  'ncyc'"
  echo "  'restraintmask'" 
  echo "  'restraint_wt'"
  echo "  'irest'" 
  echo "  'nstlim'"
  echo "  'ntb'"
  echo "  'ntc'"
  echo "  'cut'"
  echo "  'tempi'"
  echo "  'tautp'"
  echo "  'taup'"
  echo "  'mcbarint'"
  echo "  'gamma_ln'"
  echo "  'dt'"
  echo "  'nscm'"
  echo "  'ntwx'"
  echo "  'ntpr'"
  echo "  'ntwr'"
  echo "  'previousref'"
  echo "  'heavyrst'"
  echo "  'bbrst'"
  echo "  'thermo'"
  echo "  'baro'"
  echo "  'flexiblewat' : shorthand for ntc=ntf=1, jfastw=4"
}

# ------------------------------------------------------------------------------
Help() {
  echo "Command line options"
  echo " Required:"
  echo "  -p <file>            : Topology (required)"
  echo "  -c <file>            : Coordinates (required)"
  echo "  --temp <temp>        : Temperature (required)"
  echo " Optional:"
  echo "  -i <file>            : File containing custom minimization/MD steps."
  echo "                         Lines starting with 'min' are minimization steps."
  echo "                         Lines starting with 'md' are MD steps."
  echo "                         Use '--keyhelp' for recognized keys to set up steps."
  echo "  --thermo <type>      : Thermostat: berendsen, langevin (default)"
  echo "  --baro <type>        : Barostat: berendsen, montecarlo (default)"
  echo "  --finalthermo <type> : Final stage thermostat: berendsen, langevin (default)"
  echo "  --finalbaro <type>   : Final stage barostat: berendsen, montecarlo (default)"
  echo "  --nsolute <#>        : Number of solute residues."
  echo "  --type <type>        : Residues type {protein* | nucleic}; determines backbone mask."
  echo "  --mask <mask>        : Additional mask to use for restraints during steps 1-8."
  echo "  --ares <name>        : Residue name to add to heavy atom masks if present."
  echo "  --pmask <mask>       : Restraint mask to use during \"production\" (steps 9 and above)."
  echo "  --pwt <weight>       : Restraint weight to use for '--pmask'; required if '--pmask' specified."
  echo "  --pref <file>        : Optional reference structure to use for '--pmask'."
  echo "  --charmmwater        : If specified assume CHARMM water (i.e. 'TIP3')."
  echo "  --flexiblewat        : If specified assume solvent is flexible (e.g. SPC/FW)."
  echo "  --cutoff <cut>       : If specified, override default cutoffs with <cut>." 
  echo "  --test               : Test only. Do not run."
  echo "  --norestart          : Do standard Eq with no restarts."
  #echo "  --evaltype <type>    : <type = {script|cpptraj} Evaluate with EvalEquilibration.sh or cpptraj."
  #echo "                         If this is not specified, final density eq. (step 10) will be skipped."
  echo "  --skipfinaleq        : If specified, skip final eq. (step 10)."
  echo "  --nprocs <#>         : Number of CPU processes to use (default 4)."
  echo "  -O                   : Overwrite existing files, otherwise skip."
  echo "  --keyhelp            : Print help for recognized input file keys."
  echo "  --statusfile <file>  : Status file for final density eq."
  echo "Environment vars"
  echo "  PROG_MIN : Command for minimization steps."
  echo "  PROG_MD  : Command for MD steps."
}
# ==============================================================================

echo "AmberMdPrep.sh Version $VERSION"
# Parse command line options
while [ ! -z "$1" ] ; do
  case "$1" in
    '-p'            ) shift ; TOP=$1 ;;
    '-c'            ) shift ; CRD=$1 ;;
    '-i'            ) shift ; INPUT=$1 ;;
    '--nsolute'     ) shift ; S=$1 ;;
    '--type'        ) shift ; TYPE="$TYPE $1" ;;
    '--temp'        ) shift ; TEMP0=$1 ;;
    '--thermo'      ) shift ; THERMOTYPE=$1 ;;
    '--baro'        ) shift ; BAROTYPE=$1 ;;
    '--finalthermo' ) shift ; FINALTHERMO=$1 ;;
    '--finalbaro'   ) shift ; FINALBARO=$1 ;;
    #'--evaltype'    ) shift ; EVALTYPE=$1 ;;
    '--skipfinaleq' ) EVALTYPE='' ;;
    '--mask'        ) shift ; ADDITIONALMASK=$1 ;;
    '--ares'        ) shift ; ADDEDRES="$ADDEDRES $1" ;;
    '--pmask'       ) shift ; PRODUCTIONMASK=$1 ;;
    '--pwt'         ) shift ; PRODUCTIONWT=$1 ;;
    '--pref'        ) shift ; PRODUCTIONREF=$1 ;;
    '--cutoff'      ) shift ; MASTERCUT=$1 ;;
    '--charmmwater' ) CHARMMWATERFLAG=1 ;;
    '-O'            ) OVERWRITE=1 ;;
    '-h' | '--help' ) Help ; exit 0 ;;
    '--keyhelp'     ) KeyHelp ; exit 0 ;;
    '--test'        ) TEST=1 ;;
    '--norestart'   ) RUNTYPE='norestart' ;;
    '--nprocs'      ) shift ; NPROCS=$1 ;;
    '--statusfile'  ) shift ; STATUSFILE=$1 ;;
    '--flexiblewat' ) FLEXIBLEWAT='flexiblewat' ; PRODUCTION_DT='0.001' ;;
    *               ) echo "Unrecognized command line option: $1" >> /dev/stderr ; exit 1 ;;
  esac
  shift
done
REF=$CRD

echo "  TOP            : $TOP"
echo "  CRD            : $CRD"
if [ ! -f "$TOP" ] ; then
  echo "Topology $TOP not found." >> /dev/stderr
  exit 1
fi
if [ ! -f "$CRD" ] ; then
  echo "Coords $CRD not found." >> /dev/stderr
  exit 1
fi

if [ -z "$TEMP0" ] ; then
  echo "Specify temperature." >> /dev/stderr
  exit 1
fi

if [ ! -z "$PRODUCTIONMASK" -a -z "$PRODUCTIONWT" ] ; then
  echo "'--pwt <weight>' must be specified if '--pmask' specified." >> /dev/stderr
  exit 1
fi

if [ -z "$STATUSFILE" ] ; then
  STATUSFILE=/dev/stdout
else
  echo "Status file is $STATUSFILE"
fi

#if [ "$EVALTYPE" != 'script' -a "$EVALTYPE" != 'cpptraj' ] ; then
#  echo "Error: Please provide evaluation type. '--evaltype {script|cpptraj}'" >> /dev/stderr
#  exit 1
#fi

# Determine program parameters if needed
if [ -z "$NPROCS" ] ; then
  NPROCS=4
fi
if [ -z "$PROG_MIN" ] ; then
  # Program to run minimization steps
  if [ ! -z "`which pmemd.cuda_DPFP`" ] ; then
    PROG_MIN='pmemd.cuda_DPFP'
  elif [ ! -z "$MPIRUN" -a ! -z "`which sander.MPI`" ] ; then 
    PROG_MIN="mpirun -n $NPROCS sander.MPI"
  else
    PROG_MIN='sander'
  fi
fi
if [ -z "$PROG_MD" ] ; then
  # Program to run MD steps
  if [ ! -z "`which pmemd.cuda`" ] ; then
    PROG_MD='pmemd.cuda'
  elif [ ! -z "$MPIRUN" -a ! -z "`which pmemd.MPI`" ] ; then
    PROG_MD="mpirun -n $NPROCS pmemd.MPI"
  else
    PROG_MD='pmemd'
  fi
fi

# Determine system type if needed
if [ $S -lt 1 -o -z "$TYPE" ] ; then
  DetectSystemType $TOP
  if [ $NUNKNOWN -gt 0 ] ; then
    echo "Warning: Unknown residues detected; will be ignored for restraints."
  fi
  if [ $S -lt 1 ] ; then
    ((S = $NPROTEIN + $NDNA + $NRNA + $NLIPID + $NCARBO))
  fi
  if [ -z "$TYPE" ] ; then
    if [ $NPROTEIN -gt 0 ] ; then
      TYPE="$TYPE protein"
    fi
    if [ $NDNA -gt 0 -o $NRNA -gt 0 ] ; then
      TYPE="$TYPE nucleic"
    fi
    if [ $NLIPID -gt 0 ] ; then
      TYPE="$TYPE lipid"
    fi
    if [ $NCARBO -gt 0 ] ; then
      TYPE="$TYPE carbo"
    fi
  fi
  echo "  Detected types : $TYPE"
fi

# Determine box type.
BoxLine=`echo 'list trajin' | $CPPTRAJ -p $TOP -y $CRD | grep "Orthogonal box"`
#echo "$BoxLine"
if [ ! -z "$BoxLine" ] ; then
  hasOrthoBox=1
  echo "  Orthogonal box detected."
  # If lipids also present assume we need anisotropic pressure scaling
  if [ $NLIPID -gt 0 ] ; then
    echo "  Orthogonal box and lipids present; using anisotropic pressure scaling."
    NTPFLAG=2
  fi
else
  hasOrthoBox=0
fi

if [ $S -lt 1 -a -z "$ADDITIONALMASK" ] ; then
  echo "Specify number of solute residues or additional mask." >> /dev/stderr
  exit 1
fi

AssignMask() {
  if [ -z "$atommask" ] ; then
    atommask=$1
  else
    atommask=$atommask",$1"
  fi
}

# Set up solute mask
if [ $S -gt 0 ] ; then
  HEAVYMASK=":1-$S&!@H="
  atommask=''
  for rtype in $TYPE ; do
    if [ "$rtype" = 'protein' ] ; then
      AssignMask "H,N,CA,HA,C,O"
    elif [ "$rtype" = 'nucleic' ] ; then
      AssignMask "P,O5',C5',C4',C3',O3'"
    elif [ "$rtype" != 'lipid' -a "$rtype" != 'carbo' ] ; then
      echo "Unrecognized type: $rtype"
      exit 1
    fi
  done
  if [ -z "$atommask" ] ; then
    # No types. Use HEAVYMASK.
    BACKBONEMASK=$HEAVYMASK
  else
    BACKBONEMASK=":1-$S@$atommask"
    if [ ! -z "$LIPIDRESNAMES" ] ; then
      BACKBONEMASK=$BACKBONEMASK"|:$LIPIDRESNAMES&!@H="
    fi
  fi
fi

echo "  NUM SOLUTE RES : $S"
echo "  HEAVY MASK     : $HEAVYMASK"
echo "  BACKBONE MASK  : $BACKBONEMASK"
if [ ! -z "$ADDITIONALMASK" ] ; then
  echo "  ADDITIONALMASK : $ADDITIONALMASK"
fi
if [ ! -z "$PRODUCTIONMASK" ] ; then
  echo "  PRODUCTIONMASK : $PRODUCTIONMASK"
  echo "  PRODUCTIONWT   : $PRODUCTIONWT"
  if [ ! -z "$PRODUCTIONREF" ] ; then
    echo "  PRODUCTIONREF  : $PRODUCTIONREF"
  fi
fi
if [ ! -z "$ADDEDRES" ] ; then
  printf "  ADD. MASK RES. :"
  for ares in $ADDEDRES ; do
    printf " $ares"
  done
  printf "\n"
fi
echo "  TEMPERATURE    : $TEMP0"
echo "  OVERWRITE      : $OVERWRITE"
echo "  MD COMMAND     : $PROG_MD"
echo "  MIN COMMAND    : $PROG_MIN"
echo "  NPROCS         : $NPROCS"

# ==============================================================================
# Amber options
MDIN=''
RST=''
NTMIN=''         # 1 = SD+CG, 2 = steepest descent, 3 = xmin
MAXCYC=''        # Number of minimization cycles
NCYC=''          # For ntmin 1, Min will be switched from SD to CG after NCYC cycles.
RESTRAINTMASK='' # Restraint mask
RESTRAINT_WT=''  # Restraint weight
IREST=0
NSTLIM=0
NTB=0
TAUTP=''
TAUP==''
MCBARINT=''      # Monte carlo barostat interval in steps
GAMMA_LN=''
DT=''
NSCM=0
NTWX=500
NTPR=50
NTWR=500
CUT='8.0'
JFASTW=''

# ParseAmberOptions
ParseAmberOptions() {
  while [ ! -z "$1" ] ; do
    case "$1" in
      'ntmin'         ) shift ; NTMIN=$1 ;;
      'maxcyc'        ) shift ; MAXCYC=$1 ;;
      'ncyc'          ) shift ; NCYC=$1 ;;
      'restraintmask' ) shift ; RESTRAINTMASK=$1 ;;
      'restraint_wt'  ) shift ; RESTRAINT_WT=$1 ;;
      'irest'         ) shift ; IREST=$1 ;;
      'nstlim'        ) shift ; NSTLIM=$1 ;;
      'ntb'           ) shift ; NTB=$1 ;;
      'cut'           ) shift ; CUT=$1 ;;
      'ntc'           ) shift ; NTC=$1 ; NTF=$1 ;;
      'flexiblewat'   )         NTC=1 ; NTF=1 ; JFASTW='jfastw = 4,' ;;
      'tempi'         ) shift ; TEMPI=$1 ;;
      'tautp'         ) shift ; TAUTP=$1 ;;
      'taup'          ) shift ; TAUP=$1 ;;
      'mcbarint'      ) shift ; MCBARINT=$1 ;;
      'gamma_ln'      ) shift ; GAMMA_LN=$1 ;;
      'dt'            ) shift ; DT=$1 ;;
      'nscm'          ) shift ; NSCM=$1 ;;
      'ntwx'          ) shift ; NTWX=$1 ;;
      'ntpr'          ) shift ; NTPR=$1 ;;
      'ntwr'          ) shift ; NTWR=$1 ;;
      'previousref'   ) REF=$RST ;;
      'heavyrst'      ) shift ; RESTRAINT_WT=$1 ; RESTRAINTMASK="$HEAVYMASK" ;;
      'bbrst'         ) shift ; RESTRAINT_WT=$1 ; RESTRAINTMASK="$BACKBONEMASK" ;;
      'thermo'        ) shift ; THERMOTYPE=$1 ;;
      'baro'          ) shift ; BAROTYPE=$1 ;;
      *               ) echo "Unrecognized option: $1" >> /dev/stderr ; exit 1 ;;
    esac
    shift
  done
  if [ ! -z "$MASTERCUT" ] ; then
    echo "Overriding cut of $CUT with $MASTERCUT"
    CUT=$MASTERCUT
  fi
}

# RestraintLine
RestraintLine() {
  if [ ! -z "$RESTRAINTMASK" ] ; then
    if [ ! -z "$ADDITIONALMASK" ] ; then
      MASKEXP=$RESTRAINTMASK"|"$ADDITIONALMASK
    else
      MASKEXP=$RESTRAINTMASK
    fi
    echo "   ntr = 1, restraintmask = \"$MASKEXP\", restraint_wt = $RESTRAINT_WT," >> $MDIN
  elif [ ! -z "$ADDITIONALMASK" -a ! -z "$RESTRAINT_WT" ] ; then
    echo "   ntr = 1, restraintmask = \"$ADDITIONALMASK\", restraint_wt = $RESTRAINT_WT," >> $MDIN
  elif [ ! -z "$PRODUCTIONMASK" ] ; then
    echo "   ntr = 1, restraintmask = \"$PRODUCTIONMASK\", restraint_wt = $PRODUCTIONWT," >> $MDIN
    if [ ! -z "$PRODUCTIONREF" ] ; then
      if [ ! -f "$PRODUCTIONREF" -a $TEST -eq 0 ] ; then
        echo "Error: production reference $PRODUCTIONREF not found." >> /dev/stderr
        exit 1
      fi
      REF=$PRODUCTIONREF
    fi
  else
    echo "   ntr = 0," >> $MDIN
  fi 
}

# CharmmWater
CharmmWater() {
  if [ $CHARMMWATERFLAG -eq 1 ] ; then
    echo "   WATNAM = 'TIP3', OWTNM = 'OH2'," >> $MDIN
  fi
}

# Minimization input
# CreateMinInput <step> <options>
CreateMinInput() {
  RUN=$1
  MDIN="$RUN".in
  shift
  NTMIN=2
  MAXCYC=1000
  NCYC=10
  NTWX=500
  NTPR=50
  NTWR=500
  RESTRAINTMASK=''
  RESTRAINT_WT=''
  ParseAmberOptions $*
  cat > $MDIN <<EOF
Minimization: $MDIN
 &cntrl
   imin = 1, ntmin = $NTMIN, maxcyc = $MAXCYC, ncyc = $NCYC,
   ntwx = $NTWX, ioutfm = 1, ntxo = 2, ntpr = $NTPR, ntwr = $NTWR, 
   ntc = 1, ntf = 1, ntb = 1, cut = $CUT, $JFASTW
EOF
  CharmmWater
  RestraintLine
cat >> $MDIN <<EOF
 &end
EOF
  RST="$RUN".ncrst
  if [ $OVERWRITE -eq 1 -o ! -f "$RST" ] ; then
    echo "Minimization: $RUN"
    if [ $TEST -eq 0 ] ; then
      $PROG_MIN -O -i $MDIN -p $TOP -c $CRD -ref $REF -o $RUN.out -x $RUN.nc -r $RST -inf $RUN.mdinfo
    else
      echo "TEST: $PROG_MIN -O -i $MDIN -p $TOP -c $CRD -ref $REF -o $RUN.out -x $RUN.nc -r $RST -inf $RUN.mdinfo"
    fi
    if [ $? -ne 0 ] ; then
      echo "Error: Minimization failed: $MDIN" >> /dev/stderr
      exit 1
    fi
  else
    echo "Skipping $RUN"
  fi
  CRD=$RST
}

# Barostat depends on BAROTYPE
# Barostat <file>
Barostat() {
  if [ "$BAROTYPE" = 'berendsen' ] ; then
    echo "   ntp = $NTPFLAG, taup = $TAUP, pres0 = 1.0," >> $1
  elif [ "$BAROTYPE" = 'montecarlo' ] ; then
    echo "   ntp = $NTPFLAG, barostat = 2, pres0 = 1.0, mcbarint = $MCBARINT," >> $1
  else
    echo "Unrecognized BAROTYPE $BAROTYPE"
    exit 1
  fi
}

# Thermostat depends on THERMOTYPE
# Thermostat <file>
Thermostat() {
  if [ "$THERMOTYPE" = 'berendsen' ] ; then
    echo "   ntt = 1, tautp = $TAUTP, temp0 = $TEMP0, tempi = $TEMPI," >> $1
  elif [ "$THERMOTYPE" = 'langevin' ] ; then
    echo "   ntt = 3, gamma_ln = $GAMMA_LN, temp0 = $TEMP0, tempi = $TEMPI," >> $1
  else
    echo "Unrecognized THERMOTYPE $THERMOTYPE"
    exit 1
  fi
}

# MD input
# CreateMdInput <file> <options>
CreateMdInput() {
  RUN=$1
  MDIN="$RUN".in
  shift
  NSTLIM=5000
  IREST=0
  NTB=1
  TAUTP='1.0'
  TAUP='1.0'
  MCBARINT=100
  GAMMA_LN='5'
  DT='0.001'
  NSCM=0
  NTWX=500
  NTPR=50
  NTWR=500
  RESTRAINTMASK=''
  RESTRAINT_WT=''
  TEMPI=$TEMP0
  NTC=2
  NTF=2
  ParseAmberOptions $*
  if [ $IREST -eq 0 ] ; then
    NTX=1
  else
    NTX=5
  fi
  cat > $MDIN <<EOF
MD: $MDIN
 &cntrl
   imin = 0, nstlim = $NSTLIM, dt = $DT, 
   ntx = $NTX, irest = $IREST, ig = -1,
   ntwx = $NTWX, ntwv = -1, ioutfm = 1, ntxo = 2, ntpr = $NTPR, ntwr = $NTWR, 
   iwrap = 0, nscm = $NSCM,
   ntc = $NTC, ntf = $NTF, ntb = $NTB, cut = $CUT, $JFASTW 
EOF
  Thermostat $MDIN
  if [ $NTB -eq 2 ] ; then
    Barostat $MDIN
  fi
  CharmmWater
  RestraintLine
cat >> $MDIN <<EOF
 &end
EOF
  RST="$RUN".ncrst
  if [ $OVERWRITE -eq 1 -o ! -f "$RST" ] ; then
    echo "MD: $RUN"
    if [ $TEST -eq 0 ] ; then
      $PROG_MD -O -i $MDIN -p $TOP -c $CRD -ref $REF -o $RUN.out -x $RUN.nc -r $RST -inf $RUN.mdinfo
    else
      echo "TEST: $PROG_MD -O -i $MDIN -p $TOP -c $CRD -ref $REF -o $RUN.out -x $RUN.nc -r $RST -inf $RUN.mdinfo"
    fi
    if [ $? -ne 0 ] ; then
      echo "Error: MD failed: $MDIN" >> /dev/stderr
      echo "Equilibration failed" > $STATUSFILE
      exit 1
    fi
  else
    echo "Skipping $RUN"
  fi
  CRD=$RST
}

# Standard Equil for explicit solvent
StandardEq() {
  CreateMinInput step1             heavyrst 5.0 $FLEXIBLEWAT
  CreateMdInput  step2 previousref heavyrst 5.0 nstlim 15000 tautp 0.5 $FLEXIBLEWAT
  CreateMinInput step3 previousref heavyrst 2.0 $FLEXIBLEWAT
  CreateMinInput step4 previousref heavyrst 0.1 $FLEXIBLEWAT
  CreateMinInput step5 previousref $FLEXIBLEWAT
  CreateMdInput  step6 previousref heavyrst 1.0 ntb 2 $FLEXIBLEWAT
  CreateMdInput  step7             heavyrst 0.5 ntb 2 irest 1 $FLEXIBLEWAT
  CreateMdInput  step8             bbrst    0.5 ntb 2 irest 1 nstlim 10000 $FLEXIBLEWAT
  CreateMdInput  step9 ntb 2 dt $PRODUCTION_DT irest 1 nscm 1000 $FLEXIBLEWAT
  FinalEq
}

# Standard Equil for explicit solvent; reassign velocities every MD step
NoRestartEq() {
  CreateMinInput step1             heavyrst 5.0 $FLEXIBLEWAT
  CreateMdInput  step2 previousref heavyrst 5.0 nstlim 15000 tautp 0.5 $FLEXIBLEWAT
  CreateMinInput step3 previousref heavyrst 2.0 $FLEXIBLEWAT
  CreateMinInput step4 previousref heavyrst 0.1 $FLEXIBLEWAT
  CreateMinInput step5 previousref $FLEXIBLEWAT
  CreateMdInput  step6 previousref heavyrst 1.0 ntb 2 $FLEXIBLEWAT
  CreateMdInput  step7             heavyrst 0.5 ntb 2 $FLEXIBLEWAT
  CreateMdInput  step8             bbrst    0.5 ntb 2 nstlim 10000 $FLEXIBLEWAT
  CreateMdInput  step9 ntb 2 dt $PRODUCTION_DT nscm 1000
  FinalEq
}

# Final density Equil
FinalEq() {
  #if [ -z "$EVALTYPE" ] ; then
  #  echo "--evaltype not specified; skipping final density equilibration."
  #  return 0
  #fi
  echo "Starting final density equilibration."
  THERMOTYPE=$FINALTHERMO
  BAROTYPE=$FINALBARO
  # Check if we can evaluate with cpptraj/script
  useCpptrajEval=0
  if [ "$EVALTYPE" = 'cpptraj' ] ; then
    nohelpfound=`echo "help evalplateau" | $CPPTRAJ | grep "No help found"`
    if [ -z "$nohelpfound" ] ; then
      echo "Using CPPTRAJ to evaluate density plateau."
      useCpptrajEval=1
    else
      echo "$CPPTRAJ does not support density plateau evaluation."
      echo "Only running one round of final density equilibration."
      EVALTYPE=''
    fi
  elif [ "$EVALTYPE" = 'script' ] ; then
    if [ -z "`which EvalEquilibration.sh`" ] ; then
      echo "EvalEquilibration.sh not found."
      echo "Only running one round of final density equilibration."
      EVALTYPE=''
    fi
  fi
  # Loop over final rounds
  if [ "$RUNTYPE" = 'norestart' ] ; then
    finalIrest=0
  else
    finalIrest=1
  fi
  num=1
  DONE=0
  INPCRD=step9.ncrst
  OUTFILES=''
  while [ $DONE -eq 0 ] ; do
    echo "Final $num"
    CreateMdInput final.$num ntb 2 dt $PRODUCTION_DT nscm 1000 nstlim 500000 ntwx 5000 ntpr 500 ntwr 50000 cut 9.0 irest $finalIrest $FLEXIBLEWAT
    # Decide if we are done. 0 = done, 2 = error, otherwise need more.
    ERR=2
    if [ $TEST -eq 1 ] ; then
      echo "Just testing, no final density eval."
      EVALTYPE=''
    fi
    if [ "$EVALTYPE" = 'script' ] ; then
      EvalEquilibration.sh $OUTFILES > evalEquil.dat
      ERR=$?
    elif [ "$EVALTYPE" = 'cpptraj' ] ; then
      cpptraj -o Eval.out <<EOF
for FILE in final.?.out,final.??.out
  readdata \$FILE name MD
done
evalplateau *[Density] name EQ out Eval.agr resultsout Eval.results
EOF
      if [ $? -eq 0 ] ; then
        awk 'BEGIN{
          resultCol = -1;
        }{
          if (resultCol == -1) {
            for (col = 1; col <= NF; col++) {
              if ($col == "EQ[result]") {
                resultCol = col;
                break;
              }
            }
          } else {
            if ($resultCol == "yes")
              exit 0;
            else if ($resultcol == "err")
              exit 2;
            else
              exit 1;
          }
        }' Eval.results
        ERR=$?
      fi
    else
      echo "Skipping final density evaluation."
      ERR=0
    fi
    if [ $ERR -eq 0 ] ; then
      echo "Complete."
      DONE=1
      if [ $TEST -eq 0 ] ; then
        echo "Equilibration success" > $STATUSFILE
      fi
    elif [ $ERR -eq 2 ] ; then
      echo "Equlibration eval failed"
      echo "Equlibration eval failed" > $STATUSFILE
      break
    else
      INPCRD="final.$num.ncrst"
      ((num++))
      # Safety valve
      if [ $num -gt 20 ] ; then
        echo "More than 20 iterations of final density equil required. Bailing out."
        echo "Too many final iterations" > $STATUSFILE
        DONE=1
        continue
      fi
    fi
  done # END loop over final rounds
  return 0
}

# Determine run type if not already set.
if [ -z "$RUNTYPE" ] ; then
  if [ -z "$INPUT" ] ; then
    # Nothing specified. Default to standard.
    RUNTYPE='standard'
  else
    RUNTYPE='input'
  fi
fi
# ==============================================================================
echo ""
if [ "$RUNTYPE" = 'standard' ] ; then
  echo "Performing standard min/equil"
  StandardEq
elif [ "$RUNTYPE" = 'norestart' ] ; then
  echo "Performing standard min/equil with no restarts"
  NoRestartEq
elif [ "$RUNTYPE" = 'input' ] ; then
  NLINES=0
  # Read input
  echo "Reading input from file: $INPUT"
  while read OPTLINE ; do
    CMD=`echo "$OPTLINE" | awk '{print $1;}'`
    #echo "'$CMD' '$OPTLINE'"
    case "$CMD" in
      'min' ) CMDS[$NLINES]=$CMD ; OPTS[$NLINES]=${OPTLINE#min} ; ((NLINES++)) ;;
      'md'  ) CMDS[$NLINES]=$CMD ; OPTS[$NLINES]=${OPTLINE#md}  ; ((NLINES++)) ;;
      '' | '#' ) continue ;;
      *        ) echo "Unrecognized input: $CMD" >> /dev/stderr ; exit 1 ;;
    esac
  done < $INPUT
  # Execute input 
  for ((i=0; i < $NLINES; i++)) ; do
    #echo "${CMDS[$i]} ${OPTS[$i]}"
    if [ "${CMDS[$i]}" = 'min' ] ; then
      CreateMinInput ${OPTS[$i]}
    else
      CreateMdInput ${OPTS[$i]}
    fi
  done
else
  echo "Error: Unrecognized run type: $RUNTYPE"
  exit 1
fi
exit 0
