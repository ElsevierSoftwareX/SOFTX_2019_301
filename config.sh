#!/bin/bash
# config.sh
#-----------------------------------------------------------------------
#
# Change Log:
#
#-----------------------------------------------------------------------
#
date
echo "Build in Commands: $0 $*"
#
# build.sh command line argument processing.
#
argv="$*"
clean_flag=0

gnu=0
gfortran=0
g95=0
mpif90=0
debug=0

if [ $#  ==  0 ]; then
  echo " "
  echo " "
  echo " Usage: bash build.sh [host-machine-type] [parallel] [compiler] "
  echo " "
  echo "       [1] [host-machine-type] is one of the following:         "
  echo "               gnu       Linux; GNU compiler by default"
  echo "               g95       G95 compiler"
  echo "               osx       Mac OSX machines"
  echo " "
  echo "       [2] [parallel] is install flag for parallel version, which must be specified after"
  echo "           the host-machine-type argument.  You may specify any of the following: "
  echo "            0 compiles the serial version of the software                          "
  echo "            1 compiles the parallel version of the software                        "
  echo "       [3] [compiler] is one of the following:                                    "
  echo "            G95  Uses  g95/Linux for gnu (default is gfortran)."
  echo "            gfortran  Uses gfortran/Linux for gnu (default is gfortran)."
  echo "            MPIF90  Relies entirely on mpif90 wrapper for MPI compiling/linking."
  echo " DEBUG/debug Compile with debugging options to compiler                         "
  echo " "
  exit
#
elif [ $#  ==  1 ]; then
  sifm_host=$1
  gfortran=1
  debug=1
elif [ $# -ge 2 ]; then
  sifm_host=$1
  echo "SIFM host computer:    $sifm_host"
  shift
  for opt in $argv ; do 
    echo "$opt"	
    if [ "$opt"  ==  "clean" ]; then
	     clean_flag=1
	fi
    if [ "$opt"  ==  'G95' ] || [ "$opt"  ==  'g95' ]; then
      	g95=1
	fi 		 
    if [ "$opt"  ==  'GFORTRAN' ] || [ "$opt"  ==  'gfortran' ]; then
	     gfortran=1
	fi 	 
    if [ "$opt"  ==  'parallel' ] || [ "$opt"  ==  'PARALLEL' ]; then 
	    mpif90=1
    fi 
  done
fi
#
# ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  ==  == =
#                        ENVIRONMENT SETUP
# the host machine type environment variable chmhost
#
export sifmhost=$sifm_host
#
sifm_root="C:/Users/hkamb/Documents/pubs/paper_softwareX/Revised/sifmv1" 
export root=$sifm_root 

if [ ! -e $sifm_root ]; then 
     echo "root $sifm_root"
	 exit 1
fi

if [ ! -e $root/build ]; then 
     mkdir ${root}/build 
fi
if [ ! -e $root/lib ]; then 
     mkdir ${root}/lib 
fi
if [ ! -e $root/exec ]; then 
     mkdir ${root}/exec 
fi

export build=$root/build/$sifm_host
export exec=$root/exec/$sifm_host 
export lib=$root/lib/$sifm_host
export src=$root/source 

if [ ! -e $build ]; then 
     mkdir $build 
fi 

if [ ! -e $exec ]; then 
     mkdir $exec
fi 

if [ ! -e $lib ]; then 
     mkdir $lib
fi 

# for Makefiles
export LIB=${lib}
export EXEC=${exec} 

# Pseudosizes clean and distclean remove binaries
if [ $clean_flag  ==  1 ]; then
  if [ -f $exec/*.exe ]; then
    echo "Removing $exec/*.exe"
    rm $exec/*.exe
  fi
  if [ -d $lib ]; then
    echo "Removing $lib/*.o"
    rm $lib/*.a $lib/*.o
  fi 
  if [ -d $build ]; then
    echo "Removing $build/*.(o|mod|)"
    rm $build/*.o $build/*.mod
  fi
fi   

#
if [ ! -e $build ]; then 
     mkdir $build 
fi 

if [ ! -e $exec ]; then 
     mkdir $exec
fi 

if [ ! -e $lib ]; then 
     mkdir $lib
fi 
#
if [ "$sifm_host"  ==  "gnu" ]; then
  if [ $g95  ==  1 ]; then
    export GNU_G95=1 	
  else
    export GFORTRAN=YES	 
  fi 
fi 

if [ $debug ==  1 ]; then
     export DEBUG=YES  
fi   

#
#--------------------------------------------------------------#
#                                                              #   
#    up the build directory with makefiles and pref.dat    #
#                                                              #
#--------------------------------------------------------------#
if [ -e $build ]; then
  for srcmk in $root/build/mk/*.mk $root/build/mk/Makefile_*; do 
        echo "$srcmk"
		echo "$build"
        cp $srcmk $build/ 
  done 
  cp $root/source/*.f90 $build
else
  echo " "
  echo " build.sh > Directory $build does not exist."
  echo "              Creating $build ..."
  mkdir $build
  if [ ! -e $root/build/mk ]; then
    echo " build.sh> Directory $root/build/mk does not exist."
    echo "           Installation can not proceed."
    exit 1
  fi 
  cp $root/build/mk/*.mk $build
  cp $root/source/*.f90 $build
fi 


#
# Modify Makefile_template to replace F90 with gfortran
if [ "$sifm_host"  ==  "gnu" ] && [ $gfortran  ==  1 ]; then
   sed -e 's/FF90=F90/FF90 = gfortran/' \
       $build/Makefile_setup_ser > $build/Makefile_setup_$$
   mv $build/Makefile_setup_$$ $build/Makefile_setup_ser_$sifm_host
   if [ $mpif90 == 1 ]; then
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_embd_mpi > $build/Makefile_embd_$$
       mv $build/Makefile_embd_$$ $build/Makefile_embd_mpi_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_symb_mpi > $build/Makefile_symb_$$
       mv $build/Makefile_symb_$$ $build/Makefile_symb_mpi_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_te_mpi > $build/Makefile_te_$$
       mv $build/Makefile_te_$$ $build/Makefile_te_mpi_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_lte_mpi > $build/Makefile_lte_$$
       mv $build/Makefile_lte_$$ $build/Makefile_lte_mpi_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_mi_mpi > $build/Makefile_mi_$$
       mv $build/Makefile_mi_$$ $build/Makefile_mi_mpi_$sifm_host   
   else 
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_embd_ser > $build/Makefile_embd_$$
       mv $build/Makefile_embd_$$ $build/Makefile_embd_ser_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_symb_ser > $build/Makefile_symb_$$
       mv $build/Makefile_symb_$$ $build/Makefile_symb_ser_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_te_ser > $build/Makefile_te_$$
       mv $build/Makefile_te_$$ $build/Makefile_te_ser_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_lte_ser > $build/Makefile_lte_$$
       mv $build/Makefile_lte_$$ $build/Makefile_lte_ser_$sifm_host
       sed -e 's/FF90=F90/FF90 = gfortran/' \
          $build/Makefile_mi_ser > $build/Makefile_mi_$$
       mv $build/Makefile_mi_$$ $build/Makefile_mi_ser_$sifm_host      
   fi 
fi   

echo "Build:  $build"
   
goto="cd $build"
echo "goto:   $goto"
$goto
cdir=`pwd`
echo "current directory:   $cdir"

go_make="make -f Makefile_setup_ser_$sifm_host"
echo "go-make:    $go_make"
$go_make 
  
if [ $mpif90 == 1 ]; then 

  go_make="make -f Makefile_embd_mpi_$sifm_host"
  echo "go-make MPI:    $go_make"
  $go_make 

  go_make="make -f Makefile_symb_mpi_$sifm_host"
  echo "go-make MPI:    $go_make"
  $go_make 

  go_make="make -f Makefile_te_mpi_$sifm_host"
  echo "go-make MPI:    $go_make"
  $go_make 
  go_make="make -f Makefile_lte_mpi_$sifm_host"
  echo "go-make MPI:    $go_make"
  $go_make 
  go_make="make -f Makefile_mi_mpi_$sifm_host"
  echo "go-make MPI:    $go_make"
  $go_make 
  
else

  go_make="make -f Makefile_embd_ser_$sifm_host"
  echo "go-make:    $go_make"
  $go_make 

  go_make="make -f Makefile_symb_ser_$sifm_host"
  echo "go-make:    $go_make"
  $go_make 

  go_make="make -f Makefile_te_ser_$sifm_host"
  echo "go-make:    $go_make"
  $go_make 
  go_make="make -f Makefile_lte_ser_$sifm_host"
  echo "go-make:    $go_make"
  $go_make 
  go_make="make -f Makefile_mi_ser_$sifm_host"
  echo "go-make:    $go_make"
  $go_make 
  
fi 

mv $build/*.exe $exec/.
   
echo "Installation of SIFM Completed"

date 
exit 1
