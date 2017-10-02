#!/bin/bash
#############################################
#
# script To compute Boltzmann code without me! :)
#
#############################################
echo "Welcome in the shell script for busy (or lazy!) guys!";
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`

echo
echo "Ready to Start to change the power of the Laser ? $DATE";
echo

# **** Incrementation de la pression et pression initiale (bar)
# **** I use "awk" to compute float values
initPr=`awk "BEGIN{ print 0.3 }" `; maxPr=`awk "BEGIN{ print 30 }" `  ; NPrsimu=`awk "BEGIN{ print 8 }" `
initIl=`awk "BEGIN{ print 1.0 }" `; maxIl=`awk "BEGIN{ print 10 }" `  ; NIlsimu=`awk "BEGIN{ print 2 }" `
initDt=`awk "BEGIN{ print 3000}" `; maxDt=`awk "BEGIN{ print 3020}" ` ; NDtsimu=`awk "BEGIN{ print 8 }" `
# **** r is the geometric ratio (raison geometrique in Fr)
rPr=`awk "BEGIN{ print ($maxPr/$initPr)^(1/($NPrsimu-1)) }" `
rIl=`awk "BEGIN{ print ($maxIl/$initIl)^(1/($NIlsimu-1)) }" `
rDt=`awk "BEGIN{ print ($maxDt/$initDt)^(1/($NPrsimu-1)) }" `

# **** Changement de la variable lasr%Is pour l'initialisation
#sed -i "s/lasr%Is = .*/lasr%Is = $initPw/" SRC_F90/evolution.f90

Pr=$(awk "BEGIN{print $initPr }")
Dt=$(awk "BEGIN{print $initDt }")

# **** Loop sur les valeurs de la puissance laser
for valPr in `seq 1 $NPrsimu`
do
    echo "Valeur de Pression  : $Pr "
    # **** Update Dir File name
    ln=$(sed -n "/MEOP/=" SRC_F90/param.f90)
    ligne=`awk -F "\/" 'NR==77 { print $4 }' SRC_F90/param.f90`
    int=$(printf '%.1f' "$Pr")
    sed -i "$ln s/$ligne/$int\_Torr\_0.1W/" SRC_F90/param.f90

    # **** Changement de variables pour l'initialisation       
    # **** Read the pressure in file                          
    ligne=$(sed -n "/Torr/=" datFile/input_he)
    Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
    sed -i "$ligne s/$Torr2/$Pr/" datFile/input_he
    # **** Change : restart or not!!!!                        
    ligne=$(sed -n "/Restart/=" datFile/input_he)
    Torr2=`awk 'NR==2 { print $1}' ./datFile/input_he`
    sed -i "$ligne s/$Torr2/0/" datFile/input_he
    # **** Change Absorbed Power!!!!                        
    ligne=$(sed -n "/Absorbed/=" datFile/input_he)
    Torr2=`awk 'NR==12 { print $1}' ./datFile/input_he`
    sed -i "$ligne s/$Torr2/0.100000/" datFile/input_he
    # **** Read the Simulation Time in file                          
    #ligne=$(sed -n "/mic-sec/=" datFile/input_he)
    #Torr2=`awk 'NR==13 { print $1}' ./datFile/input_he`
    #sed -i "$ligne s/$Torr2/$Dt/" datFile/input_he

    make clean
    make
    ./run_BOD

    echo "Update Pressure and simulation time" $valPr, $Pr, $Dt
    echo 
    # **** Update the pressure (pres => (torr)) 
    Pr=$(awk "BEGIN{print $Pr*$rPr}")
    # **** Update the Time of simulation 
    #Dt=$(awk "BEGIN{print $Dt*$rDt}")
done
 

echo
echo "Started at : $DATE"
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`
echo "The End at : $DATE"
echo
