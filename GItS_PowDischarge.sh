#!/bin/bash
#############################################
#
# script To compute Boltzmann code without me! :)
#
#############################################
echo "Welcome in the shell script for busy (or lazy!) guys!";
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`

echo
echo "Ready to Start to change the polarization ? $DATE";
echo

# **** Incrementation de la pression et pression initiale (bar)
# **** I use "awk" to compute float values
initPr=`awk "BEGIN{ print 0.3 }" `; maxPr=`awk "BEGIN{ print 8.0 }" `; NPrsimu=`awk "BEGIN{ print 6 }" `
initPw=`awk "BEGIN{ print 0.2 }" `; maxPw=`awk "BEGIN{ print 10.0}" `; NPwsimu=`awk "BEGIN{ print 6}" `
initP2=`awk "BEGIN{ print 2e-3}" `; maxP2=`awk "BEGIN{ print 1.0}"  `; NP2simu=`awk "BEGIN{ print 6}" `
initAr=`awk "BEGIN{ print 1000.1}" `; maxAr=`awk "BEGIN{ print 300.1}" `; NArsimu=`awk "BEGIN{ print 5 }" `
NPow=$(awk "BEGIN{print 5 }")

# **** r is the geometric ratio (raison geometrique in Fr)
rPr=`awk "BEGIN{ print ($maxPr/$initPr)^(1/($NPrsimu-1)) }" `
rPw=`awk "BEGIN{ print ($maxPw/$initPw)^(1/($NPwsimu-1)) }" `
rP2=`awk "BEGIN{ print ($maxP2/$initP2)^(1/($NP2simu-1)) }" `
rAr=`awk "BEGIN{ print ($maxAr/$initAr)^(1/($NArsimu-1)) }" `

Pr=$(awk "BEGIN{print $initPr }")

# **** Define power at the beginning of the simulation
Pw=$(awk "BEGIN{print $initPw }")
P2=$(awk "BEGIN{print $initP2 }")
    
# **** Loop sur les valeurs de la puissance laser
for valPr in `seq 1 $NPrsimu`
do
    echo "**** Valeur de Pression (Torr)      : $Pr"
    echo "**** Valeur de Puissances (Min,Max) : $P2 | $Pw"
    
    # **** Power of Discharge for every pressure
    rPow=`awk "BEGIN{ print ($Pw/$P2)^(1/($NPow-1)) }" `
    Pow=$(awk "BEGIN{print $P2 }")

    Ar=$(awk "BEGIN{print $initAr }")

    if [ $valPr -eq 4 ]
    then

	for valPw in `seq 1 $NPow`
	do

	    if [ $valPw -eq 2 ]
	    then

		# **** Changement de variables pour l'initialisation       
		# **** Read the pressure in file                          
		ligne=$(sed -n "/Torr/=" datFile/input_he)
		Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
		sed -i "$ligne s/$Torr2/$Pr/" datFile/input_he

		echo "#*#*#*#*# Valeur de Puissance (W/cm3)  : $Pow "
		echo "#*#*#*#*# Valeur Arret (sec)           : $Ar  "

		# **** Update Dir File name
		ln=$(sed -n "/MEOP/=" SRC_F90/param.f90)
		ligne=`awk -F "\/" 'NR==77 { print $4 }' SRC_F90/param.f90`
		int=$(printf '%.1f' "$Pr")
		int2=$(printf '%.4f' "$Pow")
		sed -i "$ln s/$ligne/$int\_Torr\_$int2\_W/" SRC_F90/param.f90
    		echo $ligne, $int2

		# **** Changement des variable dans le fichier Rstart ************************************************
		# **** Power discharge
		Torr=`awk 'NR==12 { print $1}' ./datFile/input_he`
		int3=$(printf '%.7f' "$Pow")
		ligne=`awk "BEGIN{ print 12 }" `
		sed -i "$ligne s/$Torr/$int3/" datFile/input_he
		# **** End Time of simulation 
		Torr=`awk 'NR==13 { print $1}' ./datFile/input_he`
		Torr2=$(awk "BEGIN{print $Ar}")
		ligne=`awk "BEGIN{ print 13 }" `
		sed -i "$ligne s/$Torr/$Torr2/" datFile/input_he
		# ****************************************************************************************************
		make clean
		make
		./run_BOD
		# ****************************************************************************************************
	    fi
	    # **** Update the polarization 
	    Pow=$(awk "BEGIN{print $Pow*$rPow}")
	    # **** Update the time of end of simulation 
	    Ar=$(awk "BEGIN{print $Ar*$rAr}")

	done
    fi
    # **** Update the pressure (pres => (torr)) 
    Pr=$(awk "BEGIN{print $Pr*$rPr}")

    # **** Changement de variables pour l'initialisation       
    # **** Read the pressure in file                          
    ligne=$(sed -n "/Torr/=" datFile/input_he)
    Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
    sed -i "$ligne s/$Torr2/$Pr/" datFile/input_he

    echo 
    # **** Update the power (Min(P2) & Max(Pw)) 
    Pw=$(awk "BEGIN{print $Pw*$rPw}")
    P2=$(awk "BEGIN{print $P2*$rP2}")
done
 

echo
echo "Started at : $DATE"
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`
echo "The End at : $DATE"
echo
