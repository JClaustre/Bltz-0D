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
initPr=`awk "BEGIN{ print 0.3 }" ` ; maxPr=`awk "BEGIN{ print 30 }" ` ; NPrsimu=`awk "BEGIN{ print 8 }" `
initPo=`awk "BEGIN{ print 0.002}" `; maxPo=`awk "BEGIN{ print 1.0}" ` ; NPosimu=`awk "BEGIN{ print 8 }" `
initIl=`awk "BEGIN{ print 0.25 }" ` ; maxIl=`awk "BEGIN{ print 2.0}" ` ; NIlsimu=`awk "BEGIN{ print 3 }" `
initAr=`awk "BEGIN{ print 30E-6}" `; maxAr=`awk "BEGIN{ print 5E-6}" `; NArsimu=`awk "BEGIN{ print 8 }" `
#initDt=`awk "BEGIN{ print 3000}" `; maxDt=`awk "BEGIN{ print 3020}" ` ; NDtsimu=`awk "BEGIN{ print 8 }" `

# **** r is the geometric ratio (raison geometrique in Fr)
rPr=`awk "BEGIN{ print ($maxPr/$initPr)^(1/($NPrsimu-1)) }" `
rPo=`awk "BEGIN{ print ($maxPo/$initPo)^(1/($NPosimu-1)) }" `
rAr=`awk "BEGIN{ print ($maxAr/$initAr)^(1/($NArsimu-1)) }" `
#rIl=`awk "BEGIN{ print ($maxIl/$initIl)^(1/($NIlsimu-1)) }" `
#rDt=`awk "BEGIN{ print ($maxDt/$initDt)^(1/($NPrsimu-1)) }" `

# **** Changement de la variable lasr%Is pour l'initialisation
#sed -i "s/lasr%Is = .*/lasr%Is = $initPw/" SRC_F90/evolution.f90

Pr=$(awk "BEGIN{print $initPr }")
Il=$(awk "BEGIN{print $initIl }")
Ar=$(awk "BEGIN{print $initAr }")
#Dt=$(awk "BEGIN{print $initDt }")

# **** Loop sur les valeurs de la puissance laser
for valPr in `seq 1 $NPrsimu`
do
    Po=$(awk "BEGIN{print $initPo }")
    echo "Valeur de Pression  : $Pr"
    # **** Update Dir File name
    ln=$(sed -n "/MEOP/Calcul_Ai_P/=" SRC_F90/param.f90)
    ligne=`awk -F "\/" 'NR==77 { print $5 }' SRC_F90/param.f90`
    int=$(printf '%.1f' "$Pr")
    int2=$(printf '%.1f' "$Il")
    sed -i "$ln s/$ligne/$int\_Torr\_1mW\_$int2\_W/" SRC_F90/param.f90

    #    if [ $valPr -ge 7 ]
    #    then

    for valPo in `seq 1 $NPosimu`
    do
	echo "Valeur de Pression  : $Pr et Critere d'arret : $Ar"
	echo "########## Valeur de Polarisation  : $Po ##########"
	# **** Update Polarization
	ln=$(sed -n "/Tagada/=" SRC_F90/pumping.f90)
	ligne=`awk 'NR==124 { print $3 }' SRC_F90/pumping.f90`
	sed -i "$ln s/$ligne/$Po/" SRC_F90/pumping.f90

	# **** Copie des fichier Steady State vers le repertoire pour le calcul des Ai:
	mkdir ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/
	cp -r ./datFile/MEOP/$int\_Torr\_1mW/* ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/.
	# **** Changement des variable dans le fichier Rstart ************************************************
	# **** Laser On
	Torr=`awk 'NR==18 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	ligne=`awk "BEGIN{ print 18 }" `
	sed -i "$ligne s/$Torr/1/" datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he
	# **** Laser Power
	Torr=`awk 'NR==21 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	ligne=`awk "BEGIN{ print 21 }" `
	sed -i "$ligne s/$Torr/$int2\e+00/" datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he
	# **** Laser Section (m2)
	Torr=`awk 'NR==22 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	ligne=`awk "BEGIN{ print 22 }" `
	sed -i "$ligne s/$Torr/0.001963/" datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he
	# **** Laser Time Start
	Torr=`awk 'NR==26 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	int3=$(printf '%.4f' "$Torr")
	Torr2=`awk 'NR==23 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	ligne=`awk "BEGIN{ print 23 }" `
	sed -i "$ligne s/$Torr2/$int3/" datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he
	# **** end Time of simulation 
	Torr=`awk 'NR==13 { print $1}' ./datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he`
	Torr2=$(awk "BEGIN{print $int3+$Ar}")
	ligne=`awk "BEGIN{ print 13 }" `
	sed -i "$ligne s/$Torr/$Torr2/" datFile/MEOP/Calcul_Ai_P/$int\_Torr\_1mW\_$int2\_W/Rstart/Rs_input_he
	# ****************************************************************************************************
	make clean
	make
	./run_BOD
	# ****************************************************************************************************

	# **** Update the polarization 
	Po=$(awk "BEGIN{print $Po*$rPo}")
    done
    #    fi

    # **** Changement de variables pour l'initialisation       
    # **** Read the pressure in file                          
    ligne=$(sed -n "/Torr/=" datFile/input_he)
    Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
    sed -i "$ligne s/$Torr2/$Pr/" datFile/input_he

    echo "Update Pressure and simulation Stop time" $valPr, $Pr, $Ar
    echo 
    # **** Update the pressure (pres => (torr)) 
    Pr=$(awk "BEGIN{print $Pr*$rPr}")
    Ar=$(awk "BEGIN{print $Ar*$rAr}")
    # **** Update the Laser Intensity 
    #Il=$(awk "BEGIN{print $Il*$rIl}")
done
 

echo
echo "Started at : $DATE"
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`
echo "The End at : $DATE"
echo
