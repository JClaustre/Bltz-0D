#!/bin/bash
#############################################
#
# script To compute Boltzmann code without me! :)
#
#############################################
echo "Welcome in the shell script for busy (or lazy!) guys!";
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`

echo
echo "Ready to Start ? $DATE";
echo
make clean ; make

# **** Incrementation de la pression et pression initiale (bar)
# **** I use "awk" to compute float values
initP=`awk "BEGIN{ print 3 }" `; maxP=`awk "BEGIN{ print 7500 }" ` ; NsimuP=`awk "BEGIN{ print 11 }" `
initE=`awk "BEGIN{ print 5 }" `; maxE=`awk "BEGIN{ print 100 }" `  ; NsimuE=`awk "BEGIN{ print 11 }" `
initT=`awk "BEGIN{ print 8e-11 }" `; maxT=`awk "BEGIN{ print 2e-14 }" `
# **** r is the geometric ratio (raison geometrique in Fr)
rE=`awk "BEGIN{ print ($maxE/$initE)^(1/($NsimuE-1)) }" `
rP=`awk "BEGIN{ print ($maxP/$initP)^(1/($NsimuP-1)) }" `
rT=`awk "BEGIN{ print ($maxT/$initT)^(1/($NsimuP-1)) }" `

# **** Changement de la variable start_a pour ecriture dans le fichier
sed -i "s/start_a=.*/start_a=0/" SRC_F90/evolution.f90
# **** Changement de la variable ETownsd pour l'initialisation
sed -i "s/ETownsd=.*/ETownsd=$initE/" SRC_F90/evolution.f90
# **** Changement de la variable MxDt pour l'initialisation
sed -i "s/MxDt = .*/MxDt = $initT/" SRC_F90/evolution.f90

Torr=$(awk "BEGIN{print $initP }")
Dt=$(awk "BEGIN{print $initT }")

# **** loop on pressure
for valP in `seq 1 $NsimuP`
do

    E=$(awk "BEGIN{print $initE }")
    # **** Loop sur les valeurs du champ E
    for valE in `seq 1 $NsimuE`
    do
	if [ $valE -eq 1 ] ; then
	    sed -i "s/ETownsd=.*/ETownsd=$initE/" SRC_F90/evolution.f90
	fi
	echo "Valeur de ETownsd : $E "
	make
	./run_BOD
	# **** Update E values (Td)
	E=$(awk "BEGIN{print $E*$rE}")
	sed -i "s/ETownsd=.*/ETownsd=$E/" SRC_F90/evolution.f90
	# **** Be sure that the value of start_a has not changed
	sed -i "s/start_a=.*/start_a=1/" SRC_F90/evolution.f90
	# **** Be sure that the value of pressure has not changed
	# **** Read the pressure in file
	ligne=$(sed -n "/Torr/=" datFile/input_he)
	Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
	sed -i "$ligne s/$Torr2/$Torr/" datFile/input_he

    done
    
    if [ $valP -lt $NsimuP ] ; then
	# **** Update the pressure (pres => (torr))
	Torr=$(awk "BEGIN{print $Torr*$rP}")
	Dt=$(awk "BEGIN{print $Dt*$rT}")
	echo "Pressure (Torr): $Torr"

	# **** Recupere la ligne pour le changement de pression
	ligne=$(sed -n "/Torr/=" datFile/input_he)
	# **** Lit la pression dans le fichier
	Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
	# **** Change the pressure in the input file
	sed -i "$ligne s/$Torr2/$Torr/" datFile/input_he
	# **** Change the Max Dt in the main loop file
	sed -i "s/MxDt = .*/MxDt = $Dt/" SRC_F90/evolution.f90

    else
	# **** Recupere la ligne pour le changement de pression
	ligne=$(sed -n "/Torr/=" datFile/input_he)
	# **** Lit la pression dans le fichier
	Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
	# **** Change the pressure in the input file
	sed -i "$ligne s/$Torr2/$initP/" datFile/input_he
    fi
done

echo
echo "Started at : $DATE"
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`
echo "The End at : $DATE"
echo
