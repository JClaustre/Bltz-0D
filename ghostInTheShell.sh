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
# **** Using of "awk" to compute float values
pres=`awk "BEGIN{ print 2 }" ` ; raison=`awk "BEGIN{ print 200^0.1 }" `
# **** Incrementation du Champ Electrique max (Td)
Emx=101 ; Emn=5 ; incr_E=5


# **** Taille du tableau E_tab (len-1)
lenE=$(awk "BEGIN{print int($Emx/$incr_E+1)}")

# **** Changement de la variable start_a pour ecriture dans le fichier
sed -i "s/start_a=.*/start_a=0/" SRC_F90/evolution.f90
# **** Changement de la variable ETownsd pour l'initialisation
sed -i "s/ETownsd=.*/ETownsd=$Emn/" SRC_F90/evolution.f90

# **** loop on pressure
for valP in `seq 1 11`
do
    Torr=$(awk "BEGIN{print $pres*750.061}")
    E=5
    # **** Loop sur les valeurs du champ E
    for valE in `seq 1 $lenE`
    do
	if [ $valE -eq 1 ] ; then
	    sed -i "s/ETownsd=.*/ETownsd=$Emn/" SRC_F90/evolution.f90
	fi
	echo "Valeur de ETownsd : $E "
	make
	./run_BOD
	# **** Update E values (Td)
	E=$(awk "BEGIN{print $E+$incr_E}")
	sed -i "s/ETownsd=.*/ETownsd=$E/" SRC_F90/evolution.f90
	# **** Be sure that the value of start_a has not changed
	sed -i "s/start_a=.*/start_a=1/" SRC_F90/evolution.f90
	# **** Be sure that the value of pressure has not changed
	# **** Lit la pression dans le fichier
	ligne=$(sed -n "/Torr/=" datFile/input_he)
	Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
	sed -i "$ligne s/$Torr2/$Torr/" datFile/input_he

    done
    # **** Update of pressure (pres => (bar) ; Torr => (torr))
    pres=$(awk "BEGIN{print $pres/$raison}")
    Torr=$(awk "BEGIN{print $pres*750.061}")
    echo "Pressure (Torr): $Torr"

    # **** Recupere la ligne pour le changement de pression
    ligne=$(sed -n "/Torr/=" datFile/input_he)
    # **** Lit la pression dans le fichier
    Torr2=`awk 'NR==6 { print $1}' ./datFile/input_he`
    # **** Change the pressure in the input file
    sed -i "$ligne s/$Torr2/$Torr/" datFile/input_he

done



echo
echo "Started at : $DATE"
DATE=`date +%Y-%m-%d:%Hh:%Mmin:%Ss`
echo "The End at : $DATE";
echo
