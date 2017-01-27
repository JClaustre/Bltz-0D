#!/bin/bash
#############################################
#
# script To compute Boltzmann code without me! :)
#
#############################################
echo "Welcome in the shell script for busy (or lazy!) guys!";

echo
echo "Ready to Start ?";
echo

make clean
make
# **** Incrementation de la pression et pression initiale
let "incr_pres = 100" ; let "pres = 7" ; let "pres_mx=760"
# **** Champ Electrique max (Td)
let "Emx=101" ; let "Emn=1"
let "incr_E = ($Emx-$Emn) / 30"

# **** Taille du tableau E_tab (len-1)
let "lenE = $Emx/$incr_E + 1"

echo "$lenE, $incr_E ... entier?"

# **** Changement de la variable start_a pour ecriture dans le fichier
sed -i "s/start_a=.*/start_a=0/" SRC_F90/evolution.f90
# **** Changement de la variable ETownsd pour l'initialisation
sed -i "s/ETownsd=.*/ETownsd=$Emn/" SRC_F90/evolution.f90

# **** Loop sur les valeurs de la pression 
while [ $pres -le $pres_mx ] ; do

    echo "pression : $pres"
    let "a=0" ; let "pres2 = pres"
    # **** Loop sur les valeurs du champ E
    for val in `seq 1 $lenE`
    do
	if [ $val -eq 1 ] ; then
	    sed -i "s/ETownsd=.*/ETownsd=$Emn/" SRC_F90/evolution.f90
	fi
	make
	./run_BOD
	let "a +=$incr_E"
	sed -i "s/ETownsd=.*/ETownsd=$a/" SRC_F90/evolution.f90
	echo "nouvelle valeur de ETownsd : $a "
	if [ $val -eq 1 ] ; then
	    sed -i "s/start_a=.*/start_a=1/" SRC_F90/evolution.f90
	fi
    done
    let "pres += $incr_pres"

    # **** Recupere la ligne pour le changement de pression
    ligne=$(sed -n "/Torr/=" datFile/input_he)

    if [ $pres -gt $pres_mx ] ; then
	sed -i "$ligne s/$pres2/7/" datFile/input_he
    else
	sed -i "$ligne s/$pres2/$pres/" datFile/input_he
    fi
done
