#!/bin/bash

# Ce script automatise le téléchargement et l'installation de la bibliothèque Eigen

VERSION="3.4.0"
FILENAME="eigen-$VERSION.zip"
DIRNAME="eigen-$VERSION"
LINKNAME="eigen"
SPLIT="\n##############################"

DO_DOWNLOAD=true
DO_EXTRACT=true
DO_LINK=true

echo -e $SPLIT
echo "Etape 1 : téléchargement."
if [ -f ./$FILENAME ]
then
	echo "Le fichier est déjà téléchargé. Si une erreur se produit, vous devriez essayer de le supprimer pour le retélécharger."
	DO_DOWNLOAD=false
fi

if $DO_DOWNLOAD
then
	echo "Téléchargement de la dernière version"
	URL="https://gitlab.com/libeigen/eigen/-/archive/$VERSION/$FILENAME"

	wget $URL

	if [ ! -f ./$FILENAME ]
	then
		echo "Il y a eu un problème lors du téléchargement : le fichier $FILENAME n'existe pas !"
		exit -2
	fi
fi # Le fichier existe bien

echo -e $SPLIT
echo "Etape 2 : extraction de l'archive"
if [ -d $DIRNAME ]
then
	echo "Le répertoire $DIRNAME existe déjà. Peut-être que la bibliothèque est déjà installée ?"
	QUESTION="Souhaitez-vous supprimer ce répertoire et extraire de nouveau l'archive ? [o/N]"
	if [[ $(read -e -p "${QUESTION}" ; echo $REPLY) == [Oo]* ]]
	then
		rm -rf $DIRNAME
	else
		DO_EXTRACT=false
	fi
fi

if $DO_EXTRACT
then

	unzip $FILENAME
fi

echo -e $SPLIT
echo "Etape 3 : création du lien"
DO_RM_LINK=false
if [ -f $LINKNAME ]
then
	echo "Le fichier $LINKNAME existe déjà et n'est pas un lien. Ceci est un problème. Merci de le corriger à la main."
	exit -3
fi

if [ -e $LINKNAME -a ! -h $LINKNAME ]
then
	echo "Cette erreur ne devrait pas subvenir... Merci de supprimer le fichier/répertoire $LINKNAME."
	exit -4
fi

if [ -h $LINKNAME ]
then
	if [ $(readlink $LINKNAME) != $DIRNAME ]
	then
		echo "Le lien n'est pas correct, nous allons le recréer."
		DO_RM_LINK=true
	else
		echo "Le lien existe déjà et est correct."
		DO_LINK=false
	fi
fi

if $DO_LINK
then
	if $DO_RM_LINK
	then
		#echo "Suppression du lien."
		rm -f $LINKNAME
	fi
	ln -s $DIRNAME $LINKNAME
fi

echo -e $SPLIT
echo "Installation terminée. Vous pouvez maintenant lier vos projets avec cette bibliothèque en vérifiant bien que ce chemin est bien inclus dans les options de compilation."
