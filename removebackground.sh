for file in subjects/*; do
newfile=`echo "${file/.jpg/.png}" |cut -d'/' -f2`
#echo "./segment $file subjects/nobackground/$newfile"
./segment $file subjects/nobackground/$newfile
done
