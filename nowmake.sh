wget https://twiki.cern.ch/twiki/pub/Main/Geant4Tutorial/FindROOT.cmake
source ~/setupG.sh
cmake -G"Eclipse CDT4 - Unix Makefiles" -DGeant4_DIR=$G4LIB . 
make


