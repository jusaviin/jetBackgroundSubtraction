#include "JetBackgroundCard.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)

/*
 * Macro for printing JetBackgroundCard information from a given file
 */ 
 void checkCard(const char *fileName){
  TFile* file = TFile::Open(fileName);
  JetBackgroundCard* card = new JetBackgroundCard(file);
  card->Print();
  file->Close();
  delete card;
 }
