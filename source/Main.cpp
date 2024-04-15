#include <LRAlgo.h>

#include "Data.h"
#include "SAA.h"

/// MAIN PROGRAM ///
int main() {
  srand(0);

  // I : the number of suppliers
  // J: the number of retailers
  // K: the number of cross-docks
  // P: the number of products
  // S: the number of samples
  OCOData Data(5, 30, 10, 10, 100);

  cout << "I = " << Data.getNumSuppliers() << endl;
  cout << "J = " << Data.getNumRetailers() << endl;
  cout << "K = " << Data.getNumCrossDocks() << endl;
  cout << "P = " << Data.getNumProducts() << endl;
  cout << "S = " << Data.getNumSamples() << endl;

  // benchmark: SAA
  SAAModel SAA(Data);

  // input:
  //  data			: Data
  //  algo_type	: 1(directly solve with solver) / 2(solve with OCO)
  // LRAlgo LR1(Data, 1);
  LRAlgo LR2(Data, 2);

  return 0;
}