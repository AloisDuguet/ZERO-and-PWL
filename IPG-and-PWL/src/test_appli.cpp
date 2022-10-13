#include <zero.h>
#include <iostream>
#include <fstream>
#include "read_from_armadillo.cpp"
using namespace std;

int main(int argc, char *argv[]) {
  cout << "starting test_instance" << endl;
  GRBEnv GurobiEnv;
  try {

     // parsing here the arguments argc, argv
     string filename_instance;
	  string filename_output;
     if (argc >= 3) {
        filename_instance = argv[1];
		  filename_output = argv[2];
     } else if (argc == 2) {
		  filename_instance = argv[1];
		  filename_output = filename_instance.substr(0,filename_instance.length()-5)+"outputs/output.txt";
     } else {
		 cout << "at least one argument giving the path to the csv files describing the instance AND 'model' needs to be given.\nExample: path_to_instance/model" << endl;
		 cout << "the 'model' is because all CSV files starts by 'model', or 'test_csv' for the older ones." << endl;
		 cout << "a second argument giving the path and name of the file to write the solution is possible." << endl;
		 return 0;
	  }
     
     arma::vec sizes1 = readIntegerIndexes_CSV(filename_instance+"_sizes1.csv");
     // declaring sizes of matrices and arrays building the MILP model
     Models::IPG::IPGInstance IPG_Instance;
     int                      numBaseVars1 = static_cast<int>(sizes1.at(0));
     int                      numBaseConstraints1 = static_cast<int>(sizes1.at(1));

     int n_players = static_cast<int>(sizes1.at(2));
     int n_markets = static_cast<int>(sizes1.at(3));
     cout << n_players << " players and " << n_markets << " markets" << endl;
     cout << numBaseVars1 << " variables for player 1\n" << numBaseConstraints1 << " constraints for player 1" << endl;


     // Objective terms
     arma::vec    c1(numBaseVars1);
     arma::sp_mat C1(numBaseVars1*(n_players-1), numBaseVars1);
     // Constraints
     arma::sp_mat A1(numBaseConstraints1, numBaseVars1);
     arma::vec    b1(numBaseConstraints1);


     // PLAYER 1
     c1 = readSparseVector_CSV(filename_instance+"_c1.csv",numBaseVars1);
     C1 = readMatrix_CSV(filename_instance+"_C1.csv",numBaseVars1*(n_players-1),numBaseVars1);
     A1 = readMatrix_CSV(filename_instance+"_A1.csv",numBaseConstraints1,numBaseVars1);
     b1 = readVector_CSV(filename_instance+"_b1.csv",numBaseConstraints1);
     // The index of the integer variables
     arma::vec IntegerIndexes1 = readIntegerIndexes_CSV(filename_instance+"_IntegerIndexes1.csv");
     
     // write bounds on variables {0, 1e50} for each variable, then change those of binary variables to {0,1} according to IntegerIndexes1
     VariableBounds VarBounds1 = {};
     for (unsigned int i = 0; i < numBaseVars1; ++i)
        VarBounds1.push_back({0,pow(10.0,50.0)});
     for (int i = 0; i < IntegerIndexes1.size(); i++)
        VarBounds1.at(IntegerIndexes1.at(i)).second = 1;
     cout << "size of VarBounds1 is " << VarBounds1.size() << endl;

     cout << "before creating PlayerOne the different sizes are:\n" << "C1 " << size(C1) << "\nc1 " << size(c1) << "\nA1 " << size(A1) << "\nb1 " << size(b1) << endl;
     cout << "IntegerIndexes1 " << size(IntegerIndexes1) << "\nVarBounds1 " << size(VarBounds1) << endl;
     C1.print("C1:\n");
     IntegerIndexes1.print("IntegerIndexes1:\n");
     
     // Create a parametrized Integer Program
     MathOpt::IP_Param PlayerOne(C1, A1, b1, c1, IntegerIndexes1, VarBounds1, &GurobiEnv);
     // missing constant term in the objective function, add it after optimization for comparison purpose, as well as linear terms in spi




     cout << "PLAYER 2" << endl;
     // PLAYER 2
     arma::vec sizes2 = readIntegerIndexes_CSV(filename_instance+"_sizes2.csv");
     int                      numBaseVars2 = static_cast<int>(sizes2.at(0));
     int                      numBaseConstraints2 = static_cast<int>(sizes2.at(1));
     cout << numBaseVars2 << " variables for player 2\n" << numBaseConstraints2 << " constraints for player 2" << endl;

    // Objective terms
     arma::sp_mat C2(numBaseVars2*(n_players-1), numBaseVars2);
     arma::sp_mat A2(numBaseConstraints2, numBaseVars2);
     // Constraints
     arma::vec    b2(numBaseConstraints2);
     arma::vec    c2(numBaseVars2);

     c2 = readSparseVector_CSV(filename_instance+"_c2.csv",numBaseVars2);
     C2 = readMatrix_CSV(filename_instance+"_C2.csv",numBaseVars2*(n_players-1),numBaseVars2);
     A2 = readMatrix_CSV(filename_instance+"_A2.csv",numBaseConstraints2,numBaseVars2);
     b2 = readVector_CSV(filename_instance+"_b2.csv",numBaseConstraints2);
     // The index of the integer variables
     arma::vec IntegerIndexes2 = readIntegerIndexes_CSV(filename_instance+"_IntegerIndexes2.csv");

     VariableBounds VarBounds2 = {};
     for (unsigned int i = 0; i < numBaseVars2; ++i)
       VarBounds2.push_back({0,pow(10.0,50.0)});
     for (int i = 0; i < IntegerIndexes2.size(); i++)
		 VarBounds2.at(IntegerIndexes2.at(i)).second = 1;


     C2.print("C2:\n");
     IntegerIndexes2.print("IntegerIndexes2:\n");

     cout << "before creating PlayerTwo the different sizes are:\n" << "C2 " << size(C2) << "\nc2 " << size(c2) << "\nA2 " << size(A2) << "\nb2 " << size(b2) << endl;
     cout << "IntegerIndexes2 " << size(IntegerIndexes2) << "\nVarBounds2 " << size(VarBounds2) << endl;

     
     MathOpt::IP_Param PlayerTwo(C2, A2, b2, c2, IntegerIndexes2, VarBounds2, &GurobiEnv);
     cout << "after PLAYER 2" << endl;



     // Add the players to the instance. We also specify a file path to write the instance
     IPG_Instance.addIPParam(PlayerOne, "A_Parametrized_KnapsackProblem1");
     IPG_Instance.addIPParam(PlayerTwo, "A_Parametrized_KnapsackProblem2");
     cout << "Players IP_Param added to model" << endl;
     // Save the instance with the standardize format
     IPG_Instance.save(filename_instance);
     // Create a model from the instance
     Models::IPG::IPG IPG_Model(&GurobiEnv, IPG_Instance);
     cout << "IPG_Model defined" << endl;
     // Select the algorithm to compute a Nash Equilibrium
     IPG_Model.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
     // Extra parameters
     IPG_Model.setDeviationTolerance(3e-4);
     IPG_Model.setNumThreads(4);
     IPG_Model.setLCPAlgorithm(Data::LCP::Algorithms::PATH);
     IPG_Model.setGameObjective(Data::IPG::Objectives::Feasibility);
     IPG_Model.setTimeLimit(60);
     cout << "before finalizing the model" << endl;
     // Lock the model
     IPG_Model.finalize();
     cout << "before starting findNashEq" << endl;
     // Run!
     IPG_Model.findNashEq();
     cout << "optimization finished!" << endl;

     bool solved = IPG_Model.isSolved();
     if(!solved) {
		 Models::IPG::IPGInstance IPG_Instance2;
		 // Add the players to the instance. We also specify a file path to write the instance
		 IPG_Instance2.addIPParam(PlayerOne, "A_Parametrized_KnapsackProblem1");
		 IPG_Instance2.addIPParam(PlayerTwo, "A_Parametrized_KnapsackProblem2");
		 cout << "Players IP_Param added to model" << endl;
		 // Save the instance with the standardize format
		 IPG_Instance2.save(filename_instance);
		 // Create a model from the instance
		 Models::IPG::IPG IPG_Model2(&GurobiEnv, IPG_Instance2);
		 cout << "IPG_Model defined" << endl;
		 // Select the algorithm to compute a Nash Equilibrium
		 IPG_Model2.setAlgorithm(Data::IPG::Algorithms::CutAndPlay);
		 // Extra parameters
		 IPG_Model2.setDeviationTolerance(3e-4);
		 IPG_Model2.setNumThreads(4);
		 IPG_Model2.setLCPAlgorithm(Data::LCP::Algorithms::MIP);
		 IPG_Model2.setGameObjective(Data::IPG::Objectives::Linear);
		 IPG_Model2.setTimeLimit(20);
		 IPG_Model2.finalize();
		 cout << "optimization with MIP algorithm starting because PATH could not optimize it" << endl;
		 IPG_Model2.findNashEq();
		 IPG_Model2.getX().at(0).print("Player 1:");
		 IPG_Model2.getX().at(1).print("\n Player 2:");

		 // compute the objective function of the players
		 cout << "compute the objective values for MIP solution" << endl;
		 arma::vec X1 = IPG_Model2.getX().at(0);
		 arma::vec X2 = IPG_Model2.getX().at(1);
		 double obj1 = 0;
		 double obj2 = 0;
		 for(int i = 0; i < c1.size(); i++) {
		 	obj1 += c1.at(i)*X1.at(i);
		 }
		 for(int i = 0; i < c2.size(); i++) {
			obj2 += c2.at(i)*X2.at(i);
		 }
		 for(int i = 0; i < size(C1)[1]; i++) {
		 	for(int j = 0; j < size(C1)[0]; j++) {
			  obj1 += X2.at(j)*C1.at(j,i)*X1.at(i);
		   }
		 }
		 for(int i = 0; i < size(C2)[1]; i++) {
			for(int j = 0; j < size(C2)[0]; j++) {
			  obj2 += X1.at(j)*C2.at(j,i)*X2.at(i);
			}
		 }
     	 // adding the constant term and the spi term
     	 arma::vec constant_vec1 = readVector_CSV(filename_instance+"_constant1.csv",1);
     	 arma::vec spi_vec1 = readVector_CSV(filename_instance+"_spi_terms1.csv",1);
     	 obj1 += constant_vec1.at(0) + spi_vec1.at(0)*X2.at(n_markets);
     	 arma::vec constant_vec2 = readVector_CSV(filename_instance+"_constant2.csv",1);
     	 arma::vec spi_vec2 = readVector_CSV(filename_instance+"_spi_terms2.csv",1);
     	 obj2 += constant_vec2.at(0) + spi_vec2.at(0)*X1.at(n_markets);
     	 //cout << "the bonus value should be " << constant_vec1.at(0) << " + " << spi_vec1.at(0) << " * " <<X2.at(n_markets) << endl;
     	 cout << "objective value of player 1 is " << obj1 << endl;
     	 cout << "objective value of player 2 is " << obj2 << endl;
     	 cout << "social welfare with constant terms is " << obj1+obj2 << endl;
     	 cout << "social welfare without them is " << IPG_Model2.getSocialWelfare() << endl;

		 // Print the solution in a file
		 cout << "printing the solution in file:\n" << filename_output << endl;
		 std::ofstream savefile(filename_output,std::ios_base::out);
		 savefile << "optimization with MIP algorithm because PATH failed" << endl;
		 savefile << "Player 1:\n" << IPG_Model2.getX().at(0) << "\n" << obj1 << endl;
		 savefile << "Player 2:\n" << IPG_Model2.getX().at(1) << "\n" << obj2 << endl;

     } else {
		// compute the objective function of the players
     	cout << "compute the objective values for PATH solution" << endl;
		IPG_Model.getX().at(0).print("Player 1:");
		IPG_Model.getX().at(1).print("\n Player 2:");
     	arma::vec X1 = IPG_Model.getX().at(0);
     	arma::vec X2 = IPG_Model.getX().at(1);

     	double obj1 = 0;
     	double obj2 = 0;
		for(int i = 0; i < c1.size(); i++) {
		  obj1 += c1.at(i)*X1.at(i);
		}
		for(int i = 0; i < c2.size(); i++) {
		  obj2 += c2.at(i)*X2.at(i);
		}
		for(int i = 0; i < size(C1)[1]; i++) {
		  for(int j = 0; j < size(C1)[0]; j++) {
			 obj1 += X2.at(j)*C1.at(j,i)*X1.at(i);
		  }
		}
		for(int i = 0; i < size(C2)[1]; i++) {
		  for(int j = 0; j < size(C2)[0]; j++) {
			 obj2 += X1.at(j)*C2.at(j,i)*X2.at(i);
		  }
		}

     	// adding the constant term and the spi term
     	arma::vec constant_vec1 = readVector_CSV(filename_instance+"_constant1.csv",1);
     	arma::vec spi_vec1 = readVector_CSV(filename_instance+"_spi_terms1.csv",1);
     	obj1 += constant_vec1.at(0) + spi_vec1.at(0)*X2.at(n_markets);
     	arma::vec constant_vec2 = readVector_CSV(filename_instance+"_constant2.csv",1);
     	arma::vec spi_vec2 = readVector_CSV(filename_instance+"_spi_terms2.csv",1);
     	obj2 += constant_vec2.at(0) + spi_vec2.at(0)*X1.at(n_markets);
     	//cout << "the bonus value should be " << constant_vec1.at(0) << " + " << spi_vec1.at(0) << " * " <<X2.at(n_markets) << endl;
     	cout << "objective value of player 1 is " << obj1 << endl;
     	cout << "objective value of player 2 is " << obj2 << endl;
     	cout << "social welfare with constant terms is " << obj1+obj2 << endl;
     	cout << "social welfare without them is " << IPG_Model.getSocialWelfare() << endl;

		// Print the solution in a file
		cout << "printing the solution in file:\n" << filename_output << endl;
		std::ofstream savefile(filename_output,std::ios_base::out);
		savefile << "optimization with PATH algorithm" << endl;
		savefile << "Player 1:\n" << IPG_Model.getX().at(0) << "\n" << obj1 << endl;
		savefile << "Player 2:\n" << IPG_Model.getX().at(1) << "\n" << obj2 << endl;
     }

  } catch (ZEROException &e) {
     throw ZEROException(e);
  }
}
