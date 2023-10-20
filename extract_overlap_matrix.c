 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <fstream>
 #include <iostream>
 #include <cctype>

 using namespace std;

 // Program extract_overlap_matrix
 // Goal: To extract Overlap Matrix froma Gaussian Output
 // Authors: Julio C. S. Da Silva (julio.silva@iqb.ufal.br
 //          José wilker (josewilker280@gmail.com)
 //          Gabriela Balbino (gabriela.balbino@iqb.ufal.br)

 int main() {
 
 string linha,texto;
 int pos0, pos1, pos;
 
 fstream arq_matriz;
 arq_matriz.open("./C60-1-NO2-pyrene-C60-orbitals.log", ios:: in | ios::out | ios::app); //open file (gaussian)

 while(getline(arq_matriz,linha)){
   if(linha.find("** Overlap ***") != std::string::npos){ //inicio da posição
      pos0=arq_matriz.tellg(); cout << linha << endl; break;
   }
 }
 arq_matriz.seekg(pos0,ios::beg);
 
 while(getline(arq_matriz,linha)){ if(linha.find( "     7     ") != std::string::npos){pos1=arq_matriz.tellg(); break;}}
 arq_matriz.seekg(pos0,ios::beg); //mudando o ponteiro para pos0
 
 while(getline(arq_matriz,linha)){if(arq_matriz.tellg() < pos1){if(linha.find("     2     ")!=  std::string::npos){arq_matriz.seekg(arq_matriz.tellg(),ios::beg);}else{
 linha.erase(0,8); // clearing markers and spaces
 texto=texto+linha+'\n';
 }}else{break;}} // cleaning up unwanted characters
 
 cout<<texto<<endl;
 arq_matriz.close();
 
 return 0;
 }
