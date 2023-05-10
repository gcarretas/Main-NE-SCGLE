#include <stdio.h>
#include <math.h>
#include "./structures/structures.h"
#include "./dynamics/dynamics.h"
#include "./math/math_aux.h"
#include <time.h>
// libraries to handle directories
#include <stdbool.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

// check if some directoriy exist
bool exists(const char *fname)
{
  //printf("into exist()\n");
  FILE *file;
  if ((file = fopen(fname, "r")))
  {
    fclose(file);
    return true;
  }
  return false;
}

void num2char(float x, char * folder){
  char buffer[10];
  sprintf(buffer, "%f", x);
  for(int i=0; i<strlen(buffer); i ++){
    if (buffer[i] == '.') buffer[i] = 'p';
  }
  strcat(folder, buffer);
}

int main(void) {
  
  clock_t start,end;
  start = clock();
  printf("Program started\n");
  // phisical parameters
  //float lambda = 1.5;
  float phi=0.07;
  //float phi; // Equilibrium
  //float phi_eq[]={0.582};
  float Tf; 
  //float Temperatures[] = {0.0195,0.019,0.0188,0.0186,0.0185,0.0184,0.0182,0.0181};
  //float Temperatures[] = {0.01808,0.01806,0.01804,0.01802,0.018,0.01798,0.01796,0.01794,0.01792,0.0179};  
  //float Temperatures[] = {0.016,0.0155,0.015};
  //float phif;
  //float phi_f[]={0.58,0.62};
  // calculate size in bytes
  //int arraySize = sizeof(Temperatures);
  //int floatSize = sizeof(Temperatures[0]);
  //int arraySize = sizeof(phi_f); //NE-Hard Sphere
  //int floatSize = sizeof(phi_f[0]); //NE-Hard Sphere
  //int arraySize = sizeof(phi_eq);
  //int floatSize = sizeof(phi_eq[0]);
  // length

  //int length = arraySize / floatSize;

  for (int i = 1; i<5; i++){
    //for (int i = 0; i< length; i++){
    //phi = phi_eq[i]; // Equilibrium - HS
    //Tf = Temperatures[i];
    Tf = 0.0195 + i*0.0001;
    //phif = phi_f[i]; // NE - Hard Sphere
    // creates storage folder
    // creates a folder to equilibrium
    //char folder[120] = "eq"; // This line makes a new folder for the process in equilibrium
    //if (!exists(folder)) mkdir(folder);
    // lambda variable
    //char folder[120] = "lambda"; //non-equilibrium
    char folder[120] = "SALR_4"; //non-equilibrium
    //char folder[120] = "NE-HS"; //non-equilibrium
    //char folder[120] = "HS"; //equilibrium
    //char folder[120] = "NE-WCA"; //non-equilibrium
    //strcat(folder, "/lambda"); // equilibrium
    //num2char(lambda, folder);
    //num2char(H_S, folder);
    if (!exists(folder)) mkdir(folder, S_IRWXU); /* This works on ubuntu*/
    //if (!exists(folder)) mkdir(folder); /* This works on windows*/
    // volume fraction variable
    strcat(folder, "/phi_f");
    num2char(phi, folder);
    if (!exists(folder)) mkdir(folder, S_IRWXU); /* This works on ubuntu*/
    //if (!exists(folder)) mkdir(folder); /* This works on windows*/
    // phi_f
    //strcat(folder, "/phi_f");
    //num2char(phif, folder);
    // final temperature
    strcat(folder, "/Tf");
    num2char(Tf, folder);
    if (!exists(folder)) mkdir(folder, S_IRWXU);  /*This works on ubuntu*/
    //if (!exists(folder)) mkdir(folder); /* This works on windows*/
    strcat(folder, "/");
    //printf("%s\n", folder);
    
    /*PARÁMETROS HSDY ATERMAL*/
  /*A: K2 repulsiva fija*/
  /*double zA=10.0; double zR=0.5; double A=0.5;*/

  /*PARÁMETROS HSDY*/
  /*A: razón entre K1 y K2*/
  double zA=10.0; double zR=4; double A=0.5;

    //liquid_params lp = liquid_params_ini_phi(phi,3,1.0); // equilibrium
    liquid_params initial_lp = liquid_params_ini_phi(phi,3,0.0);//non-equilibrium process
    liquid_params final_lp   = liquid_params_ini_phi(phi,3,3.0);//non-equilibrium process
    //lp.up[0] = lambda;
    //initial_lp.Tem = 1;
    //final_lp.Tem = Tf;
    //dynamics_mono_spherical_standard_defined_structures(lp, "HS", "VW", folder ); //equilibrium process
    //final_lp.up[0] = lambda;
    //initial_lp.Tem = 0.0001;
    final_lp.Tem = Tf;
    /*PARÁMETROS FINALES HSDY*/
    final_lp.up[0]=Tf/A; final_lp.up[1]=zA; final_lp.up[2]=zR;

    /*PARÁMETROS FINALES HSDY ATERMAL*/
    /*final_lp.up[0]=1/A; final_lp.up[1]=zA; final_lp.up[2]=zR;*/

    //instant_change_dynamics_spherical_mono_standard_defined_structures( lp, final_lp, "HS", "HSSW", "VW", "SH", "./data/tests/" );
    //instant_change_dynamics_spherical_mono_standard_defined_structures( initial_lp, final_lp, "HS", "HSSW", "VW", "SH", folder); // non-equilibrium process Square-Well Hard-Sphere
    //instant_change_dynamics_spherical_mono_standard_defined_structures( initial_lp, final_lp, "HS", "HS","VW", "VW",folder); // non-equilibrium process Hard-Sphere
    //instant_change_dynamics_spherical_mono_standard_defined_structures( initial_lp, final_lp, "WCA", "WCA","BLIP", "BLIP",folder); // non-equilibrium process Soft-Sphere
    instant_change_dynamics_spherical_mono_standard_defined_structures( initial_lp, final_lp, "HS", "HSDBLEYUK","VW", "SH",folder); // non-equilibrium process SALR
    end=clock();
    double time_taken = ((double)(end-start))/ CLOCKS_PER_SEC;
    printf( "%s%f%s\n", "Program ended, time taken : ", time_taken, " sec");
  }
  return 0;
}


//6120
//eric garcia