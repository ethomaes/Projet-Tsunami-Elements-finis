# include "tsunami.h"

int main(void)
{   

        char *meshName = "../data/PacificFine.txt";
        char *resultBaseName = "output/tsunamiFine";
        
                       
        tsunamiCompute(2.0,100,25,meshName,resultBaseName);
        tsunamiAnimate(2.0,10000,25,meshName,resultBaseName);
         
        exit(0);     
}